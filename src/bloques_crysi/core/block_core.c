#include "block_core.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "serial_win.h"
#include <time.h>
#ifdef _WIN32
#include <mmsystem.h>   /* timeBeginPeriod — enlazar con -lwinmm */
#pragma comment(lib, "winmm.lib")
#endif

#define SQ3 1.7320508075688772
#define SQ3_2 0.8660254037844386
#define CV_THIRD (1.0/3.0)
#define TWO_THIRD (2.0/3.0)

static double mna_G[128][128];
static double mna_Cmat[128][128];
static double mna_W[128][128];
static double mna_Ap[128][128];
static double mna_b[128];
static double mna_b_diode[128];

/* ---------------- transform helpers (convención: d alineado al fasor a θ=0) ---------------- */
static void f_clarke(double va, double vb, double vc, double *al, double *be) {
  *al = TWO_THIRD * (va - 0.5 * vb - 0.5 * vc);
  *be = CV_THIRD * (vb - vc) * SQ3;
}
static void f_inv_clarke(double al, double be, double *va, double *vb, double *vc) {
  *va = al;
  *vb = -0.5 * al + SQ3_2 * be;
  *vc = -0.5 * al - SQ3_2 * be;
}
static void f_park(double al, double be, double th, double *d, double *q) {
  double c = cos(th), s = sin(th);
  *d = al * c + be * s;
  *q = -al * s + be * c;
}
static void f_inv_park(double d, double q, double th, double *al, double *be) {
  double c = cos(th), s = sin(th);
  *al = d * c - q * s;
  *be = d * s + q * c;
}
/* ---------------- módulos de potencia ----------------
   Convenciones:
   - DC-DC promediado: in [vin, d], estado [iL, vC], param [L, C, R],
     salida [vC, iL] (vout, iL).
       Buck:     L iL' = d*vin - vC ;  C vC' = iL - vC/R
       Boost:    L iL' = vin - (1-d)*vC ; C vC' = (1-d)*iL - vC/R
       BuckBoost:L iL' = d*vin - (1-d)*vC ; C vC' = (1-d)*iL - vC/R
   - Rectificador trifasico (AC-DC): in [va,vb,vc], estado [vC],
     param [C, R, Rint], salida [vdc, idc].
     vrec = max(ab, bc, ca). Carga con Rint (commutacion ideal), descarga R.
   - Inversor (DC-AC): in [vdc], param [f_out, fsw, m_start, m_end, t_ramp,
     Lf, Cf, R, conmutada], estado [iLa,iLb,iLc, vCa,vCb,vCc],
     salida [vCa,vCb,vCc, iLa,iLb,iLc].
     Rampa m(t) = m_start + (m_end-m_start)*min(t/t_ramp,1) (sube o baja).
     conmutada=0: v_polek = m(t)*vdc/2 * sin(w*t - 2k*pi/3)  (promedio)
     conmutada=1: v_polek = +/-vdc/2 (SPWM bipolar, igual que el 1f):
                  vdc/2 si m(t)*sin(...) > carr, -vdc/2 si no.
     carr: triangular -1..1 de periodo 1/fsw.
     Filtro por fase: L iL' = v_pole - vC ; C vC' = iL - vC/R
*/
static void pot_dcdc_deriv(BloqueC *bl, ModeloC *m, const double *x, double *dx) {
  double L = bl->param[0], C = bl->param[1], R = bl->param[2];
  double vin = m->sig[bl->in_idx[0]], d = m->sig[bl->in_idx[1]];
  double iL = x[0], vC = x[1];
  double vpol, iC;
  if (bl->op == OP_POT_BUCK) {
    vpol = d * vin;
    iC = iL;
    dx[0] = (vpol - vC) / L;
  } else if (bl->op == OP_POT_BOOST) {
    vpol = (1.0 - d) * vC;
    iC = (1.0 - d) * iL;
    dx[0] = (vin - vpol) / L;
  } else { /* buck-boost */
    vpol = (1.0 - d) * vC;
    iC = (1.0 - d) * iL;
    dx[0] = (d * vin - vpol) / L;
  }
  dx[1] = (iC - vC / R) / C;
}
static void pot_dcdc_out(BloqueC *bl, ModeloC *m, const double *x) {
  m->sig[bl->out_idx[0]] = x[1]; /* vout = vC */
  m->sig[bl->out_idx[1]] = x[0]; /* iL */
}
static void pot_rect_deriv(BloqueC *bl, ModeloC *m, const double *x, double *dx) {
  double C = bl->param[0], R = bl->param[1], Rint = bl->param[2];
  double va = m->sig[bl->in_idx[0]], vb = m->sig[bl->in_idx[1]], vc = m->sig[bl->in_idx[2]];
  double vmax = fmax(fmax(va, vb), vc), vmin = fmin(fmin(va, vb), vc);
  double vrec = vmax - vmin; /* envolvente real del puente de 6 pulsos */
  double vC = x[0];
  double ich = (vrec - vC) / Rint;
  if (ich < 0.0) ich = 0.0; /* diodos bloquean la descarga por la entrada */
  dx[0] = (ich - vC / R) / C;
}
static void pot_rect_out(BloqueC *bl, ModeloC *m, const double *x) {
  double Rint = bl->param[2];
  double va = m->sig[bl->in_idx[0]], vb = m->sig[bl->in_idx[1]], vc = m->sig[bl->in_idx[2]];
  double vmax = fmax(fmax(va, vb), vc), vmin = fmin(fmin(va, vb), vc);
  double vrec = vmax - vmin;
  double ich = (vrec - x[0]) / Rint;
  if (ich < 0.0) ich = 0.0;
  m->sig[bl->out_idx[0]] = x[0]; /* vdc */
  m->sig[bl->out_idx[1]] = ich;  /* idc */
}
/* v_pole de la fase k del inversor (0,1,2) con los parametros dados */
static double inv_vpole(BloqueC *bl, ModeloC *m, int k, double vdc) {
  double f = bl->param[0], fsw = bl->param[1];
  double m0 = bl->param[2], m1 = bl->param[3], tramp = bl->param[4];
  int conmutada = (int)bl->param[8];
  double mt;
  if (tramp > 0.0) mt = m0 + (m1 - m0) * (m->t / tramp);
  else mt = m1;
  if (tramp > 0.0 && m->t >= tramp) mt = m1; /* fin de rampa (sube o baja) */
  if (mt < 0.0) mt = 0.0;
  double w = 2.0 * M_PI * f;
  double th = w * m->t - (2.0 * M_PI / 3.0) * k; /* potencia.py: b -2pi/3, c -4pi/3 */
  double ref = mt * sin(th);
  if (!conmutada) return 0.5 * vdc * ref;
  double Tsw = 1.0 / fsw;
  double frac = fmod(m->t, Tsw) / Tsw; /* [0,1) */
  double carr = (frac < 0.5) ? (-1.0 + 4.0 * frac) : (3.0 - 4.0 * frac);
  /* SPWM bipolar +/-vdc/2 (igual que el 1f): sin componente DC en cada fase
     (el modo promedio 0.5*vdc*ref ya es bipolar; el modo conmutado debe ser
     consistente o el filtro LC referido a 0 V integra una DC de vdc/2) */
  return (ref > carr) ? (0.5 * vdc) : (-0.5 * vdc);
}
static void pot_inv_deriv(BloqueC *bl, ModeloC *m, const double *x, double *dx) {
  double Lf = bl->param[5], Cf = bl->param[6], R = bl->param[7];
  double vdc = m->sig[bl->in_idx[0]];
  for (int k = 0; k < 3; k++) {
    double vp = inv_vpole(bl, m, k, vdc);
    dx[k] = (vp - x[3 + k]) / Lf;          /* diL/dt */
    dx[3 + k] = (x[k] - x[3 + k] / R) / Cf; /* dvC/dt */
  }
}
static void pot_inv_out(BloqueC *bl, ModeloC *m, const double *x) {
  for (int k = 0; k < 3; k++) {
    m->sig[bl->out_idx[k]] = x[3 + k]; /* vCa, vCb, vCc */
    m->sig[bl->out_idx[3 + k]] = x[k]; /* iLa, iLb, iLc */
  }
}
/* ---------------- inversor monofasico (DC-AC) ----------------
   in [vdc], param [f_out, fsw, m_start, m_end, t_ramp, Lf, Cf, R, conmutada],
   estado [iL, vC], salida [vC, iL].
   v_pole = 0.5*vdc*m(t)*sin(w t) (promedio) o SPWM bipolar +/-vdc/2
   (conmutada=1, portadora triangular [-1,1]): promedio identico.
   L iL' = v_pole - vC ; C vC' = iL - vC/R
*/
static double inv1f_vpole(BloqueC *bl, ModeloC *m, double vdc) {
  double f = bl->param[0], fsw = bl->param[1];
  double m0 = bl->param[2], m1 = bl->param[3], tramp = bl->param[4];
  int conmutada = (int)bl->param[8];
  double mt;
  if (tramp > 0.0) mt = m0 + (m1 - m0) * (m->t / tramp);
  else mt = m1;
  if (tramp > 0.0 && m->t >= tramp) mt = m1; /* fin de rampa (sube o baja) */
  if (mt < 0.0) mt = 0.0;
  double ref = mt * sin(2.0 * M_PI * f * m->t);
  if (!conmutada) return 0.5 * vdc * ref;
  double Tsw = 1.0 / fsw;
  double frac = fmod(m->t, Tsw) / Tsw; /* [0,1) */
  double carr = (frac < 0.5) ? (-1.0 + 4.0 * frac) : (3.0 - 4.0 * frac);
  return (ref > carr) ? (0.5 * vdc) : (-0.5 * vdc);
}
static void pot_inv1f_deriv(BloqueC *bl, ModeloC *m, const double *x, double *dx) {
  double Lf = bl->param[5], Cf = bl->param[6], R = bl->param[7];
  double vp = inv1f_vpole(bl, m, m->sig[bl->in_idx[0]]);
  dx[0] = (vp - x[1]) / Lf;       /* diL/dt */
  dx[1] = (x[0] - x[1] / R) / Cf; /* dvC/dt */
}
static void pot_inv1f_out(BloqueC *bl, ModeloC *m, const double *x) {
  m->sig[bl->out_idx[0]] = x[1]; /* vC */
  m->sig[bl->out_idx[1]] = x[0]; /* iL */
}
/* ---------------- carga RL trifasica (estrella sin neutro) ----------------
   in [va, vb, vc], param [R, L], estado [ia, ib], salida [ia, ib, ic].
   vn = (va+vb+vc)/3 ;  L i' = v - vn - R i ;  ic = -(ia+ib)
*/
static void carga_rl_deriv(BloqueC *bl, ModeloC *m, const double *x, double *dx) {
  double R = bl->param[0], L = bl->param[1];
  double va = m->sig[bl->in_idx[0]], vb = m->sig[bl->in_idx[1]],
         vc = m->sig[bl->in_idx[2]];
  double vn = (va + vb + vc) / 3.0;
  dx[0] = (va - vn - R * x[0]) / L;
  dx[1] = (vb - vn - R * x[1]) / L;
}
static void carga_rl_out(BloqueC *bl, ModeloC *m, const double *x) {
  m->sig[bl->out_idx[0]] = x[0];          /* ia */
  m->sig[bl->out_idx[1]] = x[1];          /* ib */
  m->sig[bl->out_idx[2]] = -(x[0] + x[1]); /* ic */
}
/* ---------------- derivadas de máquinas ---------------- */
/* Saturacion magnetica opcional del PMAC.
   Param base: [rs, Ld, Lq, lam, P, J, Bm, ext]; con LUT se agrega al
   final: [n_pts, id_0, fl_0, ..., id_{n-1}, fl_{n-1}] donde fl es el
   flujo TOTAL de eje d (incluye el iman: fl(0) = lam). Sin LUT:
   lam_d = Ld*ids + lam, Lds = Ld. Con LUT: interpolacion lineal por
   tramos; fuera de rango se extrapola con el tramo extremo. */
static void pmac_lut(BloqueC *bl, double ids, double *lam_d, double *Lds) {
  double *p = bl->param;
  if (bl->n_param < 9) {
    *lam_d = p[1] * ids + p[3];
    *Lds = p[1];
    return;
  }
  int n = (int)p[8];
  double *t = p + 9; /* [id0, fl0, id1, fl1, ...] */
  if (n < 2) {       /* LUT degenerada: caer a lineal */
    *lam_d = p[1] * ids + p[3];
    *Lds = p[1];
    return;
  }
  if (ids <= t[0]) {
    double s = (t[3] - t[1]) / (t[2] - t[0]);
    *lam_d = t[1] + s * (ids - t[0]);
    *Lds = s;
    return;
  }
  for (int k = 0; k < n - 1; k++) {
    double ia = t[2 * k], fa = t[2 * k + 1];
    double ib = t[2 * k + 2], fb = t[2 * k + 3];
    if (ids <= ib) {
      double s = (fb - fa) / (ib - ia);
      *lam_d = fa + s * (ids - ia);
      *Lds = s;
      return;
    }
  }
  double ia = t[2 * (n - 2)], fa = t[2 * (n - 2) + 1];
  double ib = t[2 * (n - 1)], fb = t[2 * (n - 1) + 1];
  double s = (fb - fa) / (ib - ia);
  *lam_d = fb + s * (ids - ib);
  *Lds = s;
}
typedef struct {
  double *dx;
  double a, b, c; /* voltages abc en terminales */
  double TL;
  double vfd;
} mach_ctx;
static void maq_pmac_deriv(BloqueC *bl, ModeloC *m, const double *x, double *dx) {
  double *p = bl->param;
  double rs = p[0], Ld = p[1], Lq = p[2], lam = p[3], P = p[4], J = p[5], Bm = p[6];
  int ext = (bl->n_param > 7 && p[7] > 0.5);
  double iqs = x[0], ids = x[1];
  double wm = ext ? m->sig[bl->in_idx[3]] : x[2];
  double the = ext ? (P / 2.0) * m->sig[bl->in_idx[4]] : x[3];
  double we = (P / 2.0) * wm;
  double va = m->sig[bl->in_idx[0]], vb = m->sig[bl->in_idx[1]], vc = m->sig[bl->in_idx[2]];
  double al, be, vd, vq;
  f_clarke(va, vb, vc, &al, &be);
  f_park(al, be, the, &vd, &vq);
  double lam_d, Lds;
  pmac_lut(bl, ids, &lam_d, &Lds);
  dx[0] = (vq - rs * iqs - we * lam_d) / Lq;
  dx[1] = (vd - rs * ids + we * Lq * iqs) / Lds;
  if (!ext) {
    double TL = m->sig[bl->in_idx[3]];
    double Te = 1.5 * (P / 2.0) * (lam_d * iqs - Lq * iqs * ids);
    dx[2] = (Te - TL - Bm * wm) / J;
    dx[3] = we;
  } else {
    dx[2] = 0.0; /* ws es memoria compartida: no dejar residuos en ext */
    dx[3] = 0.0;
  }
}
static void maq_pmac_out(BloqueC *bl, ModeloC *m, const double *x) {
  double Ld = bl->param[1], Lq = bl->param[2], lam = bl->param[3], P = bl->param[4];
  int ext = (bl->n_param > 7 && bl->param[7] > 0.5);
  double iqs = x[0], ids = x[1];
  double wm = ext ? m->sig[bl->in_idx[3]] : x[2];
  double the = ext ? (P / 2.0) * m->sig[bl->in_idx[4]] : x[3];
  double al, be, ia, ib, ic;
  f_inv_park(ids, iqs, the, &al, &be);
  f_inv_clarke(al, be, &ia, &ib, &ic);
  m->sig[bl->out_idx[0]] = ia;
  m->sig[bl->out_idx[1]] = ib;
  m->sig[bl->out_idx[2]] = ic;
  m->sig[bl->out_idx[3]] = iqs;
  m->sig[bl->out_idx[4]] = ids;
  m->sig[bl->out_idx[5]] = wm;
  m->sig[bl->out_idx[6]] = the / (P / 2.0);
  m->sig[bl->out_idx[7]] = the;
  double lam_d, Lds;
  pmac_lut(bl, ids, &lam_d, &Lds);
  double Te = 1.5 * (P / 2.0) * (lam_d * iqs - Lq * iqs * ids);
  m->sig[bl->out_idx[8]] = Te;
  double P_cu = 1.5 * bl->param[0] * (iqs*iqs + ids*ids);
  if(bl->n_out > 9) m->sig[bl->out_idx[9]] = P_cu;
}
static void maq_ind_deriv(BloqueC *bl, ModeloC *m, const double *x, double *dx) {
  double *p = bl->param;
  double rs = p[0], rr = p[1], Li00 = p[2], Li01 = p[3], Li11 = p[4];
  double wf = p[5], P = p[6], J = p[7], Bm = p[8];
  int ext = (bl->n_param > 9 && p[9] > 0.5);
  double lqs = x[0], lds = x[1], lqr = x[2], ldr = x[3];
  double wm = ext ? m->sig[bl->in_idx[3]] : x[4];
  double we = (P / 2.0) * wm;
  double va = m->sig[bl->in_idx[0]], vb = m->sig[bl->in_idx[1]], vc = m->sig[bl->in_idx[2]];
  double al, be, vd, vq;
  f_clarke(va, vb, vc, &al, &be);
  f_park(al, be, 0.0, &vd, &vq);
  double iqs = Li00 * lqs + Li01 * lqr, iqr = Li01 * lqs + Li11 * lqr;
  double ids = Li00 * lds + Li01 * ldr, idr = Li01 * lds + Li11 * ldr;
  dx[0] = vq - rs * iqs - wf * lds;
  dx[1] = vd - rs * ids + wf * lqs;
  dx[2] = -rr * iqr - (wf - we) * ldr;
  dx[3] = -rr * idr + (wf - we) * lqr;
  if (!ext) {
    double TL = m->sig[bl->in_idx[3]];
    double Te = 1.5 * (P / 2.0) * (lds * iqs - lqs * ids);
    dx[4] = (Te - TL - Bm * wm) / J;
    dx[5] = wm;
  } else {
    dx[4] = 0.0;
    dx[5] = 0.0;
  }
}
static void maq_ind_out(BloqueC *bl, ModeloC *m, const double *x) {
  double P = bl->param[6];
  int ext = (bl->n_param > 9 && bl->param[9] > 0.5);
  double lqs = x[0], lds = x[1], lqr = x[2], ldr = x[3];
  double wm = ext ? m->sig[bl->in_idx[3]] : x[4];
  double thr = ext ? m->sig[bl->in_idx[4]] : x[5];
  double Li00 = bl->param[2], Li01 = bl->param[3], Li11 = bl->param[4];
  double iqs = Li00 * lqs + Li01 * lqr, iqr = Li01 * lqs + Li11 * lqr;
  double ids = Li00 * lds + Li01 * ldr, idr = Li01 * lds + Li11 * ldr;
  double al = ids, be = iqs; /* marco estacionario θ=0 */
  double ia, ib, ic;
  f_inv_clarke(al, be, &ia, &ib, &ic);
  m->sig[bl->out_idx[0]] = ia;
  m->sig[bl->out_idx[1]] = ib;
  m->sig[bl->out_idx[2]] = ic;
  m->sig[bl->out_idx[3]] = iqs;
  m->sig[bl->out_idx[4]] = ids;
  m->sig[bl->out_idx[5]] = wm;
  m->sig[bl->out_idx[6]] = thr;
  m->sig[bl->out_idx[7]] = 0.0;
  double Te = 1.5 * (P / 2.0) * (lds * iqs - lqs * ids);
  m->sig[bl->out_idx[8]] = Te;
  double rs = bl->param[0], rr = bl->param[1];
  double P_cu = 1.5 * (rs*(iqs*iqs+ids*ids) + rr*(iqr*iqr+idr*idr));
  if(bl->n_out > 9) m->sig[bl->out_idx[9]] = P_cu;
  if(bl->n_out > 12){
    double alr, ber, iar, ibr, icr;
    f_inv_park(idr, iqr, thr, &alr, &ber);
    f_inv_clarke(alr, ber, &iar, &ibr, &icr);
    m->sig[bl->out_idx[10]] = iar;
    m->sig[bl->out_idx[11]] = ibr;
    m->sig[bl->out_idx[12]] = icr;
  }
}
static void maq_sinc_deriv(BloqueC *bl, ModeloC *m, const double *x, double *dx) {
  double *p = bl->param;
  double rs = p[0], rfd = p[1], rkq1 = p[2], rkq2 = p[3], rkd = p[4];
  double P = p[5], J = p[6], Bm = p[7];
  double *Liq = &p[8], *Lid = &p[17];
  int ext = (bl->n_param > 26 && p[26] > 0.5);
  double lqs = x[0], lkq1 = x[1], lkq2 = x[2];
  double lds = x[3], lfd = x[4], lkd = x[5];
  double wm = ext ? m->sig[bl->in_idx[4]] : x[6];
  double the = ext ? (P / 2.0) * m->sig[bl->in_idx[5]] : x[7];
  double wr = (P / 2.0) * wm;
  double va = m->sig[bl->in_idx[0]], vb = m->sig[bl->in_idx[1]], vc = m->sig[bl->in_idx[2]];
  double vfd = m->sig[bl->in_idx[3]];
  double al, be, vd, vq;
  f_clarke(va, vb, vc, &al, &be);
  f_park(al, be, the, &vd, &vq);
  double iqs = Liq[0] * lqs + Liq[1] * lkq1 + Liq[2] * lkq2;
  double ikq1 = Liq[3] * lqs + Liq[4] * lkq1 + Liq[5] * lkq2;
  double ikq2 = Liq[6] * lqs + Liq[7] * lkq1 + Liq[8] * lkq2;
  double ids = Lid[0] * lds + Lid[1] * lfd + Lid[2] * lkd;
  double ifd = Lid[3] * lds + Lid[4] * lfd + Lid[5] * lkd;
  double ikd = Lid[6] * lds + Lid[7] * lfd + Lid[8] * lkd;
  dx[0] = vq - rs * iqs - wr * lds;
  dx[1] = -rkq1 * ikq1;
  dx[2] = -rkq2 * ikq2;
  dx[3] = vd - rs * ids + wr * lqs;
  dx[4] = vfd - rfd * ifd;
  dx[5] = -rkd * ikd;
  if (!ext) {
    double TL = m->sig[bl->in_idx[4]];
    double Te = 1.5 * (P / 2.0) * (lds * iqs - lqs * ids);
    dx[6] = (Te - TL - Bm * wm) / J;
    dx[7] = wr;
  } else {
    dx[6] = 0.0;
    dx[7] = 0.0;
  }
}
static void maq_sinc_out(BloqueC *bl, ModeloC *m, const double *x) {
  double P = bl->param[5];
  int ext = (bl->n_param > 26 && bl->param[26] > 0.5);
  double lqs = x[0], lds = x[3];
  double wm = ext ? m->sig[bl->in_idx[4]] : x[6];
  double the = ext ? (P / 2.0) * m->sig[bl->in_idx[5]] : x[7];
  double *Liq = &bl->param[8], *Lid = &bl->param[17];
  double iqs = Liq[0] * lqs + Liq[1] * x[1] + Liq[2] * x[2];
  double ids = Lid[0] * lds + Lid[1] * x[4] + Lid[2] * x[5];
  double al, be, ia, ib, ic;
  f_inv_park(ids, iqs, the, &al, &be);
  f_inv_clarke(al, be, &ia, &ib, &ic);
  m->sig[bl->out_idx[0]] = ia;
  m->sig[bl->out_idx[1]] = ib;
  m->sig[bl->out_idx[2]] = ic;
  m->sig[bl->out_idx[3]] = iqs;
  m->sig[bl->out_idx[4]] = ids;
  m->sig[bl->out_idx[5]] = wm;
  m->sig[bl->out_idx[6]] = the / (P / 2.0);
  m->sig[bl->out_idx[7]] = the;
  double Te = 1.5 * (P / 2.0) * (lds * iqs - lqs * ids);
  m->sig[bl->out_idx[8]] = Te;
  double rs = bl->param[0];
  double P_cu = 1.5 * rs * (iqs*iqs + ids*ids);
  if(bl->n_out > 9) m->sig[bl->out_idx[9]] = P_cu;
}
static void maq_cc_deriv(BloqueC *bl, ModeloC *m, const double *x, double *dx) {
  double *p = bl->param;
  double ra = p[0], La = p[1], rf = p[2], Lf = p[3], LAF = p[4], J = p[5], Bm = p[6];
  int ext = (bl->n_param > 7 && p[7] > 0.5);
  double ia = x[0], i_f = x[1];
  double wm = ext ? m->sig[bl->in_idx[2]] : x[2];
  double va = m->sig[bl->in_idx[0]], vf = m->sig[bl->in_idx[1]];
  double e = LAF * i_f * wm;
  dx[0] = (va - ra * ia - e) / La;
  dx[1] = (vf - rf * i_f) / Lf;
  if (!ext) {
    double TL = m->sig[bl->in_idx[2]];
    double Te = LAF * i_f * ia;
    dx[2] = (Te - TL - Bm * wm) / J;
    dx[3] = wm;
  } else {
    dx[2] = 0.0;
    dx[3] = 0.0;
  }
}
static void maq_cc_out(BloqueC *bl, ModeloC *m, const double *x) {
  int ext = (bl->n_param > 7 && bl->param[7] > 0.5);
  double ra = bl->param[0], LAF = bl->param[4];
  double wm = ext ? m->sig[bl->in_idx[2]] : x[2];
  m->sig[bl->out_idx[0]] = x[0];
  m->sig[bl->out_idx[1]] = x[1];
  m->sig[bl->out_idx[2]] = wm;
  m->sig[bl->out_idx[3]] = ext ? m->sig[bl->in_idx[3]] : x[3];
  m->sig[bl->out_idx[4]] = LAF * x[1] * x[0];
  double e = LAF * x[1] * wm;
  m->sig[bl->out_idx[5]] = e;
  m->sig[bl->out_idx[6]] = e + ra * x[0];
  double P_cu = ra*x[0]*x[0] + bl->param[2]*x[1]*x[1];
  if(bl->n_out > 7) m->sig[bl->out_idx[7]] = P_cu;
}
/* Maquina de CC de imanes permanentes (tipo Simulink "DC Machine", solo
   Ra, La y Kt = constante de par/FEM, Kt = Ke en SI):
   Param [r_a, L_a, Kt, J, Bm, ext]; estados [ia, wm, th_rm];
   entrada [va] + ([T_L] interna o [wm, th_rm] externa);
   salidas [ia, wm, th_rm, Te, Ea, V_t].
   Ea = Kt*wm ; Te = Kt*ia ; La*ia' = va - ra*ia - Ea ;
   J*wm' = Te - TL - Bm*wm ; th_rm' = wm. */
static void maq_dcpm_deriv(BloqueC *bl, ModeloC *m, const double *x, double *dx) {
  double *p = bl->param;
  double ra = p[0], La = p[1], Kt = p[2], J = p[3], Bm = p[4];
  int ext = (bl->n_param > 5 && p[5] > 0.5);
  double ia = x[0];
  double wm = ext ? m->sig[bl->in_idx[1]] : x[1];
  double e = Kt * wm;
  dx[0] = (m->sig[bl->in_idx[0]] - ra * ia - e) / La;
  if (!ext) {
    double TL = m->sig[bl->in_idx[1]];
    double Te = Kt * ia;
    dx[1] = (Te - TL - Bm * wm) / J;
    dx[2] = wm;
  } else {
    dx[1] = 0.0;
    dx[2] = 0.0;
  }
}
static void maq_dcpm_out(BloqueC *bl, ModeloC *m, const double *x) {
  int ext = (bl->n_param > 5 && bl->param[5] > 0.5);
  double ra = bl->param[0], Kt = bl->param[2];
  double wm = ext ? m->sig[bl->in_idx[1]] : x[1];
  m->sig[bl->out_idx[0]] = x[0];
  m->sig[bl->out_idx[1]] = wm;
  m->sig[bl->out_idx[2]] = ext ? m->sig[bl->in_idx[2]] : x[2];
  m->sig[bl->out_idx[3]] = Kt * x[0];
  double e = Kt * wm;
  m->sig[bl->out_idx[4]] = e;
  m->sig[bl->out_idx[5]] = e + ra * x[0];
  double P_cu = ra*x[0]*x[0];
  if(bl->n_out > 6) m->sig[bl->out_idx[6]] = P_cu;
}
static void eje_mecanico_deriv(BloqueC *bl, ModeloC *m, const double *x, double *dx) {
  double J_eq = bl->param[0];
  double Bm_eq = bl->param[1];
  double wm = x[0];
  double sum_Te = 0.0;
  for (int i=0; i<bl->n_in - 1; i++) {
    sum_Te += m->sig[bl->in_idx[i]];
  }
  double TL = m->sig[bl->in_idx[bl->n_in - 1]];
  dx[0] = (sum_Te - TL - Bm_eq * wm) / J_eq;
  dx[1] = wm;
}
static void eje_mecanico_out(BloqueC *bl, ModeloC *m, const double *x) {
  m->sig[bl->out_idx[0]] = x[0];
  m->sig[bl->out_idx[1]] = x[1];
}
/* ---------------- Bateria (Tremblay-Dessaint, tipo Simulink) ----------------
   Param [E0, K, Q, A, B, R, tau]; estados [it (Ah), i_f (A filtrada)];
   entrada [I] (A, positivo = descarga); salidas [Vbat, SOC].
   Modelo (paper EVS24 / SimPowerSystems Battery):
     descarga (i_f >= 0): E = E0 - K*Q*(it + i_f)/(Q - it) + A*e^(-B*it)
     carga (i_f < 0)    : E = E0 - K*Q*i_f/(it - 0.1Q) - K*Q*it/(Q - it)
                              + A*e^(-B*it)
     Vbat = E - R*i ;  SOC = 1 - it/Q  (valido con SOC en [10%, 100%])        */
static void bateria_deriv(BloqueC *bl, ModeloC *m, const double *x, double *dx) {
  double i = m->sig[bl->in_idx[0]];
  double Q = bl->param[2], tau = bl->param[6];
  double eta_c = bl->param[8];
  double it = x[0];
  double di = i / 3600.0;                 /* Ah/s */
  if (i < 0.0) di = i * eta_c / 3600.0;   /* eficiencia de carga (coumbica) */
  /* rango valido del modelo: it en [0, 0.9*Q] (SOC 100%..10%) */
  if ((it >= 0.9 * Q && di > 0.0) || (it <= 0.1 * Q && di < 0.0)) di = 0.0;
  dx[0] = di;
  dx[1] = (i - x[1]) / tau;
  /* histeresis (opcional): Exp' = B*|i|*(-Exp + A)  (siempre +A, como Simulink) */
  if (bl->param[9] > 0.5) {
    double A = bl->param[3], B = bl->param[4];
    dx[2] = B * fabs(i) * (-x[2] + A);
  } else {
    dx[2] = 0.0;
  }
  /* termica (opcional, C_th > 0): T' = (I^2*R(T) - (T - T_amb)/R_th)/C_th.
     R(T) = R*(1 + alpha*(T - 25 C)); p[10]=C_th [J/K], p[11]=R_th [K/W],
     p[12]=alpha [1/K], p[13]=T_amb [C]. */
  if (bl->n_param >= 14 && bl->param[10] > 0.0) {
    double R_eff = bl->param[5] * (1.0 + bl->param[12] * (x[3] - 25.0));
    dx[3] = (i * i * R_eff - (x[3] - bl->param[13]) / bl->param[11])
            / bl->param[10];
  } else {
    dx[3] = 0.0;
  }
}
static void bateria_out(BloqueC *bl, ModeloC *m, const double *x) {
  double E0 = bl->param[0], K = bl->param[1], Q = bl->param[2];
  double A = bl->param[3], B = bl->param[4];
  double Vcap = bl->param[7];
  double R = bl->param[5];
  if (bl->n_param >= 14 && bl->param[10] > 0.0)
    R = R * (1.0 + bl->param[12] * (x[3] - 25.0));
  double i = m->sig[bl->in_idx[0]];
  double it = x[0];
  if (it < 0.0) it = 0.0;
  if (it > 0.9 * Q) it = 0.9 * Q;
  double den1 = Q - it;              /* >= 0.1*Q: nunca singular */
  double den2 = it - 0.1 * Q;        /* -> 0 en SOC 90% (fin de carga) */
  if (den2 < 1e-9 * Q) den2 = 1e-9 * Q;
  double ifi = x[1];
  double exp_h = (bl->param[9] > 0.5) ? x[2] : A * exp(-B * it);
  double E;
  if (ifi >= 0.0)
    E = E0 - K * Q * (it + ifi) / den1 + exp_h;
  else {
    E = E0 - K * Q * ifi / den2 - K * Q * it / den1 + exp_h;
    if (E > Vcap) E = Vcap;          /* fin de carga limitado (BMS) */
  }
  m->sig[bl->out_idx[0]] = E - R * i;   /* Vbat */
  m->sig[bl->out_idx[1]] = 1.0 - it / Q; /* SOC */
  m->sig[bl->out_idx[2]] = x[3];        /* T */
}

/* ---------------- Bateria ECM (2do Orden + Termico) ----------------
   Param: Q_nom, V_nom, V_min_cell, V_max_cell,
          R0, R1, C1, R2, C2, N_s, N_p,
          n_ocv, [soc0, ocv0, ...],
          I_chg_c, I_dch_c, T_min_c, T_max_c, T_min_d, T_max_d,
          R_th, C_th, T_amb, deg.
   Estados: SOC, V_rc1, V_rc2, T_pack
   Entrada: I_load (A, + descarga). Salida: [V_term, SOC, T, Ploss, Ichg_lim, Idch_lim]
*/
static double bateria_ecm_ocv(double soc, double *lut, int n) {
  if (soc <= lut[0]) return lut[1];
  if (soc >= lut[2*(n-1)]) return lut[2*(n-1)+1];
  int lo=0, hi=n-1;
  while(hi-lo>1){ int mid=(lo+hi)/2; if(lut[2*mid] <= soc) lo=mid; else hi=mid; }
  double s0 = lut[2*lo], v0 = lut[2*lo+1];
  double s1 = lut[2*(lo+1)], v1 = lut[2*(lo+1)+1];
  if (s1 - s0 < 1e-12) return v0;
  return v0 + (v1 - v0) * (soc - s0) / (s1 - s0);
}
static void bateria_ecm_deriv(BloqueC *bl, ModeloC *m, const double *x, double *dx) {
  double *p = bl->param;
  double Q_nom = p[0], R1 = p[5], C1 = p[6], R2 = p[7], C2 = p[8];
  double N_p = p[10];
  int n_ocv = (int)p[11];
  double *post = &p[12 + 2*n_ocv];
  double R_th = post[6], C_th = post[7], T_amb = post[8];
  
  double I_load = m->sig[bl->in_idx[0]];
  double I_cell = I_load / N_p;
  
  /* Derivadas electricas */
  dx[0] = -I_cell / (Q_nom * 3600.0);       /* dSOC/dt */
  dx[1] = (I_cell - x[1] / R1) / C1;        /* dV_rc1/dt */
  dx[2] = (I_cell - x[2] / R2) / C2;        /* dV_rc2/dt */
  
  /* Modelo termico (paquete completo) */
  double R0 = p[4];
  double P_loss_cell = I_cell * I_cell * R0;
  if (R1 > 1e-6) P_loss_cell += x[1] * x[1] / R1;
  if (R2 > 1e-6) P_loss_cell += x[2] * x[2] / R2;
  double P_loss_pack = P_loss_cell * p[9] * p[10];
  
  dx[3] = (P_loss_pack - (x[3] - T_amb) / R_th) / C_th; /* dT/dt */
}

static void bateria_ecm_out(BloqueC *bl, ModeloC *m, const double *x) {
  double *p = bl->param;
  double R0 = p[4], R1 = p[5], R2 = p[7];
  double N_s = p[9], N_p = p[10];
  int n_ocv = (int)p[11];
  double *lut = &p[12];
  double *post = &p[12 + 2*n_ocv];
  
  double I_load = m->sig[bl->in_idx[0]];
  double I_cell = I_load / N_p;
  double soc = x[0];
  if (soc < 0.0) soc = 0.0;
  if (soc > 1.0) soc = 1.0;
  
  double ocv = bateria_ecm_ocv(soc, lut, n_ocv);
  double V_cell = ocv - I_cell * R0 - x[1] - x[2];
  double V_term = V_cell * N_s;
  
  double P_loss_cell = I_cell * I_cell * R0;
  if (R1 > 1e-6) P_loss_cell += x[1] * x[1] / R1;
  if (R2 > 1e-6) P_loss_cell += x[2] * x[2] / R2;
  double P_loss_pack = P_loss_cell * N_s * N_p;
  
  double T = x[3];
  double I_chg_c = post[0], I_dch_c = post[1];
  double Tmin_c = post[2], Tmax_c = post[3], Tmin_d = post[4], Tmax_d = post[5];
  
  double I_chg_lim = I_chg_c * N_p;
  double I_dch_lim = I_dch_c * N_p;
  
  /* Derating basico de limites */
  if (T < Tmin_c || T > Tmax_c || soc >= 1.0) I_chg_lim = 0.0;
  if (T < Tmin_d || T > Tmax_d || soc <= 0.0) I_dch_lim = 0.0;
  
  m->sig[bl->out_idx[0]] = V_term;
  m->sig[bl->out_idx[1]] = soc;
  m->sig[bl->out_idx[2]] = T;
  m->sig[bl->out_idx[3]] = P_loss_pack;
  m->sig[bl->out_idx[4]] = I_chg_lim;
  m->sig[bl->out_idx[5]] = I_dch_lim;
}


static void mutual_deriv(BloqueC *bl, ModeloC *m, const double *x, double *dx){
  double L1=bl->param[0], L2=bl->param[1], M=bl->param[2];
  double v1 = (bl->n_in>=1 && bl->in_idx[0]>=0) ? m->sig[bl->in_idx[0]] : 0.0;
  double v2 = (bl->n_in>=2 && bl->in_idx[1]>=0) ? m->sig[bl->in_idx[1]] : 0.0;
  double det = L1*L2 - M*M;
  if(fabs(det) < 1e-12) det = (det>=0?1e-12:-1e-12);
  dx[0] = ( L2*v1 - M*v2)/det;
  dx[1] = (-M*v1 + L1*v2)/det;
}
static void mutual_out(BloqueC *bl, ModeloC *m, const double *x){
  m->sig[bl->out_idx[0]] = x[0];
  m->sig[bl->out_idx[1]] = x[1];
}

/* ================== SubredMNA ==========================================
 * Contenedor MNA con tres modos:
 *  A) Pre-cálculo (ultra-rápido): 2^N matrices Bx/Bu pre-invertidas en Python
 *  B) LU disperso por paso (Dommel) : G+C/dt + b = W·u + C/dt·x_prev
 *  C) Zero-crossing: sub-pasos cuando switch/diodo conmuta dentro de dt
 * Layout param v2:
 *  p[0]=n_x, p[1]=n_u, p[2]=n_sw, p[3]=n_out, p[4]=n_dio, p[5]=n_mv,
 *  p[6..]=Vpairs, n_mi,Iidx, swCtrl, diodos(5), n_R,R*3, n_C,C*3, n_L*4, n_VS*4, n_IS*3, swTopo*4,
 *  [flag_precomp, num_states, Bx/Bu...] opcional
 */
#define MNA_MAX_NX 128
#define MNA_MAX_SW 20

// --- Helper sparse CSR (stub para KLU) ---
typedef struct { int n; int nnz; int *ia; int *ja; double *a; } csr_t;
// Por ahora wrapper denso; para KLU real enlazar SuiteSparse/KLU y reemplazar lu_solve_dense por klu_solve
static int lu_solve_dense(double LU[128][128], double bp[128], double x[128], int n){
  for(int k=0;k<n;k++){
    int piv=k; double maxv=fabs(LU[k][k]);
    for(int i=k+1;i<n;i++){ double v=fabs(LU[i][k]); if(v>maxv){maxv=v; piv=i;}}
    if(maxv<1e-14){ LU[k][k]+=1e-12; maxv=fabs(LU[k][k]);}
    if(piv!=k){ for(int j=k;j<n;j++){double t=LU[k][j]; LU[k][j]=LU[piv][j]; LU[piv][j]=t;} double t=bp[k]; bp[k]=bp[piv]; bp[piv]=t; }
    double Akk=LU[k][k]; if(fabs(Akk)<1e-18) Akk=(Akk>=0?1e-12:-1e-12);
    for(int i=k+1;i<n;i++){ double f=LU[i][k]/Akk; LU[i][k]=f; for(int j=k+1;j<n;j++) LU[i][j]-=f*LU[k][j]; bp[i]-=f*bp[k]; }
  }
  for(int i=n-1;i>=0;i--){ double s=bp[i]; for(int j=i+1;j<n;j++) s-=LU[i][j]*x[j]; double d=LU[i][i]; if(fabs(d)<1e-18) d=(d>=0?1e-12:-1e-12); x[i]=s/d; }
  return 0;
}

static void update_mna(BloqueC *bl, ModeloC *m) {
  double *p = bl->param;
  int nx   = (int)p[0];
  int nu   = (int)p[1];
  int n_sw_ctrl = (int)p[2];
  int n_out= (int)p[3];
  int n_diodos = (int)p[4];
  int n_mv = (int)p[5];
  if (nx <=0 || nx > MNA_MAX_NX) return;
  if (nu <=0) nu = 1;
  int p_mv_base = 6;
  int n_mi = (int)p[p_mv_base + 2*n_mv];
  int idx_sw_base = p_mv_base + 2*n_mv + 1 + n_mi;
  int idx_diodo_base = idx_sw_base + n_sw_ctrl;
  int pos = idx_diodo_base + 5 * n_diodos;
  
  int n_R = (int)p[pos]; pos++;
  int r_base = pos; pos += 3*n_R;
  int n_C = (int)p[pos]; pos++;
  int c_base = pos; pos += 3*n_C;
  int n_L = (int)p[pos]; pos++;
  int l_base = pos; pos += 4*n_L;
  int n_VS = (int)p[pos]; pos++;
  int vs_base = pos; pos += 4*n_VS;
  int n_IS = (int)p[pos]; pos++;
  int is_base = pos; pos += 3*n_IS;
  int n_SW_topo = (int)p[pos]; pos++;
  int sw_topo_base = pos; pos += 4*n_SW_topo;
  
  int n_vcvs = (int)p[pos]; pos++;
  int vcvs_base = pos; pos += 5*n_vcvs;
  int n_vccs = (int)p[pos]; pos++;
  int vccs_base = pos; pos += 4*n_vccs;
  int n_mut = (int)p[pos]; pos++;
  int mut_base = pos; pos += 9*n_mut;
  
  int metodo = 0;
  if(pos < bl->n_param) { pos++; } // Forzamos metodo=0 (Backward Euler) porque Python no aloca historial Trapezoidal
  
  int is_precomputed = 0;
  int num_states = 0;
  double *pre_base = NULL;
  if (pos < bl->n_param) {
    int flag = (int)p[pos]; pos++;
    if (flag==1 && pos < bl->n_param) {
      is_precomputed = 1;
      num_states = (int)p[pos]; pos++;
      pre_base = &p[pos];
    }
  }

  int sw_state[64]={0};
  for(int s=0;s<n_sw_ctrl;s++){
    int idx_u = (int)p[idx_sw_base + s];
    double cur=0;
    if(idx_u>=0 && idx_u < bl->n_in){
      long long si = bl->in_idx[idx_u];
      if(si>=0 && si < m->n_sig) cur = m->sig[si];
    }
    sw_state[s]=(cur>0.5)?1:0;
  }
  
  double x_prev[128]; for(int i=0;i<nx;i++) x_prev[i]=bl->state[i];
  double u_local[128]={0};
  for(int j=0;j<nu && j<128;j++) if(j < bl->n_in){ long long si=bl->in_idx[j]; if(si>=0 && si < m->n_sig) u_local[j]=m->sig[si]; }
  
  double dt = bl->dt; if(dt<1e-12) dt=1e-12; 
  double cur_inv = (metodo == 1) ? (2.0 / dt) : (1.0 / dt);

  if(is_precomputed && num_states>0){
    int diode_state[64]={0};
    for(int d=0;d<n_diodos;d++){
      int n1=(int)p[idx_diodo_base+5*d]; int n2=(int)p[idx_diodo_base+5*d+1]; double Vf=p[idx_diodo_base+5*d+2];
      double Va=(n1>=0 && n1<nx)? x_prev[n1]:0.0; double Vc=(n2>=0 && n2<nx)? x_prev[n2]:0.0;
      diode_state[d]=((Va-Vc)>Vf)?1:0;
    }
    
    double x_new[128]={0};
    for(int iter=0; iter<5; iter++){
      int estado=0;
      for(int s=0;s<n_sw_ctrl;s++) if(sw_state[s]) estado|=(1<<s);
      for(int d=0;d<n_diodos;d++) if(diode_state[d]) estado|=(1<<(n_sw_ctrl+d));
      if(estado<0) estado=0; if(estado>=num_states) estado=num_states-1;
      
      int Bx_sz = nx*nx; int Bu_sz = nx*nu; int blk = Bx_sz+Bu_sz;
      double *Bx = pre_base + estado*blk;
      double *Bu = pre_base + estado*blk + Bx_sz;
      double x_tmp[128]={0};
      
      for(int i=0;i<nx;i++){
        double acc=0;
        for(int j=0;j<nx;j++) acc+= Bx[i*nx+j]*x_prev[j];
        for(int j=0;j<nu;j++) acc+= Bu[i*nu+j]*u_local[j];
        x_tmp[i]=acc;
      }
      // Vf correction
      double b_diode_pre[128]={0}; int need_vf=0;
      for(int d=0;d<n_diodos;d++) if(diode_state[d]){ double Vf=p[idx_diodo_base+5*d+2]; if(fabs(Vf)>1e-12) need_vf=1; int n1=(int)p[idx_diodo_base+5*d]; int n2=(int)p[idx_diodo_base+5*d+1]; double Ron=p[idx_diodo_base+5*d+3]; double g=1.0/Ron; if(n1>=0) b_diode_pre[n1]+=g*Vf; if(n2>=0) b_diode_pre[n2]-=g*Vf; }
      if(need_vf){
#ifdef THREAD_LOCAL_FALLBACK_CALLOC
        double *Gp_data = (double*)calloc(128*128, sizeof(double));
        double *Cp_data = (double*)calloc(128*128, sizeof(double));
        double *Ap_data = (double*)calloc(128*128, sizeof(double));
        if(!Gp_data || !Cp_data || !Ap_data){ free(Gp_data); free(Cp_data); free(Ap_data); return; }
        double (*Gp)[128] = (double(*)[128])Gp_data;
        double (*Cp)[128] = (double(*)[128])Cp_data;
        double (*Ap)[128] = (double(*)[128])Ap_data;
#else
        double (*Gp)[128] = mna_G;
        double (*Cp)[128] = mna_Cmat;
        double (*Ap)[128] = mna_Ap;
        for(int _i=0;_i<nx;_i++){ for(int _j=0;_j<nx;_j++){ Gp[_i][_j]=0; Cp[_i][_j]=0; Ap[_i][_j]=0; } }
#endif
        for(int i=0;i<n_R;i++){ int n1=(int)p[r_base+3*i]; int n2=(int)p[r_base+3*i+1]; double R=p[r_base+3*i+2]; if(R<1e-12)R=1e-12; double g=1.0/R; if(n1>=0) Gp[n1][n1]+=g; if(n2>=0) Gp[n2][n2]+=g; if(n1>=0&&n2>=0){Gp[n1][n2]-=g; Gp[n2][n1]-=g;}}
        for(int i=0;i<n_C;i++){ int n1=(int)p[c_base+3*i]; int n2=(int)p[c_base+3*i+1]; double cval=p[c_base+3*i+2]; if(n1>=0) Cp[n1][n1]+=cval; if(n2>=0) Cp[n2][n2]+=cval; if(n1>=0&&n2>=0){Cp[n1][n2]-=cval; Cp[n2][n1]-=cval;}}
        for(int i=0;i<n_L;i++){ int n1=(int)p[l_base+4*i]; int n2=(int)p[l_base+4*i+1]; int id=(int)p[l_base+4*i+2]; double L=p[l_base+4*i+3]; if(L<1e-12)L=1e-12; if(n1>=0){Gp[n1][id]+=1; Gp[id][n1]+=1;} if(n2>=0){Gp[n2][id]-=1; Gp[id][n2]-=1;} if(id>=0&&id<nx) Cp[id][id]=-L; }
        for(int i=0;i<n_VS;i++){ int n1=(int)p[vs_base+4*i]; int n2=(int)p[vs_base+4*i+1]; int ii=(int)p[vs_base+4*i+3]; if(n1>=0){Gp[n1][ii]+=1; Gp[ii][n1]+=1;} if(n2>=0){Gp[n2][ii]-=1; Gp[ii][n2]-=1;} }
        for(int i=0;i<n_sw_ctrl;i++){ int n1=(int)p[sw_topo_base+4*i]; int n2=(int)p[sw_topo_base+4*i+1]; double Ron=p[sw_topo_base+4*i+2]; double Roff=p[sw_topo_base+4*i+3]; double R= (estado>>i)&1 ? Ron:Roff; if(R<1e-12)R=1e-12; double g=1.0/R; if(n1>=0) Gp[n1][n1]+=g; if(n2>=0) Gp[n2][n2]+=g; if(n1>=0&&n2>=0){Gp[n1][n2]-=g; Gp[n2][n1]-=g;}}
        for(int i=0;i<n_diodos;i++){ int n1=(int)p[idx_diodo_base+5*i]; int n2=(int)p[idx_diodo_base+5*i+1]; double Ron=p[idx_diodo_base+5*i+3]; double Roff=p[idx_diodo_base+5*i+4]; int bit=(estado>>(n_sw_ctrl+i))&1; double R= bit?Ron:Roff; if(R<1e-12)R=1e-12; double g=1.0/R; if(n1>=0) Gp[n1][n1]+=g; if(n2>=0) Gp[n2][n2]+=g; if(n1>=0&&n2>=0){Gp[n1][n2]-=g; Gp[n2][n1]-=g;}}
        for(int i=0;i<n_vcvs;i++){ int n1=(int)p[vcvs_base+5*i], n2=(int)p[vcvs_base+5*i+1], ii=(int)p[vcvs_base+5*i+3]; if(n1>=0){ Gp[n1][ii]+=1; Gp[ii][n1]+=1; } if(n2>=0){ Gp[n2][ii]-=1; Gp[ii][n2]-=1; } }
        for(int i=0;i<n_mut;i++){ int n1=(int)p[mut_base+9*i], n2=(int)p[mut_base+9*i+1], n3=(int)p[mut_base+9*i+2], n4=(int)p[mut_base+9*i+3], i1=(int)p[mut_base+9*i+4], i2=(int)p[mut_base+9*i+5]; double L1=p[mut_base+9*i+6], L2=p[mut_base+9*i+7], M=p[mut_base+9*i+8]; if(n1>=0){ Gp[n1][i1]+=1; Gp[i1][n1]+=1; } if(n2>=0){ Gp[n2][i1]-=1; Gp[i1][n2]-=1; } if(n3>=0){ Gp[n3][i2]+=1; Gp[i2][n3]+=1; } if(n4>=0){ Gp[n4][i2]-=1; Gp[i2][n4]-=1; } if(i1>=0 && i1<nx) Cp[i1][i1] = -L1; if(i2>=0 && i2<nx) Cp[i2][i2] = -L2; if(i1>=0 && i1<nx && i2>=0 && i2<nx){ Cp[i1][i2] = -M; Cp[i2][i1] = -M; } }
        
        for(int i=0;i<nx;i++) for(int j=0;j<nx;j++) Ap[i][j]=Gp[i][j]+Cp[i][j]/dt;
        for(int i=0;i<nx;i++){ double rs=0; for(int j=0;j<nx;j++) rs+=fabs(Ap[i][j]); if(rs<1e-14) Ap[i][i]+=1e-12; }
        double delta[128]={0}; lu_solve_dense(Ap,b_diode_pre,delta,nx); for(int i=0;i<nx;i++) x_tmp[i]+=delta[i];
#ifdef THREAD_LOCAL_FALLBACK_CALLOC
        free(Gp_data); free(Cp_data); free(Ap_data);
#endif
      }
      for(int i=0;i<nx;i++) x_new[i]=x_tmp[i];
      int changed=0;
      for(int d=0;d<n_diodos;d++){
        int n1=(int)p[idx_diodo_base+5*d]; int n2=(int)p[idx_diodo_base+5*d+1]; double Vf=p[idx_diodo_base+5*d+2]; double Ron=p[idx_diodo_base+5*d+3];
        double Va=(n1>=0 && n1<nx)? x_new[n1]:0.0; double Vc=(n2>=0 && n2<nx)? x_new[n2]:0.0;
        if(diode_state[d]==0){ if((Va-Vc)>Vf){ diode_state[d]=1; changed=1; } }
        else { double Id=(Va-Vc-Vf)/Ron; if(Id < -1e-9){ diode_state[d]=0; changed=1; } }
      }
      if(!changed) break;
    }
    for(int i=0;i<nx;i++) bl->state[i]=x_new[i];
    
  } else {
    // Modo B: LU disperso (MNA con VCVS, VCCS, Mutual Inductor)
    int diode_state[64]={0};
    for(int d=0;d<n_diodos;d++){ int n1=(int)p[idx_diodo_base+5*d]; int n2=(int)p[idx_diodo_base+5*d+1]; double Vf=p[idx_diodo_base+5*d+2]; double Va=(n1>=0&&n1<nx)? x_prev[n1]:0; double Vc=(n2>=0&&n2<nx)? x_prev[n2]:0; diode_state[d]=((Va-Vc)>Vf)?1:0; }
    
    double x_new[128]={0};
#ifdef THREAD_LOCAL_FALLBACK_CALLOC
    double *G_data = (double*)calloc(128*128, sizeof(double));
    double *Cmat_data = (double*)calloc(128*128, sizeof(double));
    double *W_data = (double*)calloc(128*128, sizeof(double));
    double *b_diode_data = (double*)calloc(128, sizeof(double));
    if(!G_data || !Cmat_data || !W_data || !b_diode_data){ free(G_data); free(Cmat_data); free(W_data); free(b_diode_data); return; }
    double (*G)[128] = (double(*)[128])G_data;
    double (*Cmat)[128] = (double(*)[128])Cmat_data;
    double (*W)[128] = (double(*)[128])W_data;
    double *b_diode = b_diode_data;
#else
    double (*G)[128] = mna_G;
    double (*Cmat)[128] = mna_Cmat;
    double (*W)[128] = mna_W;
    double *b_diode = mna_b_diode;
#endif

    for(int diter=0; diter<5; diter++){
      for(int _i=0;_i<nx;_i++){ b_diode[_i]=0; for(int _j=0;_j<nx;_j++){ G[_i][_j]=0; Cmat[_i][_j]=0; W[_i][_j]=0; } }
      
      for(int i=0;i<n_R;i++){ int n1=(int)p[r_base+3*i]; int n2=(int)p[r_base+3*i+1]; double R=p[r_base+3*i+2]; if(R<1e-12)R=1e-12; double g=1.0/R; if(n1>=0) G[n1][n1]+=g; if(n2>=0) G[n2][n2]+=g; if(n1>=0&&n2>=0){G[n1][n2]-=g; G[n2][n1]-=g;}}
      for(int i=0;i<n_C;i++){ int n1=(int)p[c_base+3*i]; int n2=(int)p[c_base+3*i+1]; double cval=p[c_base+3*i+2]; if(n1>=0) Cmat[n1][n1]+=cval; if(n2>=0) Cmat[n2][n2]+=cval; if(n1>=0&&n2>=0){Cmat[n1][n2]-=cval; Cmat[n2][n1]-=cval;}}
      for(int i=0;i<n_L;i++){ int n1=(int)p[l_base+4*i]; int n2=(int)p[l_base+4*i+1]; int id=(int)p[l_base+4*i+2]; double Lval=p[l_base+4*i+3]; if(Lval<1e-12)Lval=1e-12; if(n1>=0){G[n1][id]+=1; G[id][n1]+=1;} if(n2>=0){G[n2][id]-=1; G[id][n2]-=1;} if(id>=0&&id<nx) Cmat[id][id]=-Lval; }
      for(int i=0;i<n_VS;i++){ int n1=(int)p[vs_base+4*i]; int n2=(int)p[vs_base+4*i+1]; int iu=(int)p[vs_base+4*i+2]; int ii=(int)p[vs_base+4*i+3]; if(n1>=0){G[n1][ii]+=1; G[ii][n1]+=1;} if(n2>=0){G[n2][ii]-=1; G[ii][n2]-=1;} if(ii>=0&&ii<nx&&iu>=0&&iu<nu) W[ii][iu]=1.0; }
      for(int i=0;i<n_IS;i++){ int n1=(int)p[is_base+3*i]; int n2=(int)p[is_base+3*i+1]; int iu=(int)p[is_base+3*i+2]; if(n1>=0&&iu>=0&&iu<nu) W[n1][iu]+=-1.0; if(n2>=0&&iu>=0&&iu<nu) W[n2][iu]+=1.0; }
      for(int i=0;i<n_SW_topo;i++){ int n1=(int)p[sw_topo_base+4*i]; int n2=(int)p[sw_topo_base+4*i+1]; double Ron=p[sw_topo_base+4*i+2]; double Roff=p[sw_topo_base+4*i+3]; double R= sw_state[i]?Ron:Roff; if(R<1e-12)R=1e-12; double g=1.0/R; if(n1>=0) G[n1][n1]+=g; if(n2>=0) G[n2][n2]+=g; if(n1>=0&&n2>=0){G[n1][n2]-=g; G[n2][n1]-=g;}}
      for(int i=0;i<n_diodos;i++){ int n1=(int)p[idx_diodo_base+5*i]; int n2=(int)p[idx_diodo_base+5*i+1]; double Vf=p[idx_diodo_base+5*i+2]; double Ron=p[idx_diodo_base+5*i+3]; double Roff=p[idx_diodo_base+5*i+4]; double R= diode_state[i]?Ron:Roff; if(R<1e-12)R=1e-12; double g=1.0/R; if(n1>=0) G[n1][n1]+=g; if(n2>=0) G[n2][n2]+=g; if(n1>=0&&n2>=0){G[n1][n2]-=g; G[n2][n1]-=g;} if(diode_state[i]){ if(n1>=0) b_diode[n1]+=g*Vf; if(n2>=0) b_diode[n2]-=g*Vf;}}
      
      for(int i=0;i<n_vcvs;i++){
        int n1=(int)p[vcvs_base+5*i], n2=(int)p[vcvs_base+5*i+1];
        int iu=(int)p[vcvs_base+5*i+2], ii=(int)p[vcvs_base+5*i+3];
        double gain=p[vcvs_base+5*i+4];
        if(n1>=0){ G[n1][ii]+=1; G[ii][n1]+=1; }
        if(n2>=0){ G[n2][ii]-=1; G[ii][n2]-=1; }
        if(ii>=0 && ii<nx && iu>=0 && iu<nu) W[ii][iu] = gain;
      }
      for(int i=0;i<n_vccs;i++){
        int n1=(int)p[vccs_base+4*i], n2=(int)p[vccs_base+4*i+1];
        int iu=(int)p[vccs_base+4*i+2]; double gm=p[vccs_base+4*i+3];
        if(n1>=0 && iu>=0 && iu<nu) W[n1][iu] -= gm;
        if(n2>=0 && iu>=0 && iu<nu) W[n2][iu] += gm;
      }
      for(int i=0;i<n_mut;i++){
        int n1=(int)p[mut_base+9*i], n2=(int)p[mut_base+9*i+1];
        int n3=(int)p[mut_base+9*i+2], n4=(int)p[mut_base+9*i+3];
        int i1=(int)p[mut_base+9*i+4], i2=(int)p[mut_base+9*i+5];
        double L1=p[mut_base+9*i+6], L2=p[mut_base+9*i+7], M=p[mut_base+9*i+8];
        if(n1>=0){ G[n1][i1]+=1; G[i1][n1]+=1; }
        if(n2>=0){ G[n2][i1]-=1; G[i1][n2]-=1; }
        if(n3>=0){ G[n3][i2]+=1; G[i2][n3]+=1; }
        if(n4>=0){ G[n4][i2]-=1; G[i2][n4]-=1; }
        if(i1>=0 && i1<nx) Cmat[i1][i1] = -L1;
        if(i2>=0 && i2<nx) Cmat[i2][i2] = -L2;
        if(i1>=0 && i1<nx && i2>=0 && i2<nx){ Cmat[i1][i2] = -M; Cmat[i2][i1] = -M; }
      }

      double (*A)[128] = mna_Ap;
      double *b = mna_b;
#ifdef THREAD_LOCAL_FALLBACK_CALLOC
      double *A_data = (double*)calloc(128*128, sizeof(double));
      double *b_data = (double*)calloc(128, sizeof(double));
      if(!A_data || !b_data){ free(A_data); free(b_data); free(G_data); free(Cmat_data); free(W_data); free(b_diode_data); return; }
      A = (double(*)[128])A_data;
      b = b_data;
#endif
      for(int _i=0;_i<nx;_i++){ b[_i]=0; for(int _j=0;_j<nx;_j++){ A[_i][_j]=0; } }
      for(int i=0;i<nx;i++) for(int j=0;j<nx;j++) A[i][j]=G[i][j]+Cmat[i][j]*cur_inv;
      for(int i=0;i<nx;i++){ double rs=0; for(int j=0;j<nx;j++) rs+=fabs(A[i][j]); if(rs<1e-14) A[i][i]+=1e-12; }
      int n_nodes = nx - n_VS - n_L;
      for(int i=0;i<nx;i++){ double acc=b_diode[i]; for(int j=0;j<nu;j++) acc+=W[i][j]*u_local[j]; for(int j=0;j<nx;j++) acc+=Cmat[i][j]*cur_inv*x_prev[j]; if(metodo==1){ if(i < n_C){ acc += x_prev[nx + i]; } for(int _k=0; _k<n_L; _k++){ if(i == n_nodes + n_VS + _k){ double v_L_prev = x_prev[n_nodes + _k] - x_prev[n_nodes + _k + 1]; acc -= v_L_prev; break; } } } b[i]=acc; }
      lu_solve_dense(A,b,x_new,nx);
#ifdef THREAD_LOCAL_FALLBACK_CALLOC
      free(A_data); free(b_data);
#endif
      int ch=0;
      for(int d=0;d<n_diodos;d++){ int n1=(int)p[idx_diodo_base+5*d]; int n2=(int)p[idx_diodo_base+5*d+1]; double Vf=p[idx_diodo_base+5*d+2]; double Ron=p[idx_diodo_base+5*d+3]; double Va=(n1>=0&&n1<nx)?x_new[n1]:0; double Vc=(n2>=0&&n2<nx)?x_new[n2]:0; if(diode_state[d]==0){ if((Va-Vc)>Vf){diode_state[d]=1; ch=1;}} else { double Id=(Va-Vc-Vf)/Ron; if(Id < -1e-9){diode_state[d]=0; ch=1;}}}
      if(!ch) break;
    }
#ifdef THREAD_LOCAL_FALLBACK_CALLOC
    free(G_data); free(Cmat_data); free(W_data); free(b_diode_data);
#endif
    for(int i=0;i<nx;i++) bl->state[i]=x_new[i];
  }
  
  // Salidas
  for(int i=0;i<n_mv;i++){ int n1=(int)p[p_mv_base+2*i]; int n2=(int)p[p_mv_base+2*i+1]; double v1=(n1>=0&&n1<nx)?bl->state[n1]:0; double v2=(n2>=0&&n2<nx)?bl->state[n2]:0; if(i<bl->n_out) m->sig[bl->out_idx[i]]=v1-v2; }
  for(int i=0;i<n_mi;i++){ int id=(int)p[p_mv_base+2*n_mv+1+i]; double c=(id>=0&&id<nx)?bl->state[id]:0; if(n_mv+i<bl->n_out) m->sig[bl->out_idx[n_mv+i]]=c; }
}


/* ================= Helpers para LUT ================= */
static int busca_bp(const double *bp, int n, double u){
  if (u <= bp[0]) return 0;
  if (u >= bp[n-1]) return n-2;
  int lo=0, hi=n-1;
  while (hi-lo>1){
    int mid=(lo+hi)/2;
    if (bp[mid] <= u) lo=mid; else hi=mid;
  }
  return lo;
}

/* =============== Integrador =============== */
static void integrador_deriv(BloqueC *bl, ModeloC *m, const double *x, double *dx){
  (void)x;
  if (bl->n_in>=1 && bl->in_idx[0]>=0) dx[0]= m->sig[bl->in_idx[0]];
  else dx[0]=0.0;
}
static void integrador_out(BloqueC *bl, ModeloC *m, const double *x){
  m->sig[bl->out_idx[0]]= x[0];
}
static void limrap_deriv(BloqueC *bl, ModeloC *m, const double *x, double *dx){
  double up = bl->param[0], down = bl->param[1];
  double u = (bl->n_in>=1 && bl->in_idx[0]>=0) ? m->sig[bl->in_idx[0]] : 0.0;
  double pend = (u - x[0]) / (m->dt>1e-12? m->dt:1e-12);
  if (pend > up) pend = up;
  if (pend < -down) pend = -down;
  dx[0]= pend;
}
static void limrap_out(BloqueC *bl, ModeloC *m, const double *x){
  m->sig[bl->out_idx[0]]= x[0];
}
static void masa_termica_deriv(BloqueC *bl, ModeloC *m, const double *x, double *dx){
  double C = bl->param[0];
  double T_amb = bl->param[1];
  double R_amb = bl->param[2];
  double suma=0.0;
  for (int i=0;i<bl->n_in;i++) {
    int idx = (int)bl->in_idx[i];
    if (idx>=0) suma += m->sig[idx];
  }
  if (R_amb > 0.0) suma -= (x[0] - T_amb)/R_amb;
  dx[0]= suma / C;
}
static void masa_termica_out(BloqueC *bl, ModeloC *m, const double *x){
  m->sig[bl->out_idx[0]]= x[0];
}
static void eje_flexible_deriv(BloqueC *bl, ModeloC *m, const double *x, double *dx){
  (void)x;
  double w1 = (bl->n_in>=1 && bl->in_idx[0]>=0) ? m->sig[bl->in_idx[0]]:0.0;
  double w2 = (bl->n_in>=2 && bl->in_idx[1]>=0) ? m->sig[bl->in_idx[1]]:0.0;
  dx[0]= w1;
  dx[1]= w2;
}
static void eje_flexible_out(BloqueC *bl, ModeloC *m, const double *x){
  double K = bl->param[0], B = bl->param[1];
  double w1 = (bl->n_in>=1 && bl->in_idx[0]>=0) ? m->sig[bl->in_idx[0]]:0.0;
  double w2 = (bl->n_in>=2 && bl->in_idx[1]>=0) ? m->sig[bl->in_idx[1]]:0.0;
  m->sig[bl->out_idx[0]] = K*(x[0]-x[1]) + B*(w1 - w2);
}
static void vehiculo_deriv(BloqueC *bl, ModeloC *m, const double *x, double *dx){
  double mass=bl->param[0], Cd=bl->param[1], A=bl->param[2], Crr=bl->param[3], rho=bl->param[4], g=bl->param[5];
  double gr=bl->param[6], r_w=bl->param[7], eff_g=bl->param[8], eff_r=bl->param[9];
  double T_mot = (bl->n_in>=2 && bl->in_idx[1]>=0)? m->sig[bl->in_idx[1]]:0.0;
  double grade = 0.0;
  if (bl->n_in>=3 && bl->in_idx[2]>=0) grade = m->sig[bl->in_idx[2]];
  double v = x[0];
  double eff = (T_mot >=0.0)? eff_g : eff_r;
  double F_drive = T_mot * gr * eff / r_w;
  double F_aero = 0.5 * rho * Cd * A * v * fabs(v);
  double sign_v = (v > 1e-3) ? 1.0 : ((v < -1e-3) ? -1.0 : 0.0);
  double F_rod = mass * g * (Crr * sign_v * cos(grade) + sin(grade));
  dx[0]= (F_drive - F_aero - F_rod)/mass;
  if (v <= 1e-9 && v >= -1e-9 && dx[0] < 0.0 && F_drive >= -1e-9 && fabs(grade) < 1e-9) dx[0]=0.0;
}
static void vehiculo_out(BloqueC *bl, ModeloC *m, const double *x){
  double mass=bl->param[0], Cd=bl->param[1], A=bl->param[2], Crr=bl->param[3], rho=bl->param[4], g=bl->param[5];
  double gr=bl->param[6], r_w=bl->param[7], eff=bl->param[8];
  double v = x[0];
  double omega_m = (bl->n_in>=1 && bl->in_idx[0]>=0)? m->sig[bl->in_idx[0]]:0.0;
  double grade = 0.0;
  if (bl->n_in>=3 && bl->in_idx[2]>=0) grade = m->sig[bl->in_idx[2]];
  double F_aero = 0.5 * rho * Cd * A * v * fabs(v);
  double sign_v = (v > 1e-3) ? 1.0 : ((v < -1e-3) ? -1.0 : 0.0);
  double F_rod = mass * g * (Crr * sign_v * cos(grade) + sin(grade));
  double T_load = (F_aero + F_rod) * r_w / (gr * eff);
  m->sig[bl->out_idx[0]] = T_load;
  m->sig[bl->out_idx[1]] = v;
  m->sig[bl->out_idx[2]] = v * 3.6;
  m->sig[bl->out_idx[3]] = omega_m;
  double v_ref = omega_m * r_w / gr;
  double F_aero_ref = 0.5 * rho * Cd * A * v_ref * fabs(v_ref);
  double T_ff = (F_aero_ref + F_rod) * r_w / (gr * eff);
  m->sig[bl->out_idx[4]] = T_ff;
  m->sig[bl->out_idx[5]] = grade;
  m->sig[bl->out_idx[6]] = v * gr / r_w;
}
static void inductor_deriv(BloqueC *bl, ModeloC *m, const double *x, double *dx){
  (void)x;
  double L = bl->param[0];
  double v = (bl->n_in>=1 && bl->in_idx[0]>=0)? m->sig[bl->in_idx[0]]:0.0;
  dx[0]= v / L;
}
static void inductor_out(BloqueC *bl, ModeloC *m, const double *x){
  m->sig[bl->out_idx[0]]= x[0];
}
static void capacitor_deriv(BloqueC *bl, ModeloC *m, const double *x, double *dx){
  (void)x;
  double C = bl->param[0];
  double i = (bl->n_in>=1 && bl->in_idx[0]>=0)? m->sig[bl->in_idx[0]]:0.0;
  dx[0]= i / C;
}
static void capacitor_out(BloqueC *bl, ModeloC *m, const double *x){
  m->sig[bl->out_idx[0]]= x[0];
}

/* =============== TF y PID =============== */
static void update_tf(BloqueC *bl, ModeloC *m){
  int n = (int)bl->param[0];
  if (n<=0){ m->sig[bl->out_idx[0]]=0.0; return;}
  double *x = bl->state;
  double u = (bl->n_in>=1 && bl->in_idx[0]>=0)? m->sig[bl->in_idx[0]]:0.0;
  double *bd = &bl->param[1];
  double *ad = &bl->param[1+n+1];
  double y = bd[0]*u;
  for (int k=1;k<=n;k++){
    y += bd[k]* x[k-1] - ad[k-1]* x[n + k -1];
  }
  for (int k=n-1;k>=1;k--) x[k]= x[k-1];
  if (n>0) x[0]= u;
  for (int k=2*n-1;k>n;k--) x[k]= x[k-1];
  if (n>0) x[n]= y;
  m->sig[bl->out_idx[0]]= y;
}
static void update_pid(BloqueC *bl, ModeloC *m){
  double Kp=bl->param[0], Ki=bl->param[1], Kd=bl->param[2], Tf=bl->param[3], umin=bl->param[4], umax=bl->param[5];
  double *x = bl->state;
  double e = (bl->n_in>=1 && bl->in_idx[0]>=0)? m->sig[bl->in_idx[0]]:0.0;
  double dt = m->dt;
  double ud = (Kd*(e - x[1]) + Tf*x[2])/(Tf + dt);
  double u = Kp*e + Ki*x[0] + ud;
  if (! (u > umax && e > 0) && ! (u < umin && e < 0)){
    x[0] += dt*e;
  }
  x[1]= e;
  x[2]= ud;
  if (u > umax) u=umax;
  if (u < umin) u=umin;
  m->sig[bl->out_idx[0]]= u;
}

/* =============== Relay, Diodo, Retenedor, MaqEstados, PLL =============== */
static void update_relay(BloqueC *bl, ModeloC *m){
  double u = (bl->n_in>=1 && bl->in_idx[0]>=0)? m->sig[bl->in_idx[0]]:0.0;
  double on = bl->param[0], off = bl->param[1], y_on=bl->param[2], y_off=bl->param[3];
  double y = bl->state[0];
  if (u >= on) y = y_on;
  else if (u <= off) y = y_off;
  bl->state[0]= y;
  m->sig[bl->out_idx[0]]= y;
}
static void update_diodo(BloqueC *bl, ModeloC *m){
  double va = (bl->n_in>=1 && bl->in_idx[0]>=0)? m->sig[bl->in_idx[0]]:0.0;
  double vc = (bl->n_in>=2 && bl->in_idx[1]>=0)? m->sig[bl->in_idx[1]]:0.0;
  double r_on=bl->param[0], r_off=bl->param[1], Vf=bl->param[2], h=bl->param[3];
  double v = va - vc;
  int on = bl->state[0] > 0.5;
  if (v >= Vf + h) on=1;
  else if (v <= Vf - h) on=0;
  bl->state[0]= on?1.0:0.0;
  double r = on? r_on : r_off;
  if (r < 1e-12) r=1e-12;
  m->sig[bl->out_idx[0]]= v / r;
  m->sig[bl->out_idx[1]]= v;
}
static void update_retenedor(BloqueC *bl, ModeloC *m){
  double umb = bl->param[0];
  double u = (bl->n_in>=1 && bl->in_idx[0]>=0)? m->sig[bl->in_idx[0]]:0.0;
  double trig = (bl->n_in>=2 && bl->in_idx[1]>=0)? m->sig[bl->in_idx[1]]:0.0;
  double prev = bl->state[1];
  if (trig > umb && prev <= umb) bl->state[0]= u;
  bl->state[1]= trig;
  m->sig[bl->out_idx[0]]= bl->state[0];
}
static void update_maq_estados(BloqueC *bl, ModeloC *m){
  int n_trans = (int)bl->param[1];
  double *tr = &bl->param[2];
  int estado = (int)bl->state[0];
  for (int k=0;k<n_trans;k++){
    double from = tr[5*k], to=tr[5*k+1], sig=tr[5*k+2], opc=tr[5*k+3], umb=tr[5*k+4];
    if ((int)from != estado) continue;
    int sig_idx = (int)sig;
    double u = 0.0;
    if (sig_idx>=0 && sig_idx < bl->n_in){
      int gidx = (int)bl->in_idx[sig_idx];
      if (gidx>=0) u = m->sig[gidx];
    }
    int ok=0;
    int iopc=(int)opc;
    if (iopc==0) ok = (u < umb);
    else if (iopc==1) ok = (u <= umb);
    else if (iopc==2) ok = (u > umb);
    else if (iopc==3) ok = (u >= umb);
    else if (iopc==4) ok = (fabs(u - umb) <= 1e-12);
    else ok = (fabs(u - umb) > 1e-12);
    if (ok){ estado = (int)to; break; }
  }
  bl->state[0]= (double)estado;
  m->sig[bl->out_idx[0]]= (double)estado;
}
static void update_pll(BloqueC *bl, ModeloC *m){
  double h = m->dt;
  double Kp=bl->param[0], Ki=bl->param[1], wff=bl->param[2];
  double th = bl->state[0], xi = bl->state[1];
  double va = (bl->n_in>=1 && bl->in_idx[0]>=0)? m->sig[bl->in_idx[0]]:0.0;
  double vb = (bl->n_in>=2 && bl->in_idx[1]>=0)? m->sig[bl->in_idx[1]]:0.0;
  double vc = (bl->n_in>=3 && bl->in_idx[2]>=0)? m->sig[bl->in_idx[2]]:0.0;
  double al,be,d,q;
  f_clarke(va,vb,vc,&al,&be);
  f_park(al,be,th,&d,&q);
  double w = wff + Kp*q + Ki*xi;
  bl->state[1]= xi + h*q;
  bl->state[0]= th + h*w;
  m->sig[bl->out_idx[0]]= w;
  m->sig[bl->out_idx[1]]= bl->state[0];
}

/* =============== Estatico helpers =============== */
static void estatico_gain(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  double y = bl->param[0] * m->sig[bl->in_idx[0]];
  double o = m->sig[bl->out_idx[0]];
  double nv = o + w*(y - o);
  m->sig[bl->out_idx[0]]= nv;
  double d=fabs(nv - o); if (d>*maxdelta) *maxdelta=d;
}
static void estatico_sum(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w=m->w_opt;
  double y=0.0;
  for (int i=0;i<bl->n_in;i++) y += bl->param[i]* m->sig[bl->in_idx[i]];
  double o=m->sig[bl->out_idx[0]];
  double nv=o + w*(y - o);
  m->sig[bl->out_idx[0]]=nv;
  double d=fabs(nv - o); if(d>*maxdelta)*maxdelta=d;
}
static void pwm_1f_eval(BloqueC *bl, ModeloC *m, double *maxdelta){
  double fsw=bl->param[0], dead=bl->param[1];
  double td_frac = dead * fsw;
  double d = m->sig[bl->in_idx[0]];
  double d_eff = d - td_frac;
  if (d_eff <0.0) d_eff=0.0; if (d_eff>1.0) d_eff=1.0;
  double phi = fmod(m->t * fsw, 1.0); if (phi<0) phi+=1.0;
  double car = (phi <0.5)? (2.0*phi) : (2.0*(1.0-phi));
  double val = (car < d_eff) ? 1.0:0.0;
  double o = m->sig[bl->out_idx[0]];
  m->sig[bl->out_idx[0]]= val;
  double dlt=fabs(val - o); if(dlt>*maxdelta)*maxdelta=dlt;
}
static void pwm_spwm_eval(BloqueC *bl, ModeloC *m, double *maxdelta){
  double f_out=bl->param[0], fsw=bl->param[1], fase=bl->param[2], dead=bl->param[3];
  double td_frac = dead * fsw;
  double mod = m->sig[bl->in_idx[0]];
  if (mod<0.0) mod=0.0; if(mod>1.0) mod=1.0;
  double omega = 2.0*M_PI*f_out*m->t;
  double phi = fmod(m->t * fsw,1.0); if(phi<0) phi+=1.0;
  double car = (phi<0.5)? 2.0*phi : 2.0*(1.0-phi);
  for(int k=0;k<3;k++){
    double ref = 0.5 + 0.5*mod*sin(omega - k*2.0*M_PI/3.0 + fase);
    double d_eff = ref - td_frac;
    if(d_eff<0.0) d_eff=0.0; if(d_eff>1.0) d_eff=1.0;
    double val = (car < d_eff)?1.0:0.0;
    double o = m->sig[bl->out_idx[k]];
    m->sig[bl->out_idx[k]]= val;
    double dlt=fabs(val - o); if(dlt>*maxdelta)*maxdelta=dlt;
  }
}
static void pwm_svpwm_eval(BloqueC *bl, ModeloC *m, double *maxdelta){
  double Vdc=bl->param[0], fsw=bl->param[1], dead=bl->param[2];
  double td_frac = dead * fsw;
  double Va = m->sig[bl->in_idx[0]], Vb = m->sig[bl->in_idx[1]];
  double Vmax = Vdc / 1.7320508075688772;
  if (Vmax < 1e-9) { // Vdc=0: evita 0/0, fuerza salidas a 0 y no contamina con NaN
    for(int k=0;k<3;k++){ double o=m->sig[bl->out_idx[k]]; m->sig[bl->out_idx[k]]=0.0; double d=fabs(o); if(d>*maxdelta)*maxdelta=d; }
    return;
  }
  double Vmag = sqrt(Va*Va + Vb*Vb);
  double scale = 0.0;
  if (Vmag>1e-12) scale = fmin(Vmag, Vmax)/Vmag; else scale=0.0;
  double Vn_a = Va * scale / Vmax;
  double Vn_b = Vb * scale / Vmax;
  double theta = atan2(Vn_b, Vn_a);
  if (theta<0) theta+=2.0*M_PI;
  int sector = (int)(theta / (M_PI/3.0)) +1;
  if (sector<1) sector=1; if(sector>6) sector=6;
  double theta_s = theta - (sector-1)*(M_PI/3.0);
  double mag = sqrt(Vn_a*Vn_a + Vn_b*Vn_b);
  double T1 = mag * sin(M_PI/3.0 - theta_s);
  double T2 = mag * sin(theta_s);
  double T0 = 1.0 - T1 - T2; if(T0<0) T0=0.0;
  double h = T0*0.5;
  double da,db,dc;
  if (sector==1){ da=T1+T2+h; db=T2+h; dc=h; }
  else if(sector==2){ da=T1+h; db=T1+T2+h; dc=h; }
  else if(sector==3){ da=h; db=T1+T2+h; dc=T2+h; }
  else if(sector==4){ da=h; db=T1+h; dc=T1+T2+h; }
  else if(sector==5){ da=T2+h; db=h; dc=T1+T2+h; }
  else { da=T1+T2+h; db=h; dc=T1+h; }
  da -= td_frac; if(da<0) da=0; if(da>1) da=1;
  db -= td_frac; if(db<0) db=0; if(db>1) db=1;
  dc -= td_frac; if(dc<0) dc=0; if(dc>1) dc=1;
  double phi = fmod(m->t * fsw,1.0); if(phi<0) phi+=1.0;
  double car = (phi<0.5)? 2.0*phi : 2.0*(1.0-phi);
  double vals[3]={ (car<da)?1.0:0.0, (car<db)?1.0:0.0, (car<dc)?1.0:0.0 };
  for(int k=0;k<3;k++){
    double o=m->sig[bl->out_idx[k]];
    m->sig[bl->out_idx[k]]=vals[k];
    double dlt=fabs(vals[k]-o); if(dlt>*maxdelta)*maxdelta=dlt;
  }
}
static void hw_serial_eval(BloqueC *bl, ModeloC *m, double *maxdelta){
#ifdef _WIN32
  double w=m->w_opt;
  for(int k=0;k<bl->n_out;k++){
    double y=0.0;
    double o=m->sig[bl->out_idx[k]];
    double nv=o + w*(y - o);
    m->sig[bl->out_idx[k]]=nv;
    double d=fabs(nv - o); if(d>*maxdelta)*maxdelta=d;
  }
#else
  (void)bl;(void)m;(void)maxdelta;
  for(int k=0;k<bl->n_out;k++) m->sig[bl->out_idx[k]]=0.0;
#endif
}

/* =============== eval estatico general =============== */

/* =============== Static Block Implementations =============== */

static void estatico_mux(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  for(int k=0; k<bl->n_out; k++){
    double y = m->sig[bl->in_idx[k]];
    double o = m->sig[bl->out_idx[k]];
    double nv = o + w * (y - o);
    m->sig[bl->out_idx[k]] = nv;
    double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
  }
}

static void estatico_demux(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  for(int k=0; k<bl->n_out; k++){
    double y = m->sig[bl->in_idx[k]];
    double o = m->sig[bl->out_idx[k]];
    double nv = o + w * (y - o);
    m->sig[bl->out_idx[k]] = nv;
    double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
  }
}



static void estatico_lut1d(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  int n = (int)bl->param[0];
  double *bx = &bl->param[1];
  double *dy = &bl->param[1+n];
  double u = m->sig[bl->in_idx[0]];
  int i = busca_bp(bx, n, u);
  double f = (u - bx[i]) / (bx[i+1] - bx[i]);
  if (f < 0) f = 0; if (f > 1) f = 1;
  double y = dy[i] + f * (dy[i+1] - dy[i]);
  double o = m->sig[bl->out_idx[0]];
  double nv = o + m->w_opt * (y - o);
  m->sig[bl->out_idx[0]] = nv;
  double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
}

static void estatico_lut2d(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  int nx = (int)bl->param[0], ny = (int)bl->param[1];
  double *bx = &bl->param[2];
  double *by = &bl->param[2+nx];
  double *z = &bl->param[2+nx+ny];
  double u1 = m->sig[bl->in_idx[0]], u2 = m->sig[bl->in_idx[1]];
  int i = busca_bp(bx, nx, u1);
  int j = busca_bp(by, ny, u2);
  double fx = (u1 - bx[i]) / (bx[i+1] - bx[i]);
  double fy = (u2 - by[j]) / (by[j+1] - by[j]);
  if (fx < 0) fx = 0; if (fx > 1) fx = 1;
  if (fy < 0) fy = 0; if (fy > 1) fy = 1;
  double z00 = z[j*nx+i], z01 = z[j*nx+i+1];
  double z10 = z[(j+1)*nx+i], z11 = z[(j+1)*nx+i+1];
  double y = ((z00*(1-fx)+z01*fx)*(1-fy) + (z10*(1-fx)+z11*fx)*fy);
  double o = m->sig[bl->out_idx[0]];
  double nv = o + m->w_opt * (y - o);
  m->sig[bl->out_idx[0]] = nv;
  double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
}

static void estatico_lut3d(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  int nx = (int)bl->param[0], ny = (int)bl->param[1], nz = (int)bl->param[2];
  double *bx = &bl->param[3];
  double *by = &bl->param[3+nx];
  double *bz = &bl->param[3+nx+ny];
  double *z = &bl->param[3+nx+ny+nz];
  double u1 = m->sig[bl->in_idx[0]], u2 = m->sig[bl->in_idx[1]], u3 = m->sig[bl->in_idx[2]];
  int i = busca_bp(bx, nx, u1);
  int j = busca_bp(by, ny, u2);
  int k = busca_bp(bz, nz, u3);
  double fx = (u1 - bx[i]) / (bx[i+1] - bx[i]);
  double fy = (u2 - by[j]) / (by[j+1] - by[j]);
  double fz = (u3 - bz[k]) / (bz[k+1] - bz[k]);
  if (fx < 0) fx = 0; if (fx > 1) fx = 1;
  if (fy < 0) fy = 0; if (fy > 1) fy = 1;
  if (fz < 0) fz = 0; if (fz > 1) fz = 1;
  double y = 0.0;
  for(int kk=0; kk<2; kk++){
    double wz = (kk==0)?1-fz:fz;
    for(int jj=0; jj<2; jj++){
      double wy = (jj==0)?1-fy:fy;
      for(int ii=0; ii<2; ii++){
        double wx = (ii==0)?1-fx:fx;
        int idx = (k+kk)*ny*nx + (j+jj)*nx + (i+ii);
        y += wx*wy*wz*z[idx];
      }
    }
  }
  double o = m->sig[bl->out_idx[0]];
  double nv = o + m->w_opt * (y - o);
  m->sig[bl->out_idx[0]] = nv;
  double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
}

static void estatico_logico(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  int opc = (int)bl->param[0];
  double umb = bl->param[1];
  int r;
  if (opc == 6) { r = !(m->sig[bl->in_idx[0]] > bl->param[1]); }
  else {
    r = m->sig[bl->in_idx[0]] > bl->param[1];
    for(int k=1; k<bl->n_in; k++){
      int v = m->sig[bl->in_idx[k]] > bl->param[1];
      if (opc==0 || opc==2) r = r && v;
      else if (opc==1 || opc==3) r = r || v;
      else if (opc==4) r = r != v;
      else r = !(r != v);
    }
    if (opc==2 || opc==3) r = !r;
  }
  double y = r ? 1.0 : 0.0;
  double o = m->sig[bl->out_idx[0]];
  double nv = o + m->w_opt * (y - o);
  m->sig[bl->out_idx[0]] = nv;
  double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
}

static void estatico_relacional(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  int opc = (int)bl->param[0];
  double tol = bl->param[1];
  double a = m->sig[bl->in_idx[0]], b2 = m->sig[bl->in_idx[1]];
  int r = 0;
  if (opc==0) r = fabs(a-b2) <= tol;
  else if (opc==1) r = fabs(a-b2) > tol;
  else if (opc==2) r = a < b2;
  else if (opc==3) r = a <= b2;
  else if (opc==4) r = a > b2;
  else r = a >= b2;
  double y = r ? 1.0 : 0.0;
  double o = m->sig[bl->out_idx[0]];
  double nv = o + m->w_opt * (y - o);
  m->sig[bl->out_idx[0]] = nv;
  double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
}

static void estatico_clarke(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  double al, be;
  f_clarke(m->sig[bl->in_idx[0]], m->sig[bl->in_idx[1]], m->sig[bl->in_idx[2]], &al, &be);
  for(int j=0; j<2; j++){
    double val = (j==0)?al:be;
    double o = m->sig[bl->out_idx[j]];
    double nv = o + m->w_opt * (val - o);
    m->sig[bl->out_idx[j]] = nv;
    double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
  }
}

static void estatico_inv_clarke(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  double va, vb, vc;
  f_inv_clarke(m->sig[bl->in_idx[0]], m->sig[bl->in_idx[1]], &va, &vb, &vc);
  for(int j=0; j<3; j++){
    double val = (j==0)?va:((j==1)?vb:vc);
    double o = m->sig[bl->out_idx[j]];
    double nv = o + m->w_opt * (val - o);
    m->sig[bl->out_idx[j]] = nv;
    double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
  }
}

static void estatico_park(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  double d, q;
  f_park(m->sig[bl->in_idx[0]], m->sig[bl->in_idx[1]], m->sig[bl->in_idx[2]], &d, &q);
  for(int j=0; j<2; j++){
    double val = (j==0)?d:q;
    double o = m->sig[bl->out_idx[j]];
    double nv = o + m->w_opt * (val - o);
    m->sig[bl->out_idx[j]] = nv;
    double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
  }
}

static void estatico_inv_park(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  double al, be;
  f_inv_park(m->sig[bl->in_idx[0]], m->sig[bl->in_idx[1]], m->sig[bl->in_idx[2]], &al, &be);
  for(int j=0; j<2; j++){
    double val = (j==0)?al:be;
    double o = m->sig[bl->out_idx[j]];
    double nv = o + m->w_opt * (val - o);
    m->sig[bl->out_idx[j]] = nv;
    double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
  }
}

static void estatico_saturar(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  double u = m->sig[bl->in_idx[0]];
  double lo = bl->param[0], hi = bl->param[1];
  double y = (u < lo) ? lo : ((u > hi) ? hi : u);
  double o = m->sig[bl->out_idx[0]];
  double nv = o + m->w_opt * (y - o);
  m->sig[bl->out_idx[0]] = nv;
  double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
}

static void estatico_transformador(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  double a = bl->param[0];
  int fases = (bl->n_in > 3) ? 3 : 1;
  for(int k=0; k<fases; k++){
    double v1 = m->sig[bl->in_idx[k]];
    double i2 = m->sig[bl->in_idx[fases + k]];
    double y1 = a * v1;
    double y2 = -a * i2;
    double o1 = m->sig[bl->out_idx[k]];
    double nv1 = o1 + m->w_opt * (y1 - o1);
    m->sig[bl->out_idx[k]] = nv1;
    double d1 = fabs(nv1 - o1); if (d1 > *maxdelta) *maxdelta = d1;
    double o2 = m->sig[bl->out_idx[fases + k]];
    double nv2 = o2 + m->w_opt * (y2 - o2);
    m->sig[bl->out_idx[fases + k]] = nv2;
    double d2 = fabs(nv2 - o2); if (d2 > *maxdelta) *maxdelta = d2;
  }
}

static void estatico_panel(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  double G = m->sig[bl->in_idx[0]], T = m->sig[bl->in_idx[1]], V = m->sig[bl->in_idx[2]];
  double *p = bl->param;
  double Ns = p[0], Np = p[1], Voc = p[2], Isc = p[3], ki = p[6];
  double Rs = p[7], Rsh = p[8], n = p[9];
  double k_voc = (bl->n_param > 10) ? bl->param[10] : 0.0;
  double Vt = 0.02585 * (T + 273.15) / 298.15;
  double a = Ns * n * Vt;
  double Iph = Np * (Isc + ki * (T - 25.0)) * (G / 1000.0);
  double Voc_T = Voc * (1.0 + k_voc * (T - 25.0));
  double I0 = (Iph - Voc_T / Rsh) / (exp(Voc_T / a) - 1.0);
  double I = Iph - V / Rsh;
  for(int iter=0; iter<10; iter++){
    double u = (V + I * Rs) / a;
    if (u > 700.0) u = 700.0;
    double e = exp(u);
    double g = Iph - I0 * (e - 1.0) - (V + I * Rs) / Rsh - I;
    double gp = -I0 * (Rs / a) * e - Rs / Rsh - 1.0;
    double dI = g / gp;
    double lam = 1.0, g1 = fabs(g);
    for(int j=0; j<4; j++){
      double In = I - lam * dI;
      double u2 = (V + In * Rs) / a;
      if (u2 > 700.0) u2 = 700.0;
      double g2 = fabs(Iph - I0 * (exp(u2) - 1.0) - (V + In * Rs) / Rsh - In);
      if (g2 <= g1 || lam <= 0.125) break;
      lam *= 0.5;
    }
    I -= lam * dI;
  }
  double o = m->sig[bl->out_idx[0]];
  double nv = o + m->w_opt * (I - o);
  m->sig[bl->out_idx[0]] = nv;
  double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
}

static void estatico_medidor_potencia(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  int fases = (bl->n_in - 2) / 2;
  double Pe = 0.0;
  for(int k=0; k<fases; k++) Pe += m->sig[bl->in_idx[k]] * m->sig[bl->in_idx[fases + k]];
  double vals[3] = {Pe, 0.0, 0.0};
  int nvals = (bl->n_out == 3) ? 3 : 2;
  if (bl->n_out == 3){
    double v0=m->sig[bl->in_idx[0]], v1=m->sig[bl->in_idx[1]], v2=m->sig[bl->in_idx[2]];
    double i0=m->sig[bl->in_idx[3]], i1=m->sig[bl->in_idx[4]], i2=m->sig[bl->in_idx[5]];
    double vn=(v0+v1+v2)/3.0; double v0p=v0-vn, v1p=v1-vn, v2p=v2-vn;
    double in_= (i0+i1+i2)/3.0; double i0p=i0-in_, i1p=i1-in_, i2p=i2-in_;
    double SQRT23=0.8164965809277260, INV_SQRT2=0.7071067811865475;
    double va = SQRT23*(v0p -0.5*v1p -0.5*v2p); double vb = INV_SQRT2*(v1p - v2p);
    double ia = SQRT23*(i0p -0.5*i1p -0.5*i2p); double ib = INV_SQRT2*(i1p - i2p);
    vals[1] = vb*ia - va*ib; // Akagi Q
  }
  int idx_te = 2*fases;
  int idx_wm = 2*fases+1;
  double te = (bl->in_idx[idx_te] >= 0) ? m->sig[bl->in_idx[idx_te]] : 0.0;
  double wm = (bl->in_idx[idx_wm] >= 0) ? m->sig[bl->in_idx[idx_wm]] : 0.0;
  double Pm = te * wm;
  vals[nvals-1] = Pm;
  for(int j=0; j<nvals; j++){
    double o = m->sig[bl->out_idx[j]];
    double nv = o + m->w_opt * (vals[j] - o);
    m->sig[bl->out_idx[j]] = nv;
    double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
  }
}

static void estatico_interruptor(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  int modo = (int)bl->param[2];
  double r_on = bl->param[0], r_off = bl->param[1];
  double g = (bl->n_in==3) ? m->sig[bl->in_idx[0]] : 0.0;
  double va = m->sig[bl->in_idx[1]], vc = m->sig[bl->in_idx[2]];
  double r = r_off;
  if ((modo==0 && g>0.5) || (modo==1 && va-vc>0.0)) r = r_on;
  if (r<1e-12) r=1e-12;
  double v = va - vc;
  for(int k=0; k<2; k++){
    double val = (k==0) ? v/r : v;
    double o = m->sig[bl->out_idx[k]];
    double nv = o + m->w_opt * (val - o);
    m->sig[bl->out_idx[k]] = nv;
    double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
  }
}

static void estatico_puente_inv_3f(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  int prom = (bl->param[0] > 0.5);
  double vdc = m->sig[bl->in_idx[0]];
  double a = m->sig[bl->in_idx[1]], b2 = m->sig[bl->in_idx[2]], c2 = m->sig[bl->in_idx[3]];
  double vals[3];
  if (prom){
    vals[0] = vdc*(a - 0.5*b2 - 0.5*c2)/3.0;
    vals[1] = vdc*(b2 - 0.5*a - 0.5*c2)/3.0;
    vals[2] = vdc*(c2 - 0.5*a - 0.5*b2)/3.0;
  }else{
    int sa = (a>0.5), sb = (b2>0.5), sc = (c2>0.5);
    vals[0] = vdc*(2.0*sa - sb - sc)/3.0;
    vals[1] = vdc*(2.0*sb - sa - sc)/3.0;
    vals[2] = vdc*(2.0*sc - sa - sb)/3.0;
  }
  for(int k=0;k<3;k++){
    double o = m->sig[bl->out_idx[k]];
    double nv = o + m->w_opt * (vals[k] - o);
    m->sig[bl->out_idx[k]] = nv;
    double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
  }
}

static void estatico_puente_inv_1f(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  int prom = (bl->param[0] > 0.5);
  double vdc = m->sig[bl->in_idx[0]];
  double a = m->sig[bl->in_idx[1]], b2 = m->sig[bl->in_idx[2]];
  double y = prom ? vdc*(a - b2)/2.0 : vdc*(a - b2);
  double o = m->sig[bl->out_idx[0]];
  double nv = o + m->w_opt * (y - o);
  m->sig[bl->out_idx[0]] = nv;
  double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
}

static void estatico_resistencia_termica(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  double y = (m->sig[bl->in_idx[0]] - m->sig[bl->in_idx[1]]) / bl->param[0];
  double o = m->sig[bl->out_idx[0]];
  double nv = o + m->w_opt * (y - o);
  m->sig[bl->out_idx[0]] = nv;
  double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
}

static void estatico_engranaje(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  double a = bl->param[0];
  for(int k=0;k<2;k++){
    double y = (k==0) ? a * m->sig[bl->in_idx[0]] : m->sig[bl->in_idx[1]] / a;
    double o = m->sig[bl->out_idx[k]];
    double nv = o + m->w_opt * (y - o);
    m->sig[bl->out_idx[k]] = nv;
    double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
  }
}

static void estatico_embrague(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  double Tmax = bl->param[0], umb = bl->param[1];
  double T = m->sig[bl->in_idx[0]];
  double y = 0.0;
  if (m->sig[bl->in_idx[1]] > umb){
    y = T; if(y> Tmax) y=Tmax; if(y< -Tmax) y=-Tmax;
  } else y=0.0;
  double o = m->sig[bl->out_idx[0]];
  double nv = o + m->w_opt * (y - o);
  m->sig[bl->out_idx[0]] = nv;
  double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
}

static void estatico_fallo_prog(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  double t_fallo = bl->param[0], valor = bl->param[1], modo = bl->param[2];
  double u = m->sig[bl->in_idx[0]];
  double y = (m->t >= t_fallo) ? ((modo>0.5)? u+valor : valor) : u;
  double o = m->sig[bl->out_idx[0]];
  double nv = o + m->w_opt * (y - o);
  m->sig[bl->out_idx[0]] = nv;
  double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
}

static void estatico_fallo_evento(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  double umb = bl->param[0], valor = bl->param[1], modo = bl->param[2];
  double u = m->sig[bl->in_idx[0]];
  double y = (m->sig[bl->in_idx[1]] > umb) ? ((modo>0.5)? u+valor:valor) : u;
  double o = m->sig[bl->out_idx[0]];
  double nv = o + m->w_opt * (y - o);
  m->sig[bl->out_idx[0]] = nv;
  double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
}

static void estatico_multiplicador(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  double y = m->sig[bl->in_idx[0]] * m->sig[bl->in_idx[1]];
  double o = m->sig[bl->out_idx[0]];
  double nv = o + m->w_opt * (y - o);
  m->sig[bl->out_idx[0]] = nv;
  double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
}

static void estatico_sat_vectorial(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  double Vmax = bl->param[0];
  double vd = m->sig[bl->in_idx[0]], vq = m->sig[bl->in_idx[1]];
  double mag = sqrt(vd*vd + vq*vq);
  double scale = 1.0;
  if (mag > Vmax && mag > 1e-12) scale = Vmax / mag;
  for(int k=0; k<2; k++){
    double val = (k==0)? vd*scale : vq*scale;
    double o = m->sig[bl->out_idx[k]];
    double nv = o + m->w_opt * (val - o);
    m->sig[bl->out_idx[k]] = nv;
    double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
  }
}

static void estatico_calculo_idc(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  double vd = m->sig[bl->in_idx[0]], vq = m->sig[bl->in_idx[1]];
  double id = m->sig[bl->in_idx[2]], iq = m->sig[bl->in_idx[3]];
  double vdc = m->sig[bl->in_idx[4]];
  double eff = (bl->n_param>=1) ? bl->param[0] : 1.0;
  double idc = 0.0;
  if (vdc > 1e-3) idc = 1.5 * (vd * id + vq * iq) / (vdc * eff);
  double o = m->sig[bl->out_idx[0]];
  double nv = o + m->w_opt * (idc - o);
  m->sig[bl->out_idx[0]] = nv;
  double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
}

static void estatico_resistencia(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  double y = m->sig[bl->in_idx[0]] / bl->param[0];
  double o = m->sig[bl->out_idx[0]];
  double nv = o + m->w_opt * (y - o);
  m->sig[bl->out_idx[0]] = nv;
  double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
}

static void estatico_qd(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  double va = m->sig[bl->in_idx[0]], vb = m->sig[bl->in_idx[1]], vc = m->sig[bl->in_idx[2]];
  double ia = m->sig[bl->in_idx[3]], ib = m->sig[bl->in_idx[4]], ic = m->sig[bl->in_idx[5]], th = m->sig[bl->in_idx[6]];
  double alv, bev, ali, bei, d, q;
  f_clarke(va, vb, vc, &alv, &bev);
  f_clarke(ia, ib, ic, &ali, &bei);
  f_park(alv, bev, th, &d, &q); double vqs = q, vds = d;
  f_park(ali, bei, th, &d, &q); double iqs = q, ids = d;
  double vals[4] = {vqs, vds, iqs, ids};
  for(int k=0; k<4; k++){
    double o = m->sig[bl->out_idx[k]];
    double nv = o + m->w_opt * (vals[k] - o);
    m->sig[bl->out_idx[k]] = nv;
    double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
  }
}

static void estatico_carga_pq_3f(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  double va = m->sig[bl->in_idx[0]], vb = m->sig[bl->in_idx[1]], vc = m->sig[bl->in_idx[2]];
  double vn = (va + vb + vc) / 3.0;
  double v0 = va - vn, v1 = vb - vn, v2 = vc - vn;
  double SQRT23 = 0.8164965809277260, INV_SQRT2 = 0.7071067811865475, INV_SQRT6 = 0.4082482904638630;
  double valpha = SQRT23 * (v0 - 0.5 * v1 - 0.5 * v2);
  double vbeta = INV_SQRT2 * (v1 - v2);
  double vmag2 = valpha * valpha + vbeta * vbeta;
  double ia=0, ib=0, ic=0;
  if (vmag2 > 1e-6){
    double P_total = bl->param[0], Q_total = bl->param[1];
    double ialpha = (P_total * valpha + Q_total * vbeta) / vmag2;
    double ibeta  = (P_total * vbeta - Q_total * valpha) / vmag2;
    ia = SQRT23 * ialpha;
    ib = -INV_SQRT6 * ialpha + INV_SQRT2 * ibeta;
    ic = -INV_SQRT6 * ialpha - INV_SQRT2 * ibeta;
  }
  double vals[3] = {ia, ib, ic};
  for(int k=0;k<3;k++){
    double o = m->sig[bl->out_idx[k]];
    double nv = o + m->w_opt * (vals[k] - o);
    m->sig[bl->out_idx[k]] = nv;
    double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
  }
}

static void estatico_carga_pq_1f(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  double v = m->sig[bl->in_idx[0]];
  double P = bl->param[0], Q = bl->param[1];
  double mag2 = v * v;
  double i = 0.0;
  if (mag2 > 1e-6) i = P * v / mag2;
  double o = m->sig[bl->out_idx[0]];
  double nv = o + m->w_opt * (i - o);
  m->sig[bl->out_idx[0]] = nv;
  double d = fabs(nv - o); if (d > *maxdelta) *maxdelta = d;
}

static void estatico_vcvs(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  double vp = m->sig[bl->in_idx[0]];
  double vn = m->sig[bl->in_idx[1]];
  double gain = bl->param[0];
  double y = gain * (vp - vn);
  double o = m->sig[bl->out_idx[0]];
  double nv = o + m->w_opt * (y - o);
  m->sig[bl->out_idx[0]] = nv;
  if (fabs(nv - o) > *maxdelta) *maxdelta = fabs(nv - o);
}

static void estatico_vccs(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  double vp = m->sig[bl->in_idx[0]];
  double vn = m->sig[bl->in_idx[1]];
  double gm = bl->param[0];
  double y = gm * (vp - vn);
  double o = m->sig[bl->out_idx[0]];
  double nv = o + m->w_opt * (y - o);
  m->sig[bl->out_idx[0]] = nv;
  if (fabs(nv - o) > *maxdelta) *maxdelta = fabs(nv - o);
}

static void estatico_mutual_inductor(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  double L1 = bl->param[0], L2 = bl->param[1], M = bl->param[2];
  double v1 = m->sig[bl->in_idx[0]];
  double v2 = m->sig[bl->in_idx[1]];
  double det = L1 * L2 - M * M;
  double i1 = 0.0, i2 = 0.0;
  if (fabs(L1 * L2 - M * M) > 1e-12) {
    i1 = (v1 * L2 - v2 * M) / det;
    i2 = (v1 * M + v2 * L1) / det;
  }
  double o1 = m->sig[bl->out_idx[0]];
  double nv1 = o1 + m->w_opt * (i1 - o1);
  m->sig[bl->out_idx[0]] = nv1;
  if (fabs(nv1 - o1) > *maxdelta) *maxdelta = fabs(nv1 - o1);
  double o2 = m->sig[bl->out_idx[1]];
  double nv2 = o2 + m->w_opt * (i2 - o2);
  m->sig[bl->out_idx[1]] = nv2;
  if (fabs(nv2 - o2) > *maxdelta) *maxdelta = fabs(nv2 - o2);
}

static void eval_estatico(BloqueC *bl, ModeloC *m, double *maxdelta){
  double w = m->w_opt;
  int op = bl->op;
  switch(op){
    case OP_GAIN: { estatico_gain(bl,m,maxdelta); break; }
    case OP_SUM: { estatico_sum(bl,m,maxdelta); break; }
    case OP_MUX: { estatico_mux(bl,m,maxdelta); break; }
    case OP_DEMUX: { estatico_demux(bl,m,maxdelta); break; }
    case OP_LUT1D: { estatico_lut1d(bl,m,maxdelta); break; }
    case OP_LUT2D: { estatico_lut2d(bl,m,maxdelta); break; }
    case OP_LUT3D: { estatico_lut3d(bl,m,maxdelta); break; }
    case OP_LOGICO: { estatico_logico(bl,m,maxdelta); break; }
    case OP_RELACIONAL: { estatico_relacional(bl,m,maxdelta); break; }
    case OP_CLARKE: { estatico_clarke(bl,m,maxdelta); break; }
    case OP_INV_CLARKE: { estatico_inv_clarke(bl,m,maxdelta); break; }
    case OP_PARK: { estatico_park(bl,m,maxdelta); break; }
    case OP_INV_PARK: { estatico_inv_park(bl,m,maxdelta); break; }
    case OP_SATURAR: { estatico_saturar(bl,m,maxdelta); break; }
    case OP_TRANSFORMADOR: { estatico_transformador(bl,m,maxdelta); break; }
    case OP_PANEL_SOLAR: { estatico_panel(bl,m,maxdelta); break; }
    case OP_MEDIDOR_POTENCIA: { estatico_medidor_potencia(bl,m,maxdelta); break; }
    case OP_INTERRUPTOR: { estatico_interruptor(bl,m,maxdelta); break; }
    case OP_PUENTE_INV_3F: { estatico_puente_inv_3f(bl,m,maxdelta); break; }
    case OP_PUENTE_INV_1F: { estatico_puente_inv_1f(bl,m,maxdelta); break; }
    case OP_RES_TERMICA: { estatico_resistencia_termica(bl,m,maxdelta); break; }
    case OP_ENGRANAJE: { estatico_engranaje(bl,m,maxdelta); break; }
    case OP_EMBRAGUE: { estatico_embrague(bl,m,maxdelta); break; }
    case OP_FALLO_PROG: { estatico_fallo_prog(bl,m,maxdelta); break; }
    case OP_FALLO_EVENTO: { estatico_fallo_evento(bl,m,maxdelta); break; }
    case OP_MULTIPLICADOR: { estatico_multiplicador(bl,m,maxdelta); break; }
    case OP_SAT_VECTORIAL: { estatico_sat_vectorial(bl,m,maxdelta); break; }
    case OP_CALCULO_IDC: { estatico_calculo_idc(bl,m,maxdelta); break; }
    case OP_RESISTENCIA: { estatico_resistencia(bl,m,maxdelta); break; }
    case OP_QD: { estatico_qd(bl,m,maxdelta); break; }
    case OP_CARGA_PQ_3F: { estatico_carga_pq_3f(bl,m,maxdelta); break; }
    case OP_CARGA_PQ_1F: { estatico_carga_pq_1f(bl,m,maxdelta); break; }
    case OP_VCVS: { estatico_vcvs(bl,m,maxdelta); break; }
    case OP_VCCS: { estatico_vccs(bl,m,maxdelta); break; }
    case OP_MUTUAL_INDUCTOR: { estatico_mutual_inductor(bl,m,maxdelta); break; }
    default: break;
  }
}

static void actualizar_fuentes(ModeloC *m){
  for(int i=0;i<m->n_bloques;i++){
    BloqueC *bl=&m->bloques[i];
    switch(bl->op){
      case OP_SRC_CONST: {
        m->sig[bl->out_idx[0]]= bl->param[0]; break;
      }
      case OP_SRC_STEP: {
        double vf=bl->param[0], ts=bl->param[1], v0=bl->param[2];
        m->sig[bl->out_idx[0]]= (m->t >= ts)? vf : v0; break;
      }
      case OP_SRC_RAMP: {
        double k=bl->param[0], t0=bl->param[1], off=bl->param[2];
        m->sig[bl->out_idx[0]]= (m->t >= t0)? off + k*(m->t - t0) : off; break;
      }
      case OP_SRC_SIN: {
        double Am=bl->param[0], f=bl->param[1], ph=bl->param[2], off=bl->param[3];
        m->sig[bl->out_idx[0]]= off + Am*sin(2.0*M_PI*f*m->t + ph); break;
      }
      case OP_SRC_TRIF: {
        double Am=bl->param[0], f=bl->param[1], ph=bl->param[2];
        double w=2.0*M_PI*f*m->t;
        m->sig[bl->out_idx[0]]= Am*sin(w + ph);
        m->sig[bl->out_idx[1]]= Am*sin(w + ph - 2.0*M_PI/3.0);
        m->sig[bl->out_idx[2]]= Am*sin(w + ph + 2.0*M_PI/3.0);
        break;
      }
      case OP_PULSO_RECT: {
        double amp=bl->param[0], T=bl->param[1], duty=bl->param[2], ph=bl->param[3], off=bl->param[4];
        double tt = fmod(m->t + ph, T); if(tt<0) tt+=T;
        double y = off + ((tt < duty*T)? amp:0.0);
        m->sig[bl->out_idx[0]]= y; break;
      }
      case OP_SRC_CSV:
      case OP_SRC_TABLE: {
        int n=(int)bl->param[0];
        int interp = bl->param[1]>0.5;
        const double *tt=&bl->param[2];
        const double *yy=&bl->param[2+n];
        double t=m->t; double y;
        if(t <= tt[0]) y=yy[0];
        else if(t >= tt[n-1]) y=yy[n-1];
        else {
          int k=busca_bp(tt,n,t);
          if(interp){
            double f=(t - tt[k])/(tt[k+1]-tt[k]);
            y=yy[k] + f*(yy[k+1]-yy[k]);
          } else y=yy[k];
        }
        m->sig[bl->out_idx[0]]= y; break;
      }
      default: break;
    }
  }
}

/* =============== Emitir valor inicial (t=0) =============== */
static void emit_inicial(ModeloC *m){
  for(int i=0;i<m->n_bloques;i++){
    BloqueC *bl=&m->bloques[i];
    if(bl->op==OP_RELAY) update_relay(bl,m);
    else if(bl->op==OP_DIODO) update_diodo(bl,m);
    else if(bl->op==OP_RETENEDOR){
      m->sig[bl->out_idx[0]]= bl->state[0];
    }
    else if(bl->op==OP_MAQ_ESTADOS){
      m->sig[bl->out_idx[0]]= bl->state[0];
    }
    else if(bl->op==OP_PLL){
      double Kp=bl->param[0], Ki=bl->param[1], wff=bl->param[2];
      double th=bl->state[0];
      double va=m->sig[bl->in_idx[0]], vb=m->sig[bl->in_idx[1]], vc=m->sig[bl->in_idx[2]];
      double al,be,d,q;
      f_clarke(va,vb,vc,&al,&be);
      f_park(al,be,th,&d,&q);
      m->sig[bl->out_idx[0]]= wff + Kp*q + Ki*bl->state[1];
      m->sig[bl->out_idx[1]]= th;
    }
    else if(bl->out){
      bl->out(bl,m, bl->state);
    }
  }
}

/* =============== Lazo algebraico (Gauss-Seidel) =============== */
static int lazo_algebraico(ModeloC *m){
  for(int iter=0; iter< m->max_iter; iter++){
    double maxdelta=0.0;
    for(int k=0;k<m->n_alg;k++){
      int idx = (int)m->alg_list[k];
      if(idx<0 || idx>=m->n_bloques) continue;
      BloqueC *bl=&m->bloques[idx];
      if(bl->eval_estatico) bl->eval_estatico(bl,m,&maxdelta);
    }
    if(m->n_alg==0){
      double md=0.0;
      for(int i=0;i<m->n_bloques;i++){
        BloqueC *bl=&m->bloques[i];
        if(bl->eval_estatico) bl->eval_estatico(bl,m,&md);
      }
      maxdelta=md;
      if(maxdelta < m->tol) return 0;
      continue;
    }
    if(maxdelta < m->tol) return 0;
  }
  return 1;
}

/* =============== Paso Euler / RK4 =============== */
static void paso_euler(BloqueC *bl, ModeloC *m){
  int n=bl->n_state;
  if(n<=0 || !bl->deriv || !bl->out) return;
  if(n > 128) return;
  double x[128], dx[128];
  for(int i=0;i<n;i++) x[i]= bl->state[i];
  bl->deriv(bl,m,x,dx);
  for(int i=0;i<n;i++) bl->state[i] = x[i] + m->dt * dx[i];
  bl->out(bl,m, bl->state);
}
static void paso_rk4(BloqueC *bl, ModeloC *m){
  int n=bl->n_state;
  if(n<=0 || !bl->deriv || !bl->out) return;
  if(n > 128) return;
  double x0[128], k1[128], k2[128], k3[128], k4[128], xt[128];
  for(int i=0;i<n;i++) x0[i]= bl->state[i];
  bl->deriv(bl,m,x0,k1);
  for(int i=0;i<n;i++) xt[i]= x0[i] + 0.5*m->dt*k1[i];
  bl->deriv(bl,m,xt,k2);
  for(int i=0;i<n;i++) xt[i]= x0[i] + 0.5*m->dt*k2[i];
  bl->deriv(bl,m,xt,k3);
  for(int i=0;i<n;i++) xt[i]= x0[i] + m->dt*k3[i];
  bl->deriv(bl,m,xt,k4);
  for(int i=0;i<n;i++) bl->state[i]= x0[i] + (m->dt/6.0)*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
  bl->out(bl,m, bl->state);
}

/* =============== Actualizar dinamicos =============== */
static void actualizar_dinamicos(ModeloC *m){
  for(int i=0;i<m->n_bloques;i++){
    BloqueC *bl=&m->bloques[i];
    if(bl->Ts > 1e-12){
      if(m->t < bl->t_next_update - 1e-12) continue;
      bl->t_next_update += bl->Ts;
      if(bl->t_next_update < m->t) bl->t_next_update = m->t + bl->Ts;
    }
    if(bl->update){
      bl->update(bl,m);
    } else if(bl->deriv && bl->out){
      if(m->method==0) paso_euler(bl,m);
      else paso_rk4(bl,m);
    }
  }
}

/* =============== Setup block =============== */
static void setup_block(BloqueC *bl, ModeloC *m){
  (void)m;
  bl->init = NULL;
  bl->eval_estatico = NULL;
  bl->deriv = NULL;
  bl->out = NULL;
  bl->update = NULL;
  bl->t_next_update = 0.0;
  switch(bl->op){
    case OP_SRC_CONST:
    case OP_SRC_STEP:
    case OP_SRC_RAMP:
    case OP_SRC_SIN:
    case OP_SRC_TRIF:
    case OP_PULSO_RECT:
    case OP_SRC_CSV:
    case OP_SRC_TABLE:
      break;
    case OP_GAIN:
    case OP_SUM:
    case OP_CLARKE:
    case OP_INV_CLARKE:
    case OP_PARK:
    case OP_INV_PARK:
    case OP_SATURAR:
    case OP_TRANSFORMADOR:
    case OP_PANEL_SOLAR:
    case OP_MEDIDOR_POTENCIA:
    case OP_INTERRUPTOR:
    case OP_PUENTE_INV_3F:
    case OP_PUENTE_INV_1F:
    case OP_MUX:
    case OP_DEMUX:
    case OP_LUT1D:
    case OP_LUT2D:
    case OP_LUT3D:
    case OP_LOGICO:
    case OP_RELACIONAL:
      bl->eval_estatico = eval_estatico;
      break;
    case OP_PWM_1F:
      bl->eval_estatico = pwm_1f_eval;
      break;
    case OP_PWM_SPWM:
      bl->eval_estatico = pwm_spwm_eval;
      break;
    case OP_PWM_SVPWM:
      bl->eval_estatico = pwm_svpwm_eval;
      break;
    case OP_HW_SERIAL:
      bl->eval_estatico = hw_serial_eval;
      break;
    case OP_RES_TERMICA:
    case OP_ENGRANAJE:
    case OP_EMBRAGUE:
    case OP_FALLO_PROG:
    case OP_FALLO_EVENTO:
    case OP_MULTIPLICADOR:
    case OP_SAT_VECTORIAL:
    case OP_QD:
    case OP_CARGA_PQ_3F:
    case OP_CARGA_PQ_1F:
    case OP_RESISTENCIA:
    case OP_CALCULO_IDC:
      bl->eval_estatico = eval_estatico;
      break;
    case OP_INTEGRADOR:
      bl->deriv = integrador_deriv;
      bl->out = integrador_out;
      break;
    case OP_TF:
      bl->update = update_tf;
      break;
    case OP_PID:
      bl->update = update_pid;
      break;
    case OP_RELAY:
      bl->update = update_relay;
      break;
    case OP_DIODO:
      bl->update = update_diodo;
      break;
    case OP_RETENEDOR:
      bl->update = update_retenedor;
      break;
    case OP_MAQ_ESTADOS:
      bl->update = update_maq_estados;
      break;
    case OP_PLL:
      bl->update = update_pll;
      break;
    case OP_MAQ_PMAC:
      bl->deriv = maq_pmac_deriv;
      bl->out = maq_pmac_out;
      break;
    case OP_MAQ_INDUCCION:
      bl->deriv = maq_ind_deriv;
      bl->out = maq_ind_out;
      break;
    case OP_MAQ_SINCRONA:
      bl->deriv = maq_sinc_deriv;
      bl->out = maq_sinc_out;
      break;
    case OP_MAQ_CC:
      bl->deriv = maq_cc_deriv;
      bl->out = maq_cc_out;
      break;
    case OP_MAQ_DC_PM:
      bl->deriv = maq_dcpm_deriv;
      bl->out = maq_dcpm_out;
      break;
    case OP_POT_BUCK:
    case OP_POT_BOOST:
    case OP_POT_BUCKBOOST:
      bl->deriv = pot_dcdc_deriv;
      bl->out = pot_dcdc_out;
      break;
    case OP_POT_RECT_3F:
      bl->deriv = pot_rect_deriv;
      bl->out = pot_rect_out;
      break;
    case OP_POT_INV_3F:
      bl->deriv = pot_inv_deriv;
      bl->out = pot_inv_out;
      break;
    case OP_POT_INV_1F:
      bl->deriv = pot_inv1f_deriv;
      bl->out = pot_inv1f_out;
      break;
    case OP_CARGA_RL_3F:
      bl->deriv = carga_rl_deriv;
      bl->out = carga_rl_out;
      break;
    case OP_EJE_MECANICO:
      bl->deriv = eje_mecanico_deriv;
      bl->out = eje_mecanico_out;
      break;
    case OP_BATERIA:
      bl->deriv = bateria_deriv;
      bl->out = bateria_out;
      break;
    case OP_BATERIA_ECM:
      bl->deriv = bateria_ecm_deriv;
      bl->out = bateria_ecm_out;
      break;
    case OP_MASA_TERMICA:
      bl->deriv = masa_termica_deriv;
      bl->out = masa_termica_out;
      break;
    case OP_EJE_FLEXIBLE:
      bl->deriv = eje_flexible_deriv;
      bl->out = eje_flexible_out;
      break;
    case OP_VEHICULO:
      bl->deriv = vehiculo_deriv;
      bl->out = vehiculo_out;
      break;
    case OP_INDUCTOR:
      bl->deriv = inductor_deriv;
      bl->out = inductor_out;
      break;
    case OP_CAPACITOR:
      bl->deriv = capacitor_deriv;
      bl->out = capacitor_out;
      break;
    case OP_LIM_RAPIDEZ:
      bl->deriv = limrap_deriv;
      bl->out = limrap_out;
      break;
    case OP_VCVS:
    case OP_VCCS:
      bl->eval_estatico = eval_estatico; break;
    case OP_MUTUAL_INDUCTOR:
      bl->deriv = mutual_deriv; bl->out = mutual_out; break;
    case OP_MNA:
      bl->update = update_mna;
      break;
    default:
      break;
  }
}

/* =============== API =============== */
API int m_sim_iniciar(ModeloC *m){
  for(int i=0;i<m->n_bloques;i++) setup_block(&m->bloques[i], m);
  m->t = 0.0;
  m->error_flag = 0;
  for(int i=0;i<m->n_sig;i++) m->sig[i]=0.0;
  for(int i=0;i<m->n_bloques;i++) m->bloques[i].t_next_update = 0.0;
  actualizar_fuentes(m);
  emit_inicial(m);
  int err = lazo_algebraico(m);
  m->error_flag = err;
  return err;
}
API int m_sim_paso(ModeloC *m){
  m->t += m->dt;
  actualizar_fuentes(m);
  int err = lazo_algebraico(m);
  if(err){ m->error_flag=1; return 1; }
  actualizar_dinamicos(m);
  err = lazo_algebraico(m);
  m->error_flag = err;
  return err;
}
API void m_sim_run(ModeloC *m, int n_steps, int n_rec, long long *rec_idx, double *rec_buf){
  m_sim_iniciar(m);
  for(int s=0; s<n_steps; s++){
    if(s>0){
      if(m_sim_paso(m)) break;
    }
    for(int r=0; r<n_rec; r++){
      int idx = (int)rec_idx[r];
      double v = 0.0;
      if(idx>=0 && idx < m->n_sig) v = m->sig[idx];
      rec_buf[r*n_steps + s] = v;
    }
  }
}
API int m_sim_guardar(ModeloC *m, double *buf){
  int p=0;
  buf[p++]= m->t;
  for(int i=0;i<m->n_sig;i++) buf[p++]= m->sig[i];
  for(int i=0;i<m->n_bloques;i++){
    BloqueC *bl=&m->bloques[i];
    for(int j=0;j<bl->n_state;j++) buf[p++]= bl->state[j];
    for(int j=0;j<bl->n_ws;j++) buf[p++]= bl->ws[j];
    buf[p++]= bl->t_next_update;
  }
  return p;
}
API void m_sim_restaurar(ModeloC *m, const double *buf){
  int p=0;
  m->t = buf[p++];
  for(int i=0;i<m->n_sig;i++) m->sig[i]= buf[p++];
  for(int i=0;i<m->n_bloques;i++){
    BloqueC *bl=&m->bloques[i];
    for(int j=0;j<bl->n_state;j++) bl->state[j]= buf[p++];
    for(int j=0;j<bl->n_ws;j++) bl->ws[j]= buf[p++];
    bl->t_next_update = buf[p++];
  }
}
API int m_hil_ws_size(void){
  return HIL_WS_DOUBLES;
}
#ifdef _WIN32
API void m_hw_serial_cerrar(BloqueC *bl){
  if(!bl) return;
  if(bl->op != OP_HW_SERIAL) return;
  if(bl->n_ws < HIL_WS_DOUBLES) return;
  SerialHIL *s = (SerialHIL*)(bl->ws + 2);
  if(s->connected) serial_close(s);
  bl->ws[0]=0;
}
#endif
