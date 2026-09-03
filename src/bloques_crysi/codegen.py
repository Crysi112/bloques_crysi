from __future__ import annotations

import io
from pathlib import Path
from typing import List, Optional, TYPE_CHECKING

from . import opcodes as ops

if TYPE_CHECKING:
    from .modelo import Modelo

def _fmt_f(v: float) -> str:
    s = repr(float(v))
    if "e" in s or "." in s:
        return s
    return s + ".0"

def _idx(idxs: list, k: int) -> str:
    return str(int(idxs[k]))

_HEADER_MATH = """\
#include <math.h>
#include <string.h>
#include <stdio.h>
#ifndef M_PI
#define M_PI 3.141592653589793
#endif
#define MNA_MAX_NX 128

static inline void _clarke(double va, double vb, double vc,
                            double *al, double *be) {
    *al = (2.0/3.0) * (va - 0.5*vb - 0.5*vc);
    *be = (1.0/3.0) * (vb - vc) * 1.7320508075688772;
}
static inline void _inv_clarke(double al, double be,
                                double *va, double *vb, double *vc) {
    *va = al;
    *vb = -0.5*al + 0.8660254037844386*be;
    *vc = -0.5*al - 0.8660254037844386*be;
}
static inline void _park(double al, double be, double th,
                          double *d, double *q) {
    double c = cos(th), s = sin(th);
    *d =  al*c + be*s;
    *q = -al*s + be*c;
}
static inline void _inv_park(double d, double q, double th,
                              double *al, double *be) {
    double c = cos(th), s = sin(th);
    *al = d*c - q*s;
    *be = d*s + q*c;
}
static inline void _pll_step(double *st, const double *p, double dt,
                              double va, double vb, double vc,
                              double *out_w, double *out_th) {
    double Kp=p[0], Ki=p[1], w_ff=p[2];
    double al, be, vd, vq;
    _clarke(va, vb, vc, &al, &be);
    _park(al, be, st[0], &vd, &vq);
    double err = vq;
    st[1] += dt * Ki * err;
    double w = w_ff + Kp*err + st[1];
    st[0] += dt * w;
    *out_w = w; *out_th = st[0];
}
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
"""

_PLL_HELPER = """"""


_MNA_HEADER = """#define MNA_MAX_NX 128

#ifdef _MSC_VER
  #define THREAD_LOCAL __declspec(thread)
#elif defined(__STDC_NO_THREADS__) || (defined(__MINGW32__) && __GNUC__ < 8)
  #define THREAD_LOCAL_FALLBACK_CALLOC 1
#else
  #define THREAD_LOCAL __thread
#endif

#ifndef THREAD_LOCAL_FALLBACK_CALLOC
  static THREAD_LOCAL double mna_G[128][128];
  static THREAD_LOCAL double mna_Cmat[128][128];
  static THREAD_LOCAL double mna_W[128][128];
  static THREAD_LOCAL double mna_Ap[128][128];
  static THREAD_LOCAL double mna_b[128];
  static THREAD_LOCAL double mna_b_diode[128];
#endif
"""


class GeneradorCodigo:
    def __init__(self, modelo: "Modelo") -> None:
        self._m = modelo
        modelo._resolver()
        self._bloques = [b for b in modelo.bloques
                         if not getattr(b, "es_scope", False)]
        self._n_sig = modelo._nsig
        self._dt = modelo.dt
        self._metodo = modelo.metodo
        self._max_iter = modelo.max_iter
        self._tol = modelo.tol
        self._w_opt = modelo.w_opt
        self._orden_alg = modelo._orden_estatico()

    def generar(self, ruta: str, nombre: str = "modelo") -> str:
        codigo = self._construir(nombre)
        Path(ruta).write_text(codigo, encoding="utf-8")
        return ruta

    def codigo(self, nombre: str = "modelo") -> str:
        return self._construir(nombre)

    def _construir(self, nombre: str) -> str:
        buf = io.StringIO()
        w = buf.write
        w(self._gen_header(nombre))
        w(_HEADER_MATH)
        if any(b.op == ops.OP_PLL for b in self._bloques):
            w(_PLL_HELPER)
        w(_MNA_HEADER)
        w(self._gen_signals())
        w(self._gen_block_data())
        w(f"\nstatic double _t = 0.0;\n")
        w(f"static const double _DT       = {_fmt_f(self._dt)};\n")
        w(f"static const double _TOL      = {_fmt_f(self._tol)};\n")
        w(f"static const int    _MAX_ITER = {self._max_iter};\n")
        w(f"static const double _W_OPT    = {_fmt_f(self._w_opt)};\n\n")
        w(self._gen_actualizar_fuentes())
        w(self._gen_lazo_algebraico())
        w(self._gen_actualizar_dinamicos())
        w(self._gen_api(nombre))
        return buf.getvalue()

    def _gen_header(self, nombre: str) -> str:
        bl = ", ".join(b.nombre for b in self._bloques[:6])
        if len(self._bloques) > 6:
            bl += f", ... ({len(self._bloques)} total)"
        return (
            f"/*\n"
            f" * AUTO-GENERADO por bloques_crysi codegen\n"
            f" * Modelo : {nombre}\n"
            f" * Bloques: {bl}\n"
            f" * dt     : {self._dt}\n"
            f" * Metodo : {'RK4' if self._metodo else 'Euler'}\n"
            f" * Senales: {self._n_sig}\n"
            f" *\n"
            f" * gcc -O2 -o sim {nombre}.c -lm\n"
            f" * Embebido: enlaza solo modelo_init/paso/run en tu firmware\n"
            f" */\n\n"
        )

    def _gen_signals(self) -> str:
        lines = [f"/* === Senales ({self._n_sig} canales) === */",
                 f"static double _sig[{self._n_sig}];"]
        for b in self._bloques:
            lbl = getattr(b, "NOMBRES", None)
            for k, idx in enumerate(b.out_idx):
                canal = lbl[k] if lbl and k < len(lbl) else f"{b.nombre}[{k}]"
                lines.append(f"/* _sig[{idx}] = {canal} */")
        return "\n".join(lines) + "\n\n"

    def _gen_block_data(self) -> str:
        lines = ["/* === Parametros, estados y workspace por bloque === */"]
        for i, b in enumerate(self._bloques):
            lines.append(f"\n/* bloque {i}: {b.nombre} ({b.etiqueta}) */")
            for j, p in enumerate(b.param):
                lines.append(f"#define _B{i}P{j} {_fmt_f(p)}")
            if b.op in (ops.OP_MNA, ops.OP_SRC_TABLE, ops.OP_SRC_CSV, ops.OP_BATERIA_ECM):
                vals = ", ".join(_fmt_f(v) for v in b.param)
                lines.append(f"static double _b{i}_param[] = {{{vals}}};")
                lines.append(f"static const int _b{i}_nparam = {len(b.param)};")
            if b.n_state > 0:
                inits = ", ".join(_fmt_f(v) for v in b.estados_iniciales)
                lines.append(f"static double _b{i}_x[{b.n_state}] = {{{inits}}};")
                if b.op == ops.OP_MNA:
                    n_ws = max(5 * b.n_state, 200)
                else:
                    n_ws = 5 * b.n_state
                lines.append(f"static double _b{i}_ws[{max(n_ws, 1)}];")
        return "\n".join(lines) + "\n\n"

    def _gen_actualizar_fuentes(self) -> str:
        lines = ["/* === Fuentes === */",
                 "static void _actualizar_fuentes(void) {"]
        for i, b in enumerate(self._bloques):
            code = self._gen_fuente(i, b)
            if code:
                lines.append(f"    /* {b.nombre} */")
                lines.extend("    " + ln for ln in code.splitlines())
        lines.append("}\n")
        return "\n".join(lines) + "\n"

    def _gen_fuente(self, i: int, b) -> Optional[str]:
        op = b.op
        o0 = _idx(b.out_idx, 0)
        if op == ops.OP_SRC_CONST:
            return f"_sig[{o0}] = _B{i}P0;"
        if op == ops.OP_SRC_STEP:
            return f"_sig[{o0}] = (_t >= _B{i}P1) ? _B{i}P0 : _B{i}P2;"
        if op == ops.OP_SRC_RAMP:
            return (f"_sig[{o0}] = (_t < _B{i}P1) ? _B{i}P2 "
                    f": _B{i}P2 + _B{i}P0*(_t - _B{i}P1);")
        if op == ops.OP_SRC_SIN:
            return (f"_sig[{o0}] = _B{i}P0*sin(2.0*M_PI*_B{i}P1*_t + _B{i}P2)"
                    f" + _B{i}P3;")
        if op == ops.OP_SRC_TRIF:
            o1 = _idx(b.out_idx, 1)
            o2 = _idx(b.out_idx, 2)
            return (
                f"_sig[{o0}] = _B{i}P0*sin(2.0*M_PI*_B{i}P1*_t + _B{i}P2);\n"
                f"_sig[{o1}] = _B{i}P0*sin(2.0*M_PI*_B{i}P1*_t + _B{i}P2 - 2.0943951023931953);\n"
                f"_sig[{o2}] = _B{i}P0*sin(2.0*M_PI*_B{i}P1*_t + _B{i}P2 - 4.1887902047863905);"
            )
        if op == ops.OP_PULSO_RECT:
            return (f"{{ double _ph=fmod(_t+_B{i}P3,_B{i}P1); "
                    f"_sig[{o0}]=_B{i}P4+((_ph<_B{i}P2*_B{i}P1)?_B{i}P0:0.0); }}")
        if op in (ops.OP_SRC_TABLE, ops.OP_SRC_CSV):
            return (f"{{ int n=(int)_b{i}_param[0]; int interp=_b{i}_param[1]>0.5;"
                    f" const double *tt=&_b{i}_param[2]; const double *yy=&_b{i}_param[2+n];"
                    f" double y; if(_t<=tt[0]) y=yy[0]; else if(_t>=tt[n-1]) y=yy[n-1]; else {{ "
                    f" int lo=0, hi=n-1; while(hi-lo>1){{ int mid=(lo+hi)/2; if(tt[mid]<=_t)lo=mid; else hi=mid; }} "
                    f" if(interp){{ double f=(_t-tt[lo])/(tt[lo+1]-tt[lo]); y=yy[lo]+f*(yy[lo+1]-yy[lo]); }} "
                    f" else y=yy[lo]; }} _sig[{o0}]=y; }}")
        return None

    def _gen_lazo_algebraico(self) -> str:
        lines = [
            "/* === Lazo algebraico === */",
            "static void _lazo_algebraico(void) {",
            "    double _md, _old, _nv;",
            "    double _w = _W_OPT;",
            "    int _it;",
            "    for (_it = 0; _it < _MAX_ITER; _it++) {",
            "        _md = 0.0;",
        ]
        for idx in self._orden_alg:
            b = self._bloques[idx]
            code = self._gen_estatico(idx, b)
            if code:
                lines.append(f"        /* {b.nombre} */")
                lines.extend("        " + ln for ln in code.splitlines())
        lines += [
            "        if (_md < _TOL) break;",
            "    }",
            "}\n",
        ]
        return "\n".join(lines) + "\n"

    def _relax(self, out_k: int, expr: str) -> str:
        o = f"_sig[{out_k}]"
        return (f"_old={o}; _nv=_old+_w*({expr}-_old); {o}=_nv;\n"
                f"if(fabs(_nv-_old)>_md) _md=fabs(_nv-_old);")

    def _gen_estatico(self, i: int, b) -> Optional[str]:
        op = b.op

        if op == ops.OP_GAIN:
            return self._relax(int(b.out_idx[0]),
                               f"_B{i}P0*_sig[{_idx(b.in_idx,0)}]")

        if op == ops.OP_SUM:
            terms = "+".join(
                f"_B{i}P{j}*_sig[{_idx(b.in_idx,j)}]" for j in range(b.n_in))
            return self._relax(int(b.out_idx[0]), terms)

        if op == ops.OP_SATURAR:
            u = f"_sig[{_idx(b.in_idx,0)}]"
            expr = f"(({u}<_B{i}P0)?_B{i}P0:(({u}>_B{i}P1)?_B{i}P1:{u}))"
            return self._relax(int(b.out_idx[0]), expr)

        if op == ops.OP_MULTIPLICADOR:
            return self._relax(int(b.out_idx[0]),
                               f"_sig[{_idx(b.in_idx,0)}]*_sig[{_idx(b.in_idx,1)}]")

        if op == ops.OP_CLARKE:
            i0,i1,i2 = _idx(b.in_idx,0),_idx(b.in_idx,1),_idx(b.in_idx,2)
            o0,o1 = int(b.out_idx[0]),int(b.out_idx[1])
            return (f"{{ double _al,_be; _clarke(_sig[{i0}],_sig[{i1}],_sig[{i2}],&_al,&_be);\n"
                    f"  {self._relax(o0,'_al')}\n  {self._relax(o1,'_be')} }}")

        if op == ops.OP_INV_CLARKE:
            i0,i1 = _idx(b.in_idx,0),_idx(b.in_idx,1)
            o0,o1,o2 = int(b.out_idx[0]),int(b.out_idx[1]),int(b.out_idx[2])
            return (f"{{ double _va,_vb,_vc; _inv_clarke(_sig[{i0}],_sig[{i1}],&_va,&_vb,&_vc);\n"
                    f"  {self._relax(o0,'_va')}\n  {self._relax(o1,'_vb')}\n  {self._relax(o2,'_vc')} }}")

        if op == ops.OP_PARK:
            i0,i1,i2 = _idx(b.in_idx,0),_idx(b.in_idx,1),_idx(b.in_idx,2)
            o0,o1 = int(b.out_idx[0]),int(b.out_idx[1])
            return (f"{{ double _d,_q; _park(_sig[{i0}],_sig[{i1}],_sig[{i2}],&_d,&_q);\n"
                    f"  {self._relax(o0,'_d')}\n  {self._relax(o1,'_q')} }}")

        if op == ops.OP_INV_PARK:
            i0,i1,i2 = _idx(b.in_idx,0),_idx(b.in_idx,1),_idx(b.in_idx,2)
            o0,o1 = int(b.out_idx[0]),int(b.out_idx[1])
            return (f"{{ double _alp,_bep; _inv_park(_sig[{i0}],_sig[{i1}],_sig[{i2}],&_alp,&_bep);\n"
                    f"  {self._relax(o0,'_alp')}\n  {self._relax(o1,'_bep')} }}")

        if op == ops.OP_INTERRUPTOR:
            ic_,iva,ivc = _idx(b.in_idx,0),_idx(b.in_idx,1),_idx(b.in_idx,2)
            o0,o1 = int(b.out_idx[0]),int(b.out_idx[1])
            return (f"{{ double _r=_B{i}P1; if(_sig[{ic_}]>0.5) _r=_B{i}P0;\n"
                    f"  double _v=_sig[{iva}]-_sig[{ivc}];\n"
                    f"  {self._relax(o0,'_v/_r')}\n  {self._relax(o1,'_v')} }}")

        if op == ops.OP_PUENTE_INV_3F:
            i0 = _idx(b.in_idx,0)
            ia,ib_,ic_ = _idx(b.in_idx,1),_idx(b.in_idx,2),_idx(b.in_idx,3)
            o0,o1,o2 = int(b.out_idx[0]),int(b.out_idx[1]),int(b.out_idx[2])
            return (f"{{ int _prom=(_B{i}P0>0.5); double _vdc=_sig[{i0}];\n"
                    f"  double _a=_sig[{ia}],_b2=_sig[{ib_}],_c2=_sig[{ic_}];\n"
                    f"  double _va2,_vb2,_vc2;\n"
                    f"  if(_prom){{ _va2=_vdc*(_a-0.5*_b2-0.5*_c2)/3.0;"
                    f" _vb2=_vdc*(_b2-0.5*_a-0.5*_c2)/3.0;"
                    f" _vc2=_vdc*(_c2-0.5*_a-0.5*_b2)/3.0; }}\n"
                    f"  else{{ int _sa=(_a>0.5),_sb=(_b2>0.5),_sc=(_c2>0.5);\n"
                    f"    _va2=_vdc*(2*_sa-_sb-_sc)/3.0;"
                    f" _vb2=_vdc*(2*_sb-_sa-_sc)/3.0;"
                    f" _vc2=_vdc*(2*_sc-_sa-_sb)/3.0; }}\n"
                    f"  {self._relax(o0,'_va2')}\n  {self._relax(o1,'_vb2')}\n  {self._relax(o2,'_vc2')} }}")

        if op == ops.OP_MEDIDOR_POTENCIA:
            fases = (b.n_in - 2) // 2
            sumas = "+".join(
                f"_sig[{_idx(b.in_idx,k)}]*_sig[{_idx(b.in_idx,fases+k)}]"
                for k in range(fases))
            o0 = int(b.out_idx[0])
            lines = [f"{{ double _pe={sumas}; {self._relax(o0,'_pe')}"]
            if b.n_out == 3:
                i0,i1,i2 = _idx(b.in_idx,0),_idx(b.in_idx,1),_idx(b.in_idx,2)
                i3,i4,i5 = _idx(b.in_idx,3),_idx(b.in_idx,4),_idx(b.in_idx,5)
                o1 = int(b.out_idx[1])
                lines.append(
                    f"  double _qe=((_sig[{i0}]-_sig[{i1}])*_sig[{i5}]+"
                    f"(_sig[{i1}]-_sig[{i2}])*_sig[{i3}]+"
                    f"(_sig[{i2}]-_sig[{i0}])*_sig[{i4}])/1.7320508075688772;"
                )
                lines.append(f"  {self._relax(o1,'_qe')}")
            im = 2*fases
            ite = int(b.in_idx[im]) if b.in_idx[im] >= 0 else -1
            iwm = int(b.in_idx[im+1]) if b.in_idx[im+1] >= 0 else -1
            om = int(b.out_idx[b.n_out-1])
            pm_expr = f"_sig[{ite}]*_sig[{iwm}]" if ite >= 0 and iwm >= 0 else "0.0"
            lines.append(f"  {self._relax(om, pm_expr)} }}")
            return "\n".join(lines)

        if op == ops.OP_CARGA_PQ_3F:
            i0,i1,i2 = _idx(b.in_idx,0),_idx(b.in_idx,1),_idx(b.in_idx,2)
            o0,o1,o2 = int(b.out_idx[0]),int(b.out_idx[1]),int(b.out_idx[2])
            return (
                f"{{ double _va3=_sig[{i0}],_vb3=_sig[{i1}],_vc3=_sig[{i2}];\n"
                f"  double _vn3=(_va3+_vb3+_vc3)/3.0;\n"
                f"  double _v03=_va3-_vn3,_v13=_vb3-_vn3,_v23=_vc3-_vn3;\n"
                f"  double _S=0.8164965809277260,_I2=0.7071067811865475,_I6=0.4082482904638630;\n"
                f"  double _al3=_S*(_v03-0.5*_v13-0.5*_v23),_be3=_I2*(_v13-_v23);\n"
                f"  double _vm2=_al3*_al3+_be3*_be3;\n"
                f"  double _ia3=0.0,_ib3=0.0,_ic3=0.0;\n"
                f"  if(_vm2>1e-6){{ double _ial=(_B{i}P0*_al3+_B{i}P1*_be3)/_vm2;\n"
                f"    double _ibe=(_B{i}P0*_be3-_B{i}P1*_al3)/_vm2;\n"
                f"    _ia3=_S*_ial; _ib3=-_I6*_ial+_I2*_ibe; _ic3=-_I6*_ial-_I2*_ibe; }}\n"
                f"  {self._relax(o0,'_ia3')}\n  {self._relax(o1,'_ib3')}\n  {self._relax(o2,'_ic3')} }}"
            )

        return f"/* AVISO: opcode {b.op} ({b.etiqueta}) sin soporte codegen — sin cambio */"

    def _gen_actualizar_dinamicos(self) -> str:
        lines = ["/* === Actualizacion dinamica === */",
                 "static void _actualizar_dinamicos(void) {"]
        for i, b in enumerate(self._bloques):
            code = self._gen_dinamico(i, b)
            if code:
                lines.append(f"    /* {b.nombre} ({b.etiqueta}) */")
                lines.extend("    " + ln for ln in code.splitlines())
        lines.append("}\n")
        return "\n".join(lines) + "\n"

    def _gen_dinamico(self, i: int, b) -> Optional[str]:
        op = b.op
        if op not in ops.ES_DINAMICO:
            return None
        xs = f"_b{i}_x"
        dt = "_DT"

        if op == ops.OP_INTEGRADOR:
            in0 = _idx(b.in_idx, 0)
            o0  = int(b.out_idx[0])
            return f"{xs}[0] += {dt}*_sig[{in0}]; _sig[{o0}]={xs}[0];"

        if op == ops.OP_PID:
            in0 = _idx(b.in_idx, 0)
            o0  = int(b.out_idx[0])
            return (
                f"{{ double _e=_sig[{in0}];\n"
                f"  double _ud=(_B{i}P2*(_e-{xs}[1])+_B{i}P3*{xs}[2])/(_B{i}P3+{dt});\n"
                f"  double _u0=_B{i}P0*_e+_B{i}P1*{xs}[0]+_ud;\n"
                f"  if(!(_u0>_B{i}P5&&_e>0)&&!(_u0<_B{i}P4&&_e<0)) {xs}[0]+={dt}*_e;\n"
                f"  {xs}[1]=_e; {xs}[2]=_ud;\n"
                f"  if(_u0>_B{i}P5) _u0=_B{i}P5;\n"
                f"  if(_u0<_B{i}P4) _u0=_B{i}P4;\n"
                f"  _sig[{o0}]=_u0; }}"
            )

        if op == ops.OP_PLL:
            ia,ib_,ic_ = _idx(b.in_idx,0),_idx(b.in_idx,1),_idx(b.in_idx,2)
            ow,oth = int(b.out_idx[0]),int(b.out_idx[1])
            return (
                f"{{ double _pp[4]={{_B{i}P0,_B{i}P1,_B{i}P2,_B{i}P3}};\n"
                f"  double _sw,_sth;\n"
                f"  _pll_step({xs},_pp,{dt},_sig[{ia}],_sig[{ib_}],_sig[{ic_}],&_sw,&_sth);\n"
                f"  _sig[{ow}]=_sw; _sig[{oth}]=_sth; }}"
            )

        if op in (ops.OP_POT_BUCK, ops.OP_POT_BOOST, ops.OP_POT_BUCKBOOST):
            vin = _idx(b.in_idx, 0)
            d_  = _idx(b.in_idx, 1)
            o0, o1 = int(b.out_idx[0]), int(b.out_idx[1])
            if op == ops.OP_POT_BUCK:
                diL = f"(_sig[{d_}]*_sig[{vin}]-{xs}[1])/_B{i}P0"
                dvC = f"({xs}[0]-{xs}[1]/_B{i}P2)/_B{i}P1"
            elif op == ops.OP_POT_BOOST:
                diL = f"(_sig[{vin}]-(1.0-_sig[{d_}])*{xs}[1])/_B{i}P0"
                dvC = f"((1.0-_sig[{d_}])*{xs}[0]-{xs}[1]/_B{i}P2)/_B{i}P1"
            else:
                diL = f"(_sig[{d_}]*_sig[{vin}]-(1.0-_sig[{d_}])*{xs}[1])/_B{i}P0"
                dvC = f"((1.0-_sig[{d_}])*{xs}[0]-{xs}[1]/_B{i}P2)/_B{i}P1"
            return (
                f"{xs}[0]+={dt}*({diL});\n"
                f"{xs}[1]+={dt}*({dvC});\n"
                f"_sig[{o0}]={xs}[1]; _sig[{o1}]={xs}[0];"
            )

        if op == ops.OP_POT_RECT_3F:
            i0,i1,i2 = _idx(b.in_idx,0),_idx(b.in_idx,1),_idx(b.in_idx,2)
            o0,o1 = int(b.out_idx[0]),int(b.out_idx[1])
            return (
                f"{{ double _vmx=fmax(fmax(_sig[{i0}],_sig[{i1}]),_sig[{i2}]);\n"
                f"  double _vmn=fmin(fmin(_sig[{i0}],_sig[{i1}]),_sig[{i2}]);\n"
                f"  double _vrec=_vmx-_vmn, _vC={xs}[0];\n"
                f"  double _ich=(_vrec-_vC)/_B{i}P2; if(_ich<0.0)_ich=0.0;\n"
                f"  {xs}[0]+={dt}*(_ich-_vC/_B{i}P1)/_B{i}P0;\n"
                f"  _sig[{o0}]={xs}[0]; _sig[{o1}]=_ich; }}"
            )

        if op == ops.OP_CARGA_RL_3F:
            i0,i1,i2 = _idx(b.in_idx,0),_idx(b.in_idx,1),_idx(b.in_idx,2)
            o0,o1,o2 = int(b.out_idx[0]),int(b.out_idx[1]),int(b.out_idx[2])
            return (
                f"{{ double _vn=(_sig[{i0}]+_sig[{i1}]+_sig[{i2}])/3.0;\n"
                f"  {xs}[0]+={dt}*(_sig[{i0}]-_vn-_B{i}P0*{xs}[0])/_B{i}P1;\n"
                f"  {xs}[1]+={dt}*(_sig[{i1}]-_vn-_B{i}P0*{xs}[1])/_B{i}P1;\n"
                f"  _sig[{o0}]={xs}[0]; _sig[{o1}]={xs}[1]; _sig[{o2}]=-({xs}[0]+{xs}[1]); }}"
            )

        if op == ops.OP_EJE_MECANICO:
            o0, o1 = int(b.out_idx[0]), int(b.out_idx[1])
            sumas = "+".join(f"_sig[{_idx(b.in_idx,k)}]" for k in range(b.n_in-1))
            tl = _idx(b.in_idx, b.n_in-1)
            return (
                f"{{ double _sTe={sumas};\n"
                f"  {xs}[0]+={dt}*(_sTe-_sig[{tl}]-_B{i}P1*{xs}[0])/_B{i}P0;\n"
                f"  {xs}[1]+={dt}*{xs}[0];\n"
                f"  _sig[{o0}]={xs}[0]; _sig[{o1}]={xs}[1]; }}"
            )

        if op == ops.OP_MAQ_PMAC:
            ia,ib_,ic_ = _idx(b.in_idx,0),_idx(b.in_idx,1),_idx(b.in_idx,2)
            tl = _idx(b.in_idx, 3)
            o = [int(b.out_idx[k]) for k in range(9)]
            return (
                f"{{ double _rs=_B{i}P0,_Ld=_B{i}P1,_Lq=_B{i}P2,_lam=_B{i}P3;\n"
                f"  double _P=_B{i}P4,_J=_B{i}P5,_Bm=_B{i}P6;\n"
                f"  double _we=(_P/2.0)*{xs}[2];\n"
                f"  double _al,_be,_vd,_vq;\n"
                f"  _clarke(_sig[{ia}],_sig[{ib_}],_sig[{ic_}],&_al,&_be);\n"
                f"  _park(_al,_be,{xs}[3],&_vd,&_vq);\n"
                f"  double _ld={xs}[1]*_Ld+_lam;\n"
                f"  double _Te=1.5*(_P/2.0)*(_ld*{xs}[0]-_Lq*{xs}[0]*{xs}[1]);\n"
                f"  double _d0=(_vq-_rs*{xs}[0]-_we*_ld)/_Lq;\n"
                f"  double _d1=(_vd-_rs*{xs}[1]+_we*_Lq*{xs}[0])/_Ld;\n"
                f"  double _d2=(_Te-_sig[{tl}]-_Bm*{xs}[2])/_J;\n"
                f"  {xs}[0]+={dt}*_d0; {xs}[1]+={dt}*_d1;\n"
                f"  {xs}[2]+={dt}*_d2; {xs}[3]+={dt}*_we;\n"
                f"  double _a2,_b2,_ia2,_ib2,_ic2;\n"
                f"  _inv_park({xs}[1],{xs}[0],{xs}[3],&_a2,&_b2);\n"
                f"  _inv_clarke(_a2,_b2,&_ia2,&_ib2,&_ic2);\n"
                f"  _sig[{o[0]}]=_ia2; _sig[{o[1]}]=_ib2; _sig[{o[2]}]=_ic2;\n"
                f"  _sig[{o[3]}]={xs}[0]; _sig[{o[4]}]={xs}[1];\n"
                f"  _sig[{o[5]}]={xs}[2]; _sig[{o[6]}]={xs}[3]/(_P/2.0);\n"
                f"  _sig[{o[7]}]={xs}[3]; _sig[{o[8]}]=_Te; }}"
            )

        if op == ops.OP_MAQ_DC_PM:
            iva = _idx(b.in_idx, 0)
            tl  = _idx(b.in_idx, 1)
            o = [int(b.out_idx[k]) for k in range(4)]
            return (
                f"{{ double _ra=_B{i}P0,_La=_B{i}P1,_Kt=_B{i}P2,_J=_B{i}P3,_Bm=_B{i}P4;\n"
                f"  {xs}[0]+={dt}*(_sig[{iva}]-_ra*{xs}[0]-_Kt*{xs}[1])/_La;\n"
                f"  {xs}[1]+={dt}*(_Kt*{xs}[0]-_sig[{tl}]-_Bm*{xs}[1])/_J;\n"
                f"  {xs}[2]+={dt}*{xs}[1];\n"
                f"  _sig[{o[0]}]={xs}[0]; _sig[{o[1]}]={xs}[1];\n"
                f"  _sig[{o[2]}]={xs}[2]; _sig[{o[3]}]=_Kt*{xs}[0]; }}"
            )

        if op == ops.OP_BATERIA:
            in0 = _idx(b.in_idx, 0)
            o0, o1 = int(b.out_idx[0]), int(b.out_idx[1])
            return (
                f"{{ double _i=_sig[{in0}],_Q=_B{i}P2,_tau=_B{i}P6;\n"
                f"  double _di=_i/3600.0;\n"
                f"  if(_di>0.0&&{xs}[0]>=0.9*_Q) _di=0.0;\n"
                f"  if(_di<0.0&&{xs}[0]<=0.1*_Q) _di=0.0;\n"
                f"  {xs}[0]+={dt}*_di; {xs}[1]+={dt}*(_i-{xs}[1])/_tau;\n"
                f"  double _it={xs}[0];\n"
                f"  if(_it<0.0)_it=0.0; if(_it>0.9*_Q)_it=0.9*_Q;\n"
                f"  double _d1=_Q-_it; if(_d1<1e-9*_Q)_d1=1e-9*_Q;\n"
                f"  double _ev=_B{i}P3*exp(-_B{i}P4*_it);\n"
                f"  double _E=({xs}[1]>=0.0)?\n"
                f"    _B{i}P0-_B{i}P1*_Q*(_it+{xs}[1])/_d1+_ev:\n"
                f"    _B{i}P0-_B{i}P1*_Q*{xs}[1]/fmax(_it-0.1*_Q,1e-9*_Q)-_B{i}P1*_Q*_it/_d1+_ev;\n"
                f"  _sig[{o0}]=_E-_B{i}P5*_i; _sig[{o1}]=1.0-_it/_Q; }}"
            )

        if op == ops.OP_BATERIA_ECM:
            iload = _idx(b.in_idx, 0)
            o = [int(b.out_idx[k]) for k in range(6)]
            n_ocv = int(b.param[11])
            po = 12 + 2*n_ocv
            return (
                f"{{ double _Q=_B{i}P0, _R1=_B{i}P5, _C1=_B{i}P6, _R2=_B{i}P7, _C2=_B{i}P8;\n"
                f"  double _Ns=_B{i}P9, _Np=_B{i}P10, _Rth=_B{i}P{po+6}, _Cth=_B{i}P{po+7}, _Tamb=_B{i}P{po+8};\n"
                f"  double _Icell = _sig[{iload}] / _Np;\n"
                f"  {xs}[0] += {dt} * (-_Icell / (_Q * 3600.0));\n"
                f"  {xs}[1] += {dt} * ((_Icell - {xs}[1]/_R1)/_C1);\n"
                f"  {xs}[2] += {dt} * ((_Icell - {xs}[2]/_R2)/_C2);\n"
                f"  double _R0=_B{i}P4, _Ploss_c = _Icell*_Icell*_R0;\n"
                f"  if(_R1>1e-6) _Ploss_c += {xs}[1]*{xs}[1]/_R1;\n"
                f"  if(_R2>1e-6) _Ploss_c += {xs}[2]*{xs}[2]/_R2;\n"
                f"  double _Ploss = _Ploss_c * _Ns * _Np;\n"
                f"  {xs}[3] += {dt} * ((_Ploss - ({xs}[3] - _Tamb)/_Rth)/_Cth);\n"
                f"  if({xs}[0]<0.0) {xs}[0]=0.0; if({xs}[0]>1.0) {xs}[0]=1.0;\n"
                f"  double _soc={xs}[0], _ocv=0.0;\n"
                f"  /* Interpolacion OCV harcodeada inline temporalmente */\n"
                f"  double *_lut = &_b{i}_param[12]; int _n_ocv = (int)_b{i}_param[11];\n"
                f"  if (_soc <= _lut[0]) _ocv = _lut[1];\n"
                f"  else if (_soc >= _lut[2*(_n_ocv-1)]) _ocv = _lut[2*(_n_ocv-1)+1];\n"
                f"  else {{ for(int _k=0; _k<_n_ocv-1; _k++){{ \n"
                f"      if(_soc>=_lut[2*_k] && _soc<=_lut[2*_k+2]) {{\n"
                f"          _ocv = _lut[2*_k+1] + (_lut[2*_k+3]-_lut[2*_k+1])*(_soc-_lut[2*_k])/(_lut[2*_k+2]-_lut[2*_k]); break;\n"
                f"      }} }} }}\n"
                f"  double _Vcell = _ocv - _Icell*_R0 - {xs}[1] - {xs}[2];\n"
                f"  _sig[{o[0]}] = _Vcell * _Ns;\n"
                f"  _sig[{o[1]}] = _soc;\n"
                f"  _sig[{o[2]}] = {xs}[3];\n"
                f"  _sig[{o[3]}] = _Ploss;\n"
                f"  double _Ichglim = _B{i}P{po+0} * _Np, _Idchlim = _B{i}P{po+1} * _Np;\n"
                f"  if({xs}[3]<_B{i}P{po+2} || {xs}[3]>_B{i}P{po+3} || _soc>=1.0) _Ichglim = 0.0;\n"
                f"  if({xs}[3]<_B{i}P{po+4} || {xs}[3]>_B{i}P{po+5} || _soc<=0.0) _Idchlim = 0.0;\n"
                f"  _sig[{o[4]}] = _Ichglim; _sig[{o[5]}] = _Idchlim; }}"
            )

        if op == ops.OP_MNA:
            return self._gen_mna(i, b)

        if op == ops.OP_VCVS:
            in_p = _idx(b.in_idx, 0)
            in_n = _idx(b.in_idx, 1)
            o0 = int(b.out_idx[0])
            gain = b.param[0] if b.param else 1.0
            return (f"{{ double _vp = _sig[{in_p}]; double _vn = _sig[{in_n}];\n"
                    f"  double _vout = {_fmt_f(gain)} * (_vp - _vn);\n"
                    f"  {self._relax(o0, '_vout')} }}")

        if op == ops.OP_VCCS:
            in_p = _idx(b.in_idx, 0)
            in_n = _idx(b.in_idx, 1)
            o0 = int(b.out_idx[0])
            gm = b.param[0] if b.param else 1.0
            return (f"{{ double _vp = _sig[{in_p}]; double _vn = _sig[{in_n}];\n"
                    f"  double _iout = {_fmt_f(gm)} * (_vp - _vn);\n"
                    f"  {self._relax(o0, '_iout')} }}")

        return f"/* AVISO: opcode {b.op} ({b.etiqueta}) — sin generacion de codigo */"

    def _gen_mna(self, i: int, b) -> str:
        xs = f"_b{i}_x"
        in_list = ", ".join(str(int(v)) for v in b.in_idx) if b.in_idx else "-1"
        out_list = ", ".join(str(int(v)) for v in b.out_idx) if b.out_idx else "-1"
        code = []
        code.append("{")
        code.append(f"    double *p = _b{i}_param;")
        code.append(f"    int _nparam = _b{i}_nparam;")
        code.append(f"    int nx = (int)p[0];")
        code.append(f"    int nu = (int)p[1];")
        code.append(f"    int n_sw_ctrl = (int)p[2];")
        code.append(f"    int n_out = (int)p[3];")
        code.append(f"    int n_diodos = (int)p[4];")
        code.append(f"    int n_mv = (int)p[5];")
        code.append(f"    if (nx <=0 || nx > MNA_MAX_NX) {{}} else {{")
        code.append(f"        int p_mv_base = 6;")
        code.append(f"        int n_mi = (int)p[p_mv_base + 2*n_mv];")
        code.append(f"        int idx_sw_base = p_mv_base + 2*n_mv + 1 + n_mi;")
        code.append(f"        int idx_diodo_base = idx_sw_base + n_sw_ctrl;")
        code.append(f"        int pos = idx_diodo_base + 5*n_diodos;")
        code.append(f"        int n_R = (int)p[pos]; pos++;")
        code.append(f"        int r_base = pos; pos += 3*n_R;")
        code.append(f"        int n_C = (int)p[pos]; pos++;")
        code.append(f"        int c_base = pos; pos += 3*n_C;")
        code.append(f"        int n_L = (int)p[pos]; pos++;")
        code.append(f"        int l_base = pos; pos += 4*n_L;")
        code.append(f"        int n_VS = (int)p[pos]; pos++;")
        code.append(f"        int vs_base = pos; pos += 4*n_VS;")
        code.append(f"        int n_IS = (int)p[pos]; pos++;")
        code.append(f"        int is_base = pos; pos += 3*n_IS;")
        code.append(f"        int n_SW_topo = (int)p[pos]; pos++;")
        code.append(f"        int sw_topo_base = pos; pos += 4*n_SW_topo;")
        code.append(f"        int n_vcvs = (int)p[pos]; pos++;")
        code.append(f"        int vcvs_base = pos; pos += 5*n_vcvs;")
        code.append(f"        int n_vccs = (int)p[pos]; pos++;")
        code.append(f"        int vccs_base = pos; pos += 4*n_vccs;")
        code.append(f"        int n_mut = (int)p[pos]; pos++;")
        code.append(f"        int mut_base = pos; pos += 9*n_mut;")
        code.append(f"        int metodo = 0;")
        code.append(f"        int is_precomputed = 0;")
        code.append(f"        int num_states = 0;")
        code.append(f"        double *pre_base = NULL;")
        code.append(f"        if (pos < _nparam) {{ metodo = (int)p[pos]; pos++; if (pos < _nparam) {{ is_precomputed = (int)p[pos]; pos++; if (is_precomputed==1 && pos < _nparam) {{ num_states = (int)p[pos]; pos++; pre_base = &p[pos]; }} }} }}")
        code.append(f"        double dt = _DT; if (dt < 1e-12) dt = 1e-12; double inv_dt = 1.0/dt; if (metodo==1) inv_dt = 2.0/dt;")
        code.append(f"        double x_prev[128]; for(int _k=0;_k<nx;_k++) x_prev[_k] = {xs}[_k];")
        code.append(f"        double u_local[128]={{0}};")
        if b.in_idx:
            for idx_j, sig_idx in enumerate(b.in_idx):
                if sig_idx >= 0:
                    code.append(f"        if ({idx_j}<nu) u_local[{idx_j}] = _sig[{int(sig_idx)}];")
        code.append(f"        if (is_precomputed && pre_base) {{")
        code.append(f"            int diode_state[64]={{0}};")
        code.append(f"            for(int d=0; d<n_diodos; d++){{ int n1=(int)p[idx_diodo_base+5*d]; int n2=(int)p[idx_diodo_base+5*d+1]; double Vf=p[idx_diodo_base+5*d+2]; double Va=(n1>=0&&n1<nx)?x_prev[n1]:0; double Vc=(n2>=0&&n2<nx)?x_prev[n2]:0; diode_state[d]=((Va-Vc)>Vf)?1:0; }}")
        code.append(f"            int sw_state[64]={{0}};")
        code.append(f"            for(int s=0;s<n_sw_ctrl;s++){{ int idx_u=(int)p[idx_sw_base+s]; double cur=0; if(idx_u>=0 && idx_u < nu) cur=u_local[idx_u]; sw_state[s]=(cur>0.5)?1:0; }}")
        code.append(f"            double x_new[128]={{0}};")
        code.append(f"            for(int _iter=0; _iter<5; _iter++){{")
        code.append(f"                int estado=0;")
        code.append(f"                for(int s=0;s<n_sw_ctrl;s++) if(sw_state[s]) estado|=(1<<s);")
        code.append(f"                for(int d=0;d<n_diodos;d++) if(diode_state[d]) estado|=(1<<(n_sw_ctrl+d));")
        code.append(f"                if(estado<0) estado=0; if(estado>=num_states) estado=num_states-1;")
        code.append(f"                int Bx_sz = nx*nx; int Bu_sz = nx*nu; int blk = Bx_sz+Bu_sz;")
        code.append(f"                double *Bx = pre_base + estado*blk;")
        code.append(f"                double *Bu = pre_base + estado*blk + Bx_sz;")
        code.append(f"                double x_tmp[128]={{0}};")
        code.append(f"                for(int _ii=0;_ii<nx;_ii++){{ double acc=0; for(int _jj=0;_jj<nx;_jj++) acc+= Bx[_ii*nx+_jj]*x_prev[_jj]; for(int _jj=0;_jj<nu;_jj++) acc+= Bu[_ii*nu+_jj]*u_local[_jj]; x_tmp[_ii]=acc; }}")
        code.append(f"                double b_diode_pre[128]={{0}};")
        code.append(f"                int need_vf=0;")
        code.append(f"                for(int d=0;d<n_diodos;d++) if(diode_state[d]){{ double Vf=p[idx_diodo_base+5*d+2]; if(fabs(Vf)>1e-12) need_vf=1; int n1=(int)p[idx_diodo_base+5*d]; int n2=(int)p[idx_diodo_base+5*d+1]; double Ron=p[idx_diodo_base+5*d+3]; double g=1.0/Ron; if(n1>=0) b_diode_pre[n1]+=g*Vf; if(n2>=0) b_diode_pre[n2]-=g*Vf; }}")
        code.append(f"                if(need_vf){{")
        code.append(f"#ifdef THREAD_LOCAL_FALLBACK_CALLOC")
        code.append(f"                    double *Gp_data = (double*)calloc(128*128, sizeof(double));")
        code.append(f"                    double *Cp_data = (double*)calloc(128*128, sizeof(double));")
        code.append(f"                    double *Ap_data = (double*)calloc(128*128, sizeof(double));")
        code.append(f"                    if(!Gp_data || !Cp_data || !Ap_data){{ free(Gp_data); free(Cp_data); free(Ap_data); return; }}")
        code.append(f"                    double (*Gp)[128] = (double(*)[128])Gp_data;")
        code.append(f"                    double (*Cp)[128] = (double(*)[128])Cp_data;")
        code.append(f"                    double (*Ap)[128] = (double(*)[128])Ap_data;")
        code.append(f"#else")
        code.append(f"                    double (*Gp)[128] = mna_G; double (*Cp)[128] = mna_Cmat; double (*Ap)[128] = mna_Ap;")
        code.append(f"                    for(int _i=0;_i<nx;_i++) for(int _j=0;_j<nx;_j++){{ Gp[_i][_j]=0; Cp[_i][_j]=0; Ap[_i][_j]=0; }}")
        code.append(f"#endif")
        code.append(f"                    for(int _k=0;_k<n_R;_k++){{ int n1=(int)p[r_base+3*_k]; int n2=(int)p[r_base+3*_k+1]; double R=p[r_base+3*_k+2]; if(R<1e-12)R=1e-12; double g=1.0/R; if(n1>=0) Gp[n1][n1]+=g; if(n2>=0) Gp[n2][n2]+=g; if(n1>=0&&n2>=0){{Gp[n1][n2]-=g; Gp[n2][n1]-=g;}}}}")
        code.append(f"                    for(int _k=0;_k<n_C;_k++){{ int n1=(int)p[c_base+3*_k]; int n2=(int)p[c_base+3*_k+1]; double Cval=p[c_base+3*_k+2]; if(n1>=0) Cp[n1][n1]+=Cval; if(n2>=0) Cp[n2][n2]+=Cval; if(n1>=0&&n2>=0){{Cp[n1][n2]-=Cval; Cp[n2][n1]-=Cval;}}}}")
        code.append(f"                    for(int _k=0;_k<n_L;_k++){{ int n1=(int)p[l_base+4*_k]; int n2=(int)p[l_base+4*_k+1]; int id=(int)p[l_base+4*_k+2]; double Lval=p[l_base+4*_k+3]; if(Lval<1e-12)Lval=1e-12; if(n1>=0){{Gp[n1][id]+=1; Gp[id][n1]+=1;}} if(n2>=0){{Gp[n2][id]-=1; Gp[id][n2]-=1;}} if(id>=0&&id<nx) Cp[id][id]=-Lval; }}")
        code.append(f"                    for(int _k=0;_k<n_VS;_k++){{ int n1=(int)p[vs_base+4*_k]; int n2=(int)p[vs_base+4*_k+1]; int ii=(int)p[vs_base+4*_k+3]; if(n1>=0){{Gp[n1][ii]+=1; Gp[ii][n1]+=1;}} if(n2>=0){{Gp[n2][ii]-=1; Gp[ii][n2]-=1;}} }}")
        code.append(f"                    for(int _k=0;_k<n_sw_ctrl;_k++){{ int n1=(int)p[sw_topo_base+4*_k]; int n2=(int)p[sw_topo_base+4*_k+1]; double Ron=p[sw_topo_base+4*_k+2]; double Roff=p[sw_topo_base+4*_k+3]; double R= ((estado>>_k)&1)?Ron:Roff; if(R<1e-12)R=1e-12; double g=1.0/R; if(n1>=0) Gp[n1][n1]+=g; if(n2>=0) Gp[n2][n2]+=g; if(n1>=0&&n2>=0){{Gp[n1][n2]-=g; Gp[n2][n1]-=g;}}}}")
        code.append(f"                    for(int _k=0;_k<n_diodos;_k++){{ int n1=(int)p[idx_diodo_base+5*_k]; int n2=(int)p[idx_diodo_base+5*_k+1]; double Ron=p[idx_diodo_base+5*_k+3]; double Roff=p[idx_diodo_base+5*_k+4]; int bit=(estado>>(n_sw_ctrl+_k))&1; double R=bit?Ron:Roff; if(R<1e-12)R=1e-12; double g=1.0/R; if(n1>=0) Gp[n1][n1]+=g; if(n2>=0) Gp[n2][n2]+=g; if(n1>=0&&n2>=0){{Gp[n1][n2]-=g; Gp[n2][n1]-=g;}}}}")
        code.append(f"                    for(int _k=0;_k<n_vcvs;_k++){{ int n1=(int)p[vcvs_base+5*_k]; int n2=(int)p[vcvs_base+5*_k+1]; int ii=(int)p[vcvs_base+5*_k+3]; if(n1>=0){{ Gp[n1][ii]+=1; Gp[ii][n1]+=1; }} if(n2>=0){{ Gp[n2][ii]-=1; Gp[ii][n2]-=1; }} }}")
        code.append(f"                    for(int _k=0;_k<n_mut;_k++){{ int n1=(int)p[mut_base+9*_k], n2=(int)p[mut_base+9*_k+1], n3=(int)p[mut_base+9*_k+2], n4=(int)p[mut_base+9*_k+3], i1=(int)p[mut_base+9*_k+4], i2=(int)p[mut_base+9*_k+5]; double L1=p[mut_base+9*_k+6], L2=p[mut_base+9*_k+7], M=p[mut_base+9*_k+8]; if(n1>=0){{ Gp[n1][i1]+=1; Gp[i1][n1]+=1; }} if(n2>=0){{ Gp[n2][i1]-=1; Gp[i1][n2]-=1; }} if(n3>=0){{ Gp[n3][i2]+=1; Gp[i2][n3]+=1; }} if(n4>=0){{ Gp[n4][i2]-=1; Gp[i2][n4]-=1; }} if(i1>=0 && i1<nx) Cp[i1][i1] = -L1; if(i2>=0 && i2<nx) Cp[i2][i2] = -L2; if(i1>=0 && i1<nx && i2>=0 && i2<nx){{ Cp[i1][i2] = -M; Cp[i2][i1] = -M; }} }}")
        code.append(f"                    for(int ii=0;ii<nx;ii++) for(int jj=0;jj<nx;jj++) Ap[ii][jj]=Gp[ii][jj]+Cp[ii][jj]*inv_dt;")
        code.append(f"                    for(int ii=0;ii<nx;ii++){{ double rs=0; for(int jj=0;jj<nx;jj++) rs+=fabs(Ap[ii][jj]); if(rs<1e-14) Ap[ii][ii]+=1e-12; }}")
        code.append(f"                    double delta[128]={{0}}; lu_solve_dense(Ap,b_diode_pre,delta,nx); for(int ii=0;ii<nx;ii++) x_tmp[ii]+=delta[ii];")
        code.append(f"#ifdef THREAD_LOCAL_FALLBACK_CALLOC")
        code.append(f"                    free(Gp_data); free(Cp_data); free(Ap_data);")
        code.append(f"#endif")
        code.append(f"                }}")
        code.append(f"                for(int ii=0;ii<nx;ii++) x_new[ii]=x_tmp[ii];")
        code.append(f"                int changed=0;")
        code.append(f"                for(int d=0;d<n_diodos;d++){{ int n1=(int)p[idx_diodo_base+5*d]; int n2=(int)p[idx_diodo_base+5*d+1]; double Vf=p[idx_diodo_base+5*d+2]; double Ron=p[idx_diodo_base+5*d+3]; double Va=(n1>=0&&n1<nx)?x_new[n1]:0; double Vc=(n2>=0&&n2<nx)?x_new[n2]:0; if(diode_state[d]==0){{ if((Va-Vc)>Vf){{diode_state[d]=1; changed=1;}} }} else {{ double Id=(Va-Vc-Vf)/Ron; if(Id < -1e-9){{diode_state[d]=0; changed=1;}} }} }}")
        code.append(f"                if(!changed) break;")
        code.append(f"            }}")
        code.append(f"            for(int ii=0;ii<nx;ii++) {xs}[ii]=x_new[ii];")
        if b.out_idx:
            code.append(f"            {{")
            code.append(f"                // salidas V(n1)-V(n2) e I(idx_i)")
            code.append(f"                int p_mv_base2 = 6;")
            for mv_idx in range(int(b.param[5]) if len(b.param)>5 else 0):
                sig_idx = int(b.out_idx[mv_idx]) if mv_idx < len(b.out_idx) else -1
                if sig_idx >=0:
                    code.append(f"                {{ int n1=(int)p[p_mv_base2+2*{mv_idx}]; int n2=(int)p[p_mv_base2+2*{mv_idx}+1]; double v1=(n1>=0&&n1<nx)?x_new[n1]:0; double v2=(n2>=0&&n2<nx)?x_new[n2]:0; _sig[{sig_idx}] = v1-v2; }}")
            n_mv_val = int(b.param[5]) if len(b.param)>5 else 0
            n_mi_param = int(b.param[6+2*n_mv_val]) if len(b.param) > 6+2*n_mv_val else 0
            for mi_idx in range(n_mi_param):
                sig_pos = n_mv_val + mi_idx
                if sig_pos < len(b.out_idx):
                    sig_idx = int(b.out_idx[sig_pos])
                    code.append(f"                {{ int p_mi_base = p_mv_base2+2*n_mv+1; int id=(int)p[p_mi_base+{mi_idx}]; double cur=(id>=0&&id<nx)?x_new[id]:0; _sig[{sig_idx}] = cur; }}")
            code.append(f"            }}")
        code.append(f"        }} else {{")
        code.append(f"            int sw_state[64]={{0}};")
        code.append(f"            for(int s=0;s<n_sw_ctrl;s++){{ int idx_u=(int)p[idx_sw_base+s]; double cur=0; if(idx_u>=0 && idx_u < nu) cur=u_local[idx_u]; sw_state[s]=(cur>0.5)?1:0; }}")
        code.append(f"            int diode_state[64]={{0}};")
        code.append(f"            for(int d=0; d<n_diodos; d++){{ int n1=(int)p[idx_diodo_base+5*d]; int n2=(int)p[idx_diodo_base+5*d+1]; double Vf=p[idx_diodo_base+5*d+2]; double Va=(n1>=0&&n1<nx)?x_prev[n1]:0; double Vc=(n2>=0&&n2<nx)?x_prev[n2]:0; diode_state[d]=((Va-Vc)>Vf)?1:0; }}")
        code.append(f"            double x_new[128]={{0}};")
        code.append(f"            for(int _diter=0; _diter<5; _diter++){{")
        code.append(f"#ifdef THREAD_LOCAL_FALLBACK_CALLOC")
        code.append(f"                double *G_data = (double*)calloc(128*128, sizeof(double));")
        code.append(f"                double *Cmat_data = (double*)calloc(128*128, sizeof(double));")
        code.append(f"                double *W_data = (double*)calloc(128*128, sizeof(double));")
        code.append(f"                double *b_diode_data = (double*)calloc(128, sizeof(double));")
        code.append(f"                if(!G_data || !Cmat_data || !W_data || !b_diode_data){{ free(G_data); free(Cmat_data); free(W_data); free(b_diode_data); return; }}")
        code.append(f"                double (*G)[128] = (double(*)[128])G_data;")
        code.append(f"                double (*Cmat)[128] = (double(*)[128])Cmat_data;")
        code.append(f"                double (*W)[128] = (double(*)[128])W_data;")
        code.append(f"                double *b_diode = b_diode_data;")
        code.append(f"#else")
        code.append(f"                double (*G)[128] = mna_G; double (*Cmat)[128] = mna_Cmat; double (*W)[128] = mna_W; double *b_diode = mna_b_diode;")
        code.append(f"#endif")
        code.append(f"                for(int _i=0;_i<nx;_i++){{ b_diode[_i]=0; for(int _j=0;_j<nx;_j++){{ G[_i][_j]=0; Cmat[_i][_j]=0; W[_i][_j]=0; }} }}")
        code.append(f"                for(int _k=0;_k<n_R;_k++){{ int n1=(int)p[r_base+3*_k]; int n2=(int)p[r_base+3*_k+1]; double R=p[r_base+3*_k+2]; if(R<1e-12)R=1e-12; double g=1.0/R; if(n1>=0) G[n1][n1]+=g; if(n2>=0) G[n2][n2]+=g; if(n1>=0&&n2>=0){{G[n1][n2]-=g; G[n2][n1]-=g;}}}}")
        code.append(f"                for(int _k=0;_k<n_C;_k++){{ int n1=(int)p[c_base+3*_k]; int n2=(int)p[c_base+3*_k+1]; double Cval=p[c_base+3*_k+2]; if(n1>=0) Cmat[n1][n1]+=Cval; if(n2>=0) Cmat[n2][n2]+=Cval; if(n1>=0&&n2>=0){{Cmat[n1][n2]-=Cval; Cmat[n2][n1]-=Cval;}}}}")
        code.append(f"                for(int _k=0;_k<n_L;_k++){{ int n1=(int)p[l_base+4*_k]; int n2=(int)p[l_base+4*_k+1]; int id=(int)p[l_base+4*_k+2]; double Lval=p[l_base+4*_k+3]; if(Lval<1e-12)Lval=1e-12; if(n1>=0){{G[n1][id]+=1; G[id][n1]+=1;}} if(n2>=0){{G[n2][id]-=1; G[id][n2]-=1;}} if(id>=0&&id<nx) Cmat[id][id]=-Lval; }}")
        code.append(f"                for(int _k=0;_k<n_VS;_k++){{ int n1=(int)p[vs_base+4*_k]; int n2=(int)p[vs_base+4*_k+1]; int iu=(int)p[vs_base+4*_k+2]; int ii=(int)p[vs_base+4*_k+3]; if(n1>=0){{G[n1][ii]+=1; G[ii][n1]+=1;}} if(n2>=0){{G[n2][ii]-=1; G[ii][n2]-=1;}} if(ii>=0&&ii<nx&&iu>=0&&iu<nu) W[ii][iu]=1.0; }}")
        code.append(f"                for(int _k=0;_k<n_IS;_k++){{ int n1=(int)p[is_base+3*_k]; int n2=(int)p[is_base+3*_k+1]; int iu=(int)p[is_base+3*_k+2]; if(n1>=0&&iu>=0&&iu<nu) W[n1][iu]+=-1.0; if(n2>=0&&iu>=0&&iu<nu) W[n2][iu]+=1.0; }}")
        code.append(f"                for(int _k=0;_k<n_sw_ctrl;_k++){{ int n1=(int)p[sw_topo_base+4*_k]; int n2=(int)p[sw_topo_base+4*_k+1]; double Ron=p[sw_topo_base+4*_k+2]; double Roff=p[sw_topo_base+4*_k+3]; double R= sw_state[_k]?Ron:Roff; if(R<1e-12)R=1e-12; double g=1.0/R; if(n1>=0) G[n1][n1]+=g; if(n2>=0) G[n2][n2]+=g; if(n1>=0&&n2>=0){{G[n1][n2]-=g; G[n2][n1]-=g;}}}}")
        code.append(f"                for(int _k=0;_k<n_diodos;_k++){{ int n1=(int)p[idx_diodo_base+5*_k]; int n2=(int)p[idx_diodo_base+5*_k+1]; double Vf=p[idx_diodo_base+5*_k+2]; double Ron=p[idx_diodo_base+5*_k+3]; double Roff=p[idx_diodo_base+5*_k+4]; double R= diode_state[_k]?Ron:Roff; if(R<1e-12)R=1e-12; double g=1.0/R; if(n1>=0) G[n1][n1]+=g; if(n2>=0) G[n2][n2]+=g; if(n1>=0&&n2>=0){{G[n1][n2]-=g; G[n2][n1]-=g;}} if(diode_state[_k]){{ if(n1>=0) b_diode[n1]+=g*Vf; if(n2>=0) b_diode[n2]-=g*Vf; }} }}")
        code.append(f"                for(int _k=0;_k<n_vcvs;_k++){{ int n1=(int)p[vcvs_base+5*_k]; int n2=(int)p[vcvs_base+5*_k+1]; int iu=(int)p[vcvs_base+5*_k+2]; int ii=(int)p[vcvs_base+5*_k+3]; double gain=p[vcvs_base+5*_k+4]; if(n1>=0){{G[n1][ii]+=1; G[ii][n1]+=1;}} if(n2>=0){{G[n2][ii]-=1; G[ii][n2]-=1;}} if(ii>=0&&ii<nx&&iu>=0&&iu<nu) W[ii][iu]=gain; }}")
        code.append(f"                for(int _k=0;_k<n_vccs;_k++){{ int n1=(int)p[vccs_base+4*_k]; int n2=(int)p[vccs_base+4*_k+1]; int iu=(int)p[vccs_base+4*_k+2]; double gm=p[vccs_base+4*_k+3]; if(n1>=0&&iu>=0&&iu<nu) W[n1][iu]-=gm; if(n2>=0&&iu>=0&&iu<nu) W[n2][iu]+=gm; }}")
        code.append(f"                for(int _k=0;_k<n_mut;_k++){{ int n1=(int)p[mut_base+9*_k]; int n2=(int)p[mut_base+9*_k+1]; int n3=(int)p[mut_base+9*_k+2]; int n4=(int)p[mut_base+9*_k+3]; int i1=(int)p[mut_base+9*_k+4]; int i2=(int)p[mut_base+9*_k+5]; double L1=p[mut_base+9*_k+6]; double L2=p[mut_base+9*_k+7]; double M=p[mut_base+9*_k+8]; if(n1>=0){{G[n1][i1]+=1; G[i1][n1]+=1;}} if(n2>=0){{G[n2][i1]-=1; G[i1][n2]-=1;}} if(n3>=0){{G[n3][i2]+=1; G[i2][n3]+=1;}} if(n4>=0){{G[n4][i2]-=1; G[i2][n4]-=1;}} if(i1>=0&&i1<nx) Cmat[i1][i1]=-L1; if(i2>=0&&i2<nx) Cmat[i2][i2]=-L2; if(i1>=0&&i1<nx&&i2>=0&&i2<nx){{Cmat[i1][i2]=-M; Cmat[i2][i1]=-M;}} }}")
        code.append(f"#ifdef THREAD_LOCAL_FALLBACK_CALLOC")
        code.append(f"                double *A_data = (double*)calloc(128*128, sizeof(double));")
        code.append(f"                double *b_data = (double*)calloc(128, sizeof(double));")
        code.append(f"                if(!A_data || !b_data){{ free(A_data); free(b_data); free(G_data); free(Cmat_data); free(W_data); free(b_diode_data); return; }}")
        code.append(f"                double (*A)[128] = (double(*)[128])A_data; double *bvec = b_data;")
        code.append(f"#else")
        code.append(f"                double (*A)[128] = mna_Ap; double *bvec = mna_b;")
        code.append(f"                for(int _i=0;_i<nx;_i++){{ bvec[_i]=0; for(int _j=0;_j<nx;_j++){{ A[_i][_j]=0; }} }}")
        code.append(f"#endif")
        code.append(f"                for(int ii=0;ii<nx;ii++) for(int jj=0;jj<nx;jj++) A[ii][jj]=G[ii][jj]+Cmat[ii][jj]*inv_dt;")
        code.append(f"                for(int ii=0;ii<nx;ii++){{ double rs=0; for(int jj=0;jj<nx;jj++) rs+=fabs(A[ii][jj]); if(rs<1e-14) A[ii][ii]+=1e-12; }}")
        code.append(f"                for(int ii=0;ii<nx;ii++){{ double acc=b_diode[ii]; for(int jj=0;jj<nu;jj++) acc+=W[ii][jj]*u_local[jj]; for(int jj=0;jj<nx;jj++) acc+=Cmat[ii][jj]*inv_dt*x_prev[jj]; bvec[ii]=acc; }}")
        code.append(f"                lu_solve_dense(A,bvec,x_new,nx);")
        code.append(f"#ifdef THREAD_LOCAL_FALLBACK_CALLOC")
        code.append(f"                free(A_data); free(b_data);")
        code.append(f"#endif")
        code.append(f"                int ch=0;")
        code.append(f"                for(int d=0;d<n_diodos;d++){{ int n1=(int)p[idx_diodo_base+5*d]; int n2=(int)p[idx_diodo_base+5*d+1]; double Vf=p[idx_diodo_base+5*d+2]; double Ron=p[idx_diodo_base+5*d+3]; double Va=(n1>=0&&n1<nx)?x_new[n1]:0; double Vc=(n2>=0&&n2<nx)?x_new[n2]:0; if(diode_state[d]==0){{ if((Va-Vc)>Vf){{diode_state[d]=1; ch=1;}} }} else {{ double Id=(Va-Vc-Vf)/Ron; if(Id < -1e-9){{diode_state[d]=0; ch=1;}} }} }}")
        code.append(f"                if(!ch) break;")
        code.append(f"            }}")
        code.append(f"#ifdef THREAD_LOCAL_FALLBACK_CALLOC")
        code.append(f"            free(G_data); free(Cmat_data); free(W_data); free(b_diode_data);")
        code.append(f"#endif")
        code.append(f"            for(int ii=0;ii<nx;ii++) {xs}[ii]=x_new[ii];")
        if b.out_idx:
            code.append(f"            {{")
            for mv_idx in range(int(b.param[5]) if len(b.param)>5 else 0):
                sig_idx = int(b.out_idx[mv_idx]) if mv_idx < len(b.out_idx) else -1
                if sig_idx >=0:
                    code.append(f"                {{ int n1=(int)p[p_mv_base+2*{mv_idx}]; int n2=(int)p[p_mv_base+2*{mv_idx}+1]; double v1=(n1>=0&&n1<nx)?x_new[n1]:0; double v2=(n2>=0&&n2<nx)?x_new[n2]:0; _sig[{sig_idx}] = v1-v2; }}")
            n_mv_val = int(b.param[5]) if len(b.param)>5 else 0
            n_mi_param = int(b.param[6+2*n_mv_val]) if len(b.param) > 6+2*n_mv_val else 0
            for mi_idx in range(n_mi_param):
                sig_pos = n_mv_val + mi_idx
                if sig_pos < len(b.out_idx):
                    sig_idx = int(b.out_idx[sig_pos])
                    code.append(f"                {{ int p_mi_base = p_mv_base+2*n_mv+1; int id=(int)p[p_mi_base+{mi_idx}]; double cur=(id>=0&&id<nx)?x_new[id]:0; _sig[{sig_idx}] = cur; }}")
            code.append(f"            }}")
        code.append(f"        }}")
        code.append(f"    }}")
        code.append("}")
        return "\n".join(code)

    def _gen_api(self, nombre: str) -> str:
        lines = [
            "/* === API Publica === */\n",
            f"void {nombre}_init(void) {{",
            "    memset(_sig, 0, sizeof(_sig));",
            "    _t = 0.0;",
        ]
        for i, b in enumerate(self._bloques):
            if b.n_state > 0:
                for j, v in enumerate(b.estados_iniciales):
                    lines.append(f"    _b{i}_x[{j}] = {_fmt_f(v)};")
                lines.append(f"    memset(_b{i}_ws, 0, sizeof(_b{i}_ws));")
        lines += [
            "    _actualizar_fuentes();",
            "    _lazo_algebraico();",
            "}",
            "",
            f"void {nombre}_paso(void) {{",
            "    _t += _DT;",
            "    _actualizar_fuentes();",
            "    _lazo_algebraico();",
            "    _actualizar_dinamicos();",
            "    _lazo_algebraico();",
            "}",
            "",
            f"void {nombre}_run(int n_steps, double *buf, int n_out, int *idx) {{",
            f"    {nombre}_init();",
            "    for (int s = 0; s < n_steps; s++) {",
            "        for (int r = 0; r < n_out; r++)",
            "            buf[r * n_steps + s] = _sig[idx[r]];",
            f"        {nombre}_paso();",
            "    }",
            "}",
            "",
            f"double {nombre}_t(void)    {{ return _t; }}",
            f"double *{nombre}_sig(void) {{ return _sig; }}",
            "",
            "/* main() de ejemplo — COMENTAR para embeber en firmware */",
            "#ifdef CODEGEN_STANDALONE_MAIN",
            "int main(void) {",
            f"    {nombre}_init();",
            "    int n = 1000;",
            f"    printf(\"Simulando %d pasos (dt=%.2e)...\\n\", n, _DT);",
            "    for (int s = 0; s < n; s++) {",
            f"        {nombre}_paso();",
            "        if (s % 100 == 0)",
            f"            printf(\"t=%.4f  sig[0]=%.6f\\n\", {nombre}_t(), _sig[0]);",
            "    }",
            "    return 0;",
            "}",
            "#endif",
        ]
        return "\n".join(lines) + "\n"
