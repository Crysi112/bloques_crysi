import numpy as np
from . import opcodes as ops
SQ3 = np.sqrt(3.0)
def _clarke(va, vb, vc):
    al = (2.0 / 3.0) * (va - 0.5 * vb - 0.5 * vc)
    be = (1.0 / 3.0) * (vb - vc) * SQ3
    return al, be
def _inv_clarke(al, be):
    return al, -0.5 * al + SQ3 / 2.0 * be, -0.5 * al - SQ3 / 2.0 * be
def _park(al, be, th):
    c, s = np.cos(th), np.sin(th)
    return al * c + be * s, -al * s + be * c
def _inv_park(d, q, th):
    c, s = np.cos(th), np.sin(th)
    return d * c - q * s, d * s + q * c
def _busqueda_bp(bp, u):
    if u <= bp[0]:
        return 0
    if u >= bp[-1]:
        return len(bp) - 2
    lo, hi = 0, len(bp) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if bp[mid] <= u:
            lo = mid
        else:
            hi = mid
    return lo
class _Ctx:
    def __init__(self, bloques, dt, metodo, max_iter, tol, w_opt, n_sig):
        self.bloques = bloques
        self.dt = dt
        self.metodo = metodo
        self.max_iter = max_iter
        self.tol = tol
        self.w_opt = w_opt
        self.sig = np.zeros(n_sig)
        self.t = 0.0
        self.est = [np.array(b.estados_iniciales, dtype=float) for b in bloques]
def _eval_estatico(b, ctx, maxdelta):
    s = ctx.sig
    w = ctx.w_opt
    op = b.op
    if op == ops.OP_GAIN:
        y = b.param[0] * s[b.in_idx[0]]
        o = s[b.out_idx[0]]
        nv = o + w * (y - o)
        s[b.out_idx[0]] = nv
        maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif op == ops.OP_SUM:
        y = sum(b.param[i] * s[idx] for i, idx in enumerate(b.in_idx))
        o = s[b.out_idx[0]]
        nv = o + w * (y - o)
        s[b.out_idx[0]] = nv
        maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif op == ops.OP_MUX:
        for k, j in enumerate(b.out_idx):
            y = s[b.in_idx[k]]
            o = s[j]
            nv = o + w * (y - o)
            s[j] = nv
            maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif op == ops.OP_DEMUX:
        for k, j in enumerate(b.out_idx):
            y = s[b.in_idx[k]]
            o = s[j]
            nv = o + w * (y - o)
            s[j] = nv
            maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif op == ops.OP_LUT1D:
        n = int(b.param[0])
        bx = b.param[1:1 + n]
        dy = b.param[1 + n:]
        u = s[b.in_idx[0]]
        i = _busqueda_bp(bx, u)
        f = min(max((u - bx[i]) / (bx[i + 1] - bx[i]), 0.0), 1.0)
        y = dy[i] + f * (dy[i + 1] - dy[i])
        o = s[b.out_idx[0]]
        nv = o + w * (y - o)
        s[b.out_idx[0]] = nv
        maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif op == ops.OP_LUT2D:
        nx, ny = int(b.param[0]), int(b.param[1])
        bx = b.param[2:2 + nx]
        by = b.param[2 + nx:2 + nx + ny]
        z = b.param[2 + nx + ny:]
        u1, u2 = s[b.in_idx[0]], s[b.in_idx[1]]
        i = _busqueda_bp(bx, u1)
        j = _busqueda_bp(by, u2)
        fx = min(max((u1 - bx[i]) / (bx[i + 1] - bx[i]), 0.0), 1.0)
        fy = min(max((u2 - by[j]) / (by[j + 1] - by[j]), 0.0), 1.0)
        z00 = z[j * nx + i]
        z01 = z[j * nx + i + 1]
        z10 = z[(j + 1) * nx + i]
        z11 = z[(j + 1) * nx + i + 1]
        y = ((z00 * (1 - fx) + z01 * fx) * (1 - fy)
             + (z10 * (1 - fx) + z11 * fx) * fy)
        o = s[b.out_idx[0]]
        nv = o + w * (y - o)
        s[b.out_idx[0]] = nv
        maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif op == ops.OP_LUT3D:
        nx, ny, nz = int(b.param[0]), int(b.param[1]), int(b.param[2])
        bx = b.param[3:3 + nx]
        by = b.param[3 + nx:3 + nx + ny]
        bz = b.param[3 + nx + ny:3 + nx + ny + nz]
        z = b.param[3 + nx + ny + nz:]
        u1, u2, u3 = s[b.in_idx[0]], s[b.in_idx[1]], s[b.in_idx[2]]
        i = _busqueda_bp(bx, u1)
        j = _busqueda_bp(by, u2)
        k = _busqueda_bp(bz, u3)
        fx = min(max((u1 - bx[i]) / (bx[i + 1] - bx[i]), 0.0), 1.0)
        fy = min(max((u2 - by[j]) / (by[j + 1] - by[j]), 0.0), 1.0)
        fz = min(max((u3 - bz[k]) / (bz[k + 1] - bz[k]), 0.0), 1.0)
        y = 0.0
        for kk in (0, 1):
            wz = 1 - fz if kk == 0 else fz
            for jj in (0, 1):
                wy = 1 - fy if jj == 0 else fy
                for ii in (0, 1):
                    wx = 1 - fx if ii == 0 else fx
                    y += (wx * wy * wz
                          * z[(k + kk) * ny * nx + (j + jj) * nx + (i + ii)])
        o = s[b.out_idx[0]]
        nv = o + w * (y - o)
        s[b.out_idx[0]] = nv
        maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif op == ops.OP_LOGICO:
        opc, umb = int(b.param[0]), b.param[1]
        if opc == 6:
            r = not (s[b.in_idx[0]] > umb)
        else:
            r = s[b.in_idx[0]] > umb
            for k in range(1, b.n_in):
                v = s[b.in_idx[k]] > umb
                if opc in (0, 2):
                    r = r and v
                elif opc in (1, 3):
                    r = r or v
                elif opc == 4:
                    r = r != v
                else:
                    r = not (r != v)
            if opc in (2, 3):
                r = not r
        y = 1.0 if r else 0.0
        o = s[b.out_idx[0]]
        nv = o + w * (y - o)
        s[b.out_idx[0]] = nv
        maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif op == ops.OP_RELACIONAL:
        opc, tol = int(b.param[0]), b.param[1]
        a, bb = s[b.in_idx[0]], s[b.in_idx[1]]
        if opc == 0:
            r = abs(a - bb) <= tol
        elif opc == 1:
            r = abs(a - bb) > tol
        elif opc == 2:
            r = a < bb
        elif opc == 3:
            r = a <= bb
        elif opc == 4:
            r = a > bb
        else:
            r = a >= bb
        y = 1.0 if r else 0.0
        o = s[b.out_idx[0]]
        nv = o + w * (y - o)
        s[b.out_idx[0]] = nv
        maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif op == ops.OP_CLARKE:
        al, be = _clarke(s[b.in_idx[0]], s[b.in_idx[1]], s[b.in_idx[2]])
        for j, val in zip(b.out_idx, (al, be)):
            o = s[j]
            nv = o + w * (val - o)
            s[j] = nv
            maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif op == ops.OP_INV_CLARKE:
        va, vb, vc = _inv_clarke(s[b.in_idx[0]], s[b.in_idx[1]])
        for j, val in zip(b.out_idx, (va, vb, vc)):
            o = s[j]
            nv = o + w * (val - o)
            s[j] = nv
            maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif op == ops.OP_PARK:
        d, q = _park(s[b.in_idx[0]], s[b.in_idx[1]], s[b.in_idx[2]])
        for j, val in zip(b.out_idx, (d, q)):
            o = s[j]
            nv = o + w * (val - o)
            s[j] = nv
            maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif op == ops.OP_INV_PARK:
        al, be = _inv_park(s[b.in_idx[0]], s[b.in_idx[1]], s[b.in_idx[2]])
        for j, val in zip(b.out_idx, (al, be)):
            o = s[j]
            nv = o + w * (val - o)
            s[j] = nv
            maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif op == ops.OP_QD:
        va, vb, vc, ia, ib, ic, th = (s[b.in_idx[k]] for k in range(7))
        alv, bev = _clarke(va, vb, vc)
        ali, bei = _clarke(ia, ib, ic)
        d, q = _park(alv, bev, th)
        vqs, vds = q, d
        d, q = _park(ali, bei, th)
        iqs, ids = q, d
        for j, val in zip(b.out_idx, (vqs, vds, iqs, ids)):
            o = s[j]
            nv = o + w * (val - o)
            s[j] = nv
            maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif op == ops.OP_SATURAR:
        u = s[b.in_idx[0]]
        lo, hi = b.param[0], b.param[1]
        y = lo if u < lo else (hi if u > hi else u)
        o = s[b.out_idx[0]]
        nv = o + w * (y - o)
        s[b.out_idx[0]] = nv
        maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif op == ops.OP_TRANSFORMADOR:
        a = b.param[0]
        fases = 3 if b.n_in > 3 else 1
        for k in range(fases):
            v1 = s[b.in_idx[k]]
            i2 = s[b.in_idx[fases + k]]
            for j, val in ((b.out_idx[k], a * v1),
                           (b.out_idx[fases + k], -a * i2)):
                o = s[j]
                nv = o + w * (val - o)
                s[j] = nv
                maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif op == ops.OP_MEDIDOR_POTENCIA:
        fases = (b.n_in - 2) // 2
        Pe = sum(s[b.in_idx[k]] * s[b.in_idx[fases + k]] for k in range(fases))
        if b.n_out == 3:
            v0, v1, v2 = s[b.in_idx[0]], s[b.in_idx[1]], s[b.in_idx[2]]
            i0, i1, i2 = s[b.in_idx[3]], s[b.in_idx[4]], s[b.in_idx[5]]
            vn = (v0 + v1 + v2) / 3.0
            v0p, v1p, v2p = v0 - vn, v1 - vn, v2 - vn
            inn = (i0 + i1 + i2) / 3.0
            i0p, i1p, i2p = i0 - inn, i1 - inn, i2 - inn
            SQRT23 = 0.8164965809277260
            INV_SQRT2 = 0.7071067811865475
            va = SQRT23 * (v0p - 0.5 * v1p - 0.5 * v2p)
            vb = INV_SQRT2 * (v1p - v2p)
            ia = SQRT23 * (i0p - 0.5 * i1p - 0.5 * i2p)
            ib = INV_SQRT2 * (i1p - i2p)
            Qe = vb * ia - va * ib
        else:
            Qe = 0.0
        te = s[b.in_idx[2 * fases]] if b.in_idx[2 * fases] >= 0 else 0.0
        wm = s[b.in_idx[2 * fases + 1]] if b.in_idx[2 * fases + 1] >= 0 else 0.0
        Pm = te * wm
        vals = [Pe, Qe] if b.n_out == 3 else [Pe]
        vals.append(Pm)
        for j, val in zip(b.out_idx, vals):
            o = s[j]
            nv = o + w * (val - o)
            s[j] = nv
            maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif op == ops.OP_INTERRUPTOR:
        g = s[b.in_idx[0]] if b.n_in == 3 else 0.0
        va, vc = s[b.in_idx[1]], s[b.in_idx[2]]
        r_on, r_off, modo = b.param[0], b.param[1], int(b.param[2])
        r = r_off
        if (modo == 0 and g > 0.5) or (modo == 1 and va - vc > 0.0):
            r = r_on
        v = va - vc
        for j, val in zip(b.out_idx, (v / r, v)):
            o = s[j]
            nv = o + w * (val - o)
            s[j] = nv
            maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif op == ops.OP_PUENTE_INV_3F:
        vdc = s[b.in_idx[0]]
        a, b2, c2 = s[b.in_idx[1]], s[b.in_idx[2]], s[b.in_idx[3]]
        if b.param[0] > 0.5:
            vals = (vdc * (a - 0.5 * b2 - 0.5 * c2) / 3.0,
                    vdc * (b2 - 0.5 * a - 0.5 * c2) / 3.0,
                    vdc * (c2 - 0.5 * a - 0.5 * b2) / 3.0)
        else:
            vals = (vdc * (2.0 * a - b2 - c2) / 3.0,
                    vdc * (2.0 * b2 - a - c2) / 3.0,
                    vdc * (2.0 * c2 - a - b2) / 3.0)
        for j, val in zip(b.out_idx, vals):
            o = s[j]
            nv = o + w * (val - o)
            s[j] = nv
            maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif op == ops.OP_PUENTE_INV_1F:
        vdc = s[b.in_idx[0]]
        a, b2 = s[b.in_idx[1]], s[b.in_idx[2]]
        y = vdc * (a - b2) / 2.0 if b.param[0] > 0.5 else vdc * (a - b2)
        o = s[b.out_idx[0]]
        nv = o + w * (y - o)
        s[b.out_idx[0]] = nv
        maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif op == ops.OP_PANEL_SOLAR:
        G, T, V = s[b.in_idx[0]], s[b.in_idx[1]], s[b.in_idx[2]]
        Ns, Np, Voc, Isc, ki = b.param[0], b.param[1], b.param[2], \
            b.param[3], b.param[6]
        Rs, Rsh, n = b.param[7], b.param[8], b.param[9]
        k_voc = b.param[10] if len(b.param) > 10 else 0.0
        Vt = 0.02585 * (T + 273.15) / 298.15
        a = Ns * n * Vt
        Iph = Np * (Isc + ki * (T - 25.0)) * (G / 1000.0)
        Voc_T = Voc * (1.0 + k_voc * (T - 25.0))
        I0 = (Iph - Voc_T / Rsh) / (np.exp(Voc_T / a) - 1.0)
        I = Iph - V / Rsh
        for _ in range(10):
            u = (V + I * Rs) / a
            if u > 700.0:
                u = 700.0
            e = np.exp(u)
            g = Iph - I0 * (e - 1.0) - (V + I * Rs) / Rsh - I
            gp = -I0 * (Rs / a) * e - Rs / Rsh - 1.0
            dI = g / gp
            lam, g1 = 1.0, abs(g)
            for _ in range(4):
                In = I - lam * dI
                u2 = (V + In * Rs) / a
                if u2 > 700.0:
                    u2 = 700.0
                g2 = abs(Iph - I0 * (np.exp(u2) - 1.0)
                         - (V + In * Rs) / Rsh - In)
                if g2 <= g1 or lam <= 0.125:
                    break
                lam *= 0.5
            I -= lam * dI
        o = s[b.out_idx[0]]
        nv = o + w * (I - o)
        s[b.out_idx[0]] = nv
        maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif op == ops.OP_PWM_1F:
        fsw, td_frac_param = b.param
        td_frac = td_frac_param * fsw
        d = s[b.in_idx[0]]
        d_eff = max(0.0, min(1.0, d - td_frac))
        phi = (ctx.t * fsw) % 1.0
        car = 2.0 * phi if phi < 0.5 else 2.0 * (1.0 - phi)
        val = 1.0 if car < d_eff else 0.0
        o = s[b.out_idx[0]]
        s[b.out_idx[0]] = val
        maxdelta[0] = max(maxdelta[0], abs(val - o))
    elif op == ops.OP_PWM_SPWM:
        f_out, fsw, fase, td_frac_param = b.param
        td_frac = td_frac_param * fsw
        mod = max(0.0, min(1.0, s[b.in_idx[0]]))
        omega = 2.0 * np.pi * f_out * ctx.t
        phi = (ctx.t * fsw) % 1.0
        car = 2.0 * phi if phi < 0.5 else 2.0 * (1.0 - phi)
        TWO_PI_3 = 2.0 * np.pi / 3.0
        for k in range(3):
            ref = 0.5 + 0.5 * mod * np.sin(omega - k * TWO_PI_3 + fase)
            d_eff = max(0.0, min(1.0, ref - td_frac))
            val = 1.0 if car < d_eff else 0.0
            o = s[b.out_idx[k]]
            s[b.out_idx[k]] = val
            maxdelta[0] = max(maxdelta[0], abs(val - o))
    elif op == ops.OP_PWM_SVPWM:
        Vdc, fsw, td_frac_param = b.param
        td_frac = td_frac_param * fsw
        Va, Vb = s[b.in_idx[0]], s[b.in_idx[1]]
        Vmax = Vdc / SQ3
        Vmag = np.sqrt(Va*Va + Vb*Vb)
        scale = min(Vmag, Vmax) / Vmag if Vmag > 1e-12 else 0.0
        Vn_a = Va * scale / Vmax
        Vn_b = Vb * scale / Vmax
        theta = np.arctan2(Vn_b, Vn_a)
        if theta < 0.0:
            theta += 2.0 * np.pi
        sector = int(theta / (np.pi / 3.0)) + 1
        sector = min(sector, 6)
        theta_s = theta - (sector - 1) * (np.pi / 3.0)
        mag = np.sqrt(Vn_a*Vn_a + Vn_b*Vn_b)
        T1 = mag * np.sin(np.pi / 3.0 - theta_s)
        T2 = mag * np.sin(theta_s)
        T0 = max(0.0, 1.0 - T1 - T2)
        h = T0 * 0.5
        if sector == 1:
            da, db, dc = T1 + T2 + h, T2 + h, h
        elif sector == 2:
            da, db, dc = T1 + h, T1 + T2 + h, h
        elif sector == 3:
            da, db, dc = h, T1 + T2 + h, T2 + h
        elif sector == 4:
            da, db, dc = h, T1 + h, T1 + T2 + h
        elif sector == 5:
            da, db, dc = T2 + h, h, T1 + T2 + h
        else:
            da, db, dc = T1 + T2 + h, h, T1 + h
        da = max(0.0, da - td_frac)
        db = max(0.0, db - td_frac)
        dc = max(0.0, dc - td_frac)
        phi = (ctx.t * fsw) % 1.0
        car = 2.0 * phi if phi < 0.5 else 2.0 * (1.0 - phi)
        vals = [1.0 if car < da else 0.0,
                1.0 if car < db else 0.0,
                1.0 if car < dc else 0.0]
        for k, val in enumerate(vals):
            o = s[b.out_idx[k]]
            s[b.out_idx[k]] = val
            maxdelta[0] = max(maxdelta[0], abs(val - o))
    elif b.op == ops.OP_RES_TERMICA:
        y = (s[b.in_idx[0]] - s[b.in_idx[1]]) / b.param[0]
        o = s[b.out_idx[0]]
        nv = o + w * (y - o)
        s[b.out_idx[0]] = nv
        maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif b.op == ops.OP_ENGRANAJE:
        a = b.param[0]
        for k, j in enumerate(b.out_idx):
            y = a * s[b.in_idx[0]] if k == 0 else s[b.in_idx[1]] / a
            o = s[j]
            nv = o + w * (y - o)
            s[j] = nv
            maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif b.op == ops.OP_EMBRAGUE:
        Tmax, umb = b.param[0], b.param[1]
        T = s[b.in_idx[0]]
        if s[b.in_idx[1]] > umb:
            y = max(min(T, Tmax), -Tmax)
        else:
            y = 0.0
        o = s[b.out_idx[0]]
        nv = o + w * (y - o)
        s[b.out_idx[0]] = nv
        maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif b.op == ops.OP_FALLO_PROG:
        t_fallo, valor, modo = b.param[0], b.param[1], b.param[2]
        u = s[b.in_idx[0]]
        y = (u + valor if modo > 0.5 else valor) if ctx.t >= t_fallo else u
        o = s[b.out_idx[0]]
        nv = o + w * (y - o)
        s[b.out_idx[0]] = nv
        maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif b.op == ops.OP_FALLO_EVENTO:
        umb, valor, modo = b.param[0], b.param[1], b.param[2]
        u = s[b.in_idx[0]]
        y = (u + valor if modo > 0.5 else valor) if s[b.in_idx[1]] > umb else u
        o = s[b.out_idx[0]]
        nv = o + w * (y - o)
        s[b.out_idx[0]] = nv
        maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif b.op == ops.OP_MULTIPLICADOR:
        y = s[b.in_idx[0]] * s[b.in_idx[1]]
        o = s[b.out_idx[0]]
        nv = o + w * (y - o)
        s[b.out_idx[0]] = nv
        maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif b.op == ops.OP_RESISTENCIA:
        y = s[b.in_idx[0]] / b.param[0]
        o = s[b.out_idx[0]]
        nv = o + w * (y - o)
        s[b.out_idx[0]] = nv
        maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif b.op == ops.OP_SAT_VECTORIAL:
        Vmax = b.param[0]
        vd = s[b.in_idx[0]]
        vq = s[b.in_idx[1]]
        mag = np.sqrt(vd * vd + vq * vq)
        scale = 1.0
        if mag > Vmax and mag > 1e-12:
            scale = Vmax / mag
        for k, val in enumerate([vd * scale, vq * scale]):
            o = s[b.out_idx[k]]
            nv = o + w * (val - o)
            s[b.out_idx[k]] = nv
            maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif op == ops.OP_CARGA_PQ_3F:
        va, vb, vc = s[b.in_idx[0]], s[b.in_idx[1]], s[b.in_idx[2]]
        vn = (va + vb + vc) / 3.0
        v0, v1, v2 = va - vn, vb - vn, vc - vn
        SQRT23 = 0.8164965809277260
        INV_SQRT2 = 0.7071067811865475
        INV_SQRT6 = 0.4082482904638630
        valpha = SQRT23 * (v0 - 0.5 * v1 - 0.5 * v2)
        vbeta = INV_SQRT2 * (v1 - v2)
        vmag2 = valpha * valpha + vbeta * vbeta
        if vmag2 < 1e-6:
            ia, ib, ic = 0.0, 0.0, 0.0
        else:
            P_total = b.param[0]
            Q_total = b.param[1]
            ialpha = (P_total * valpha + Q_total * vbeta) / vmag2
            ibeta  = (P_total * vbeta - Q_total * valpha) / vmag2
            ia = SQRT23 * ialpha
            ib = -INV_SQRT6 * ialpha + INV_SQRT2 * ibeta
            ic = -INV_SQRT6 * ialpha - INV_SQRT2 * ibeta
        for k, val in enumerate((ia, ib, ic)):
            o = s[b.out_idx[k]]
            nv = o + w * (val - o)
            s[b.out_idx[k]] = nv
            maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif b.op == ops.OP_CALCULO_IDC:
        vd = s[b.in_idx[0]]
        vq = s[b.in_idx[1]]
        id_ = s[b.in_idx[2]]
        iq = s[b.in_idx[3]]
        vdc = s[b.in_idx[4]]
        eff = b.param[0] if len(b.param) >= 1 else 1.0
        idc = 0.0
        if vdc > 1e-3:
            idc = 1.5 * (vd * id_ + vq * iq) / (vdc * eff)
        o = s[b.out_idx[0]]
        nv = o + w * (idc - o)
        s[b.out_idx[0]] = nv
        maxdelta[0] = max(maxdelta[0], abs(nv - o))
    elif b.op == ops.OP_VCVS:
        vp, vn = s[b.in_idx[0]], s[b.in_idx[1]]
        gain = b.param[0] if b.param else 1.0
        y = gain * (vp - vn)
        o = s[b.out_idx[0]]
        s[b.out_idx[0]] = o + w * (y - o)
        maxdelta[0] = max(maxdelta[0], abs(s[b.out_idx[0]] - o))
    elif b.op == ops.OP_VCCS:
        vp, vn = s[b.in_idx[0]], s[b.in_idx[1]]
        gm = b.param[0] if b.param else 1.0
        y = gm * (vp - vn)
        o = s[b.out_idx[0]]
        s[b.out_idx[0]] = o + w * (y - o)
        maxdelta[0] = max(maxdelta[0], abs(s[b.out_idx[0]] - o))
def _fuentes(bloques, ctx):
    s, t = ctx.sig, ctx.t
    for b in bloques:
        op = b.op
        if op == ops.OP_SRC_CONST:
            s[b.out_idx[0]] = b.param[0]
        elif op == ops.OP_SRC_STEP:
            vf, ts, v0 = b.param
            s[b.out_idx[0]] = vf if t >= ts else v0
        elif op == ops.OP_SRC_RAMP:
            k, t0, off = b.param
            s[b.out_idx[0]] = off + k * (t - t0) if t >= t0 else off
        elif op == ops.OP_SRC_SIN:
            Am, f, ph, off = b.param
            s[b.out_idx[0]] = off + Am * np.sin(2 * np.pi * f * t + ph)
        elif op == ops.OP_SRC_TRIF:
            Am, f, ph = b.param
            w = 2 * np.pi * f * t
            s[b.out_idx[0]] = Am * np.sin(w + ph)
            s[b.out_idx[1]] = Am * np.sin(w + ph - 2 * np.pi / 3)
            s[b.out_idx[2]] = Am * np.sin(w + ph + 2 * np.pi / 3)
        elif op == ops.OP_PULSO_RECT:
            amp, T, duty, ph, off = b.param
            tt = (t + ph) % T
            if tt < 0.0:
                tt += T
            s[b.out_idx[0]] = off + (amp if tt < duty * T else 0.0)
        elif op == ops.OP_SRC_CSV:
            n = int(b.param[0])
            interp = b.param[1] > 0.5
            tt = b.param[2:2 + n]
            yy = b.param[2 + n:]
            if t <= tt[0]:
                y = yy[0]
            elif t >= tt[-1]:
                y = yy[-1]
            else:
                i = _busqueda_bp(tt, t)
                if interp:
                    f = (t - tt[i]) / (tt[i + 1] - tt[i])
                    y = yy[i] + f * (yy[i + 1] - yy[i])
                else:
                    y = yy[i]
            s[b.out_idx[0]] = y
        elif op == ops.OP_SRC_TABLE:
            n = int(b.param[0])
            interp = b.param[1] > 0.5
            tt = b.param[2:2 + n]
            yy = b.param[2 + n:]
            if t <= tt[0]:
                y = yy[0]
            elif t >= tt[-1]:
                y = yy[-1]
            else:
                i = _busqueda_bp(tt, t)
                if interp:
                    f = (t - tt[i]) / (tt[i + 1] - tt[i])
                    y = yy[i] + f * (yy[i + 1] - yy[i])
                else:
                    y = yy[i]
            s[b.out_idx[0]] = y
def _relay_update(b, ctx, x):
    u = ctx.sig[b.in_idx[0]]
    on, off, y_on, y_off = b.param[:4]
    y = x[0]
    if u >= on:
        y = y_on
    elif u <= off:
        y = y_off
    x[0] = y
    ctx.sig[b.out_idx[0]] = y
def _diodo_update(b, ctx, x):
    va, vc = ctx.sig[b.in_idx[0]], ctx.sig[b.in_idx[1]]
    r_on, r_off, vf, h = b.param[:4]
    v = va - vc
    on = x[0] > 0.5
    if v >= vf + h:
        on = True
    elif v <= vf - h:
        on = False
    x[0] = 1.0 if on else 0.0
    r = r_on if on else r_off
    ctx.sig[b.out_idx[0]] = v / r
    ctx.sig[b.out_idx[1]] = v
def _retenedor_update(b, ctx, x):
    umb = b.param[0]
    trig = ctx.sig[b.in_idx[1]]
    prev = x[1]
    if trig > umb and prev <= umb:
        x[0] = ctx.sig[b.in_idx[0]]
    x[1] = trig
    ctx.sig[b.out_idx[0]] = x[0]
def _maq_estados_update(b, ctx, x):
    n_trans = int(b.param[1])
    tr = b.param[2:]
    estado = int(x[0])
    for k in range(n_trans):
        from_, to, sig, opc, umb = tr[5 * k:5 * k + 5]
        if int(from_) != estado:
            continue
        u = ctx.sig[b.in_idx[int(sig)]]
        opc = int(opc)
        if opc == 0:
            ok = u < umb
        elif opc == 1:
            ok = u <= umb
        elif opc == 2:
            ok = u > umb
        elif opc == 3:
            ok = u >= umb
        elif opc == 4:
            ok = abs(u - umb) <= 1e-12
        else:
            ok = abs(u - umb) > 1e-12
        if ok:
            estado = int(to)
            break
    x[0] = estado
    ctx.sig[b.out_idx[0]] = float(estado)
def _pll_update(b, ctx, x):
    h = ctx.dt
    Kp, Ki, wff = b.param[0], b.param[1], b.param[2]
    th, xi = x[0], x[1]
    al, be = _clarke(ctx.sig[b.in_idx[0]], ctx.sig[b.in_idx[1]],
                     ctx.sig[b.in_idx[2]])
    d, q = _park(al, be, th)
    w = wff + Kp * q + Ki * xi
    x[1] = xi + h * q
    x[0] = th + h * w
    ctx.sig[b.out_idx[0]] = w
    ctx.sig[b.out_idx[1]] = x[0]
def _pmac_lut(p, ids):
    if len(p) < 9:
        return p[1] * ids + p[3], p[1]
    n = int(p[8])
    t = p[9:9 + 2 * n]
    if n < 2:
        return p[1] * ids + p[3], p[1]
    if ids <= t[0]:
        s = (t[3] - t[1]) / (t[2] - t[0])
        return t[1] + s * (ids - t[0]), s
    for k in range(n - 1):
        ia, fa, ib, fb = t[2 * k], t[2 * k + 1], t[2 * k + 2], t[2 * k + 3]
        if ids <= ib:
            s = (fb - fa) / (ib - ia)
            return fa + s * (ids - ia), s
    ia, fa = t[2 * (n - 2)], t[2 * (n - 2) + 1]
    ib, fb = t[2 * (n - 1)], t[2 * (n - 1) + 1]
    s = (fb - fa) / (ib - ia)
    return fb + s * (ids - ib), s
def _maq_deriv(b, ctx, x, dx):
    s = ctx.sig
    p = b.param
    in0 = b.in_idx
    if b.op == ops.OP_MAQ_PMAC:
        rs, Ld, Lq, lam, P, J, Bm = p[:7]
        ext = p[7] > 0.5
        iqs, ids = x[0], x[1]
        wm = s[in0[3]] if ext else x[2]
        the = (P / 2) * s[in0[4]] if ext else x[3]
        we = (P / 2) * wm
        al, be = _clarke(s[in0[0]], s[in0[1]], s[in0[2]])
        vd, vq = _park(al, be, the)
        lam_d, Lds = _pmac_lut(p, ids)
        dx[0] = (vq - rs * iqs - we * lam_d) / Lq
        dx[1] = (vd - rs * ids + we * Lq * iqs) / Lds
        Te = 1.5 * (P / 2) * (lam_d * iqs - Lq * iqs * ids)
        if not ext:
            dx[2] = (Te - s[in0[3]] - Bm * wm) / J
            dx[3] = we
    elif b.op == ops.OP_MAQ_INDUCCION:
        rs, rr, Li00, Li01, Li11, wf, P, J, Bm, ext = p
        lqs, lds, lqr, ldr = x[0], x[1], x[2], x[3]
        wm = s[in0[3]] if ext else x[4]
        we = (P / 2) * wm
        al, be = _clarke(s[in0[0]], s[in0[1]], s[in0[2]])
        vd, vq = _park(al, be, 0.0)
        iqs = Li00 * lqs + Li01 * lqr
        iqr = Li01 * lqs + Li11 * lqr
        ids = Li00 * lds + Li01 * ldr
        idr = Li01 * lds + Li11 * ldr
        dx[0] = vq - rs * iqs - wf * lds
        dx[1] = vd - rs * ids + wf * lqs
        dx[2] = -rr * iqr - (wf - we) * ldr
        dx[3] = -rr * idr + (wf - we) * lqr
        Te = 1.5 * (P / 2) * (lds * iqs - lqs * ids)
        if not ext:
            dx[4] = (Te - s[in0[3]] - Bm * wm) / J
            dx[5] = wm
    elif b.op == ops.OP_MAQ_SINCRONA:
        rs, rfd, rkq1, rkq2, rkd, P, J, Bm = p[:8]
        Liq = p[8:17]
        Lid = p[17:26]
        ext = len(p) > 26 and p[26] > 0.5
        lqs, lkq1, lkq2, lds, lfd, lkd = x[:6]
        wm = s[in0[4]] if ext else x[6]
        the = (P / 2) * s[in0[5]] if ext else x[7]
        wr = (P / 2) * wm
        al, be = _clarke(s[in0[0]], s[in0[1]], s[in0[2]])
        vd, vq = _park(al, be, the)
        iqs = Liq[0] * lqs + Liq[1] * lkq1 + Liq[2] * lkq2
        ikq1 = Liq[3] * lqs + Liq[4] * lkq1 + Liq[5] * lkq2
        ikq2 = Liq[6] * lqs + Liq[7] * lkq1 + Liq[8] * lkq2
        ids = Lid[0] * lds + Lid[1] * lfd + Lid[2] * lkd
        ifd = Lid[3] * lds + Lid[4] * lfd + Lid[5] * lkd
        ikd = Lid[6] * lds + Lid[7] * lfd + Lid[8] * lkd
        dx[0] = vq - rs * iqs - wr * lds
        dx[1] = -rkq1 * ikq1
        dx[2] = -rkq2 * ikq2
        dx[3] = vd - rs * ids + wr * lqs
        dx[4] = s[in0[3]] - rfd * ifd
        dx[5] = -rkd * ikd
        Te = 1.5 * (P / 2) * (lds * iqs - lqs * ids)
        if not ext:
            dx[6] = (Te - s[in0[4]] - Bm * wm) / J
            dx[7] = wr
    elif b.op == ops.OP_MAQ_CC:
        ra, La, rf, Lf, LAF, J, Bm, ext = p
        ia, i_f = x[0], x[1]
        wm = s[in0[2]] if ext else x[2]
        dx[0] = (s[in0[0]] - ra * ia - LAF * i_f * wm) / La
        dx[1] = (s[in0[1]] - rf * i_f) / Lf
        Te = LAF * i_f * ia
        if not ext:
            dx[2] = (Te - s[in0[2]] - Bm * wm) / J
            dx[3] = wm
    elif b.op == ops.OP_MAQ_DC_PM:
        ra, La, Kt, J, Bm, ext = p
        ia = x[0]
        wm = s[in0[1]] if ext else x[1]
        dx[0] = (s[in0[0]] - ra * ia - Kt * wm) / La
        Te = Kt * ia
        if not ext:
            dx[1] = (Te - s[in0[1]] - Bm * wm) / J
            dx[2] = wm
    elif b.op == ops.OP_EJE_MECANICO:
        J_eq, Bm_eq = p[0], p[1]
        wm = x[0]
        sum_Te = sum(s[in0[i]] for i in range(b.n_in - 1))
        TL = s[in0[b.n_in - 1]]
        dx[0] = (sum_Te - TL - Bm_eq * wm) / J_eq
        dx[1] = wm
    elif b.op == ops.OP_INTEGRADOR:
        dx[0] = s[in0[0]]
    elif b.op == ops.OP_LIM_RAPIDEZ:
        up, down = p[0], p[1]
        pend = (s[in0[0]] - x[0]) / ctx.dt
        dx[0] = max(-down, min(up, pend))
    elif b.op in (ops.OP_POT_BUCK, ops.OP_POT_BOOST, ops.OP_POT_BUCKBOOST):
        L, C, R = p
        vin, d = s[in0[0]], s[in0[1]]
        iL, vC = x
        if b.op == ops.OP_POT_BUCK:
            dx[0] = (d * vin - vC) / L
            iC = iL
        else:
            dx[0] = (d * vin - (1 - d) * vC) / L if b.op == ops.OP_POT_BUCKBOOST \
                    else (vin - (1 - d) * vC) / L
            iC = (1 - d) * iL
        dx[1] = (iC - vC / R) / C
    elif b.op == ops.OP_POT_RECT_3F:
        C, R, Rint = p
        va, vb, vc = s[in0[0]], s[in0[1]], s[in0[2]]
        vrec = max(va, vb, vc) - min(va, vb, vc)
        ich = max((vrec - x[0]) / Rint, 0.0)
        dx[0] = (ich - x[0] / R) / C
    elif b.op == ops.OP_POT_INV_3F:
        f, fsw, m0, m1, tramp, Lf, Cf, R = p[:8]
        conmutada = int(p[8])
        vdc = s[in0[0]]
        if tramp > 0.0:
            mt = m0 + (m1 - m0) * min(ctx.t / tramp, 1.0)
        else:
            mt = m1
        mt = max(mt, 0.0)
        Tsw = 1.0 / fsw
        frac = (ctx.t % Tsw) / Tsw
        carr = -1.0 + 4.0 * frac if frac < 0.5 else 3.0 - 4.0 * frac
        for k in range(3):
            th = 2 * np.pi * f * ctx.t - (2 * np.pi / 3.0) * k
            ref = mt * np.sin(th)
            if conmutada:
                vp = 0.5 * vdc * (1.0 if ref > carr else -1.0)
            else:
                vp = 0.5 * vdc * ref
            dx[k] = (vp - x[3 + k]) / Lf
            dx[3 + k] = (x[k] - x[3 + k] / R) / Cf
    elif b.op == ops.OP_POT_INV_1F:
        f, fsw, m0, m1, tramp, Lf, Cf, R = p[:8]
        conmutada = int(p[8])
        vdc = s[in0[0]]
        if tramp > 0.0:
            mt = m0 + (m1 - m0) * min(ctx.t / tramp, 1.0)
        else:
            mt = m1
        mt = max(mt, 0.0)
        ref = mt * np.sin(2 * np.pi * f * ctx.t)
        if conmutada:
            Tsw = 1.0 / fsw
            frac = (ctx.t % Tsw) / Tsw
            carr = -1.0 + 4.0 * frac if frac < 0.5 else 3.0 - 4.0 * frac
            vp = 0.5 * vdc * (1.0 if ref > carr else -1.0)
        else:
            vp = 0.5 * vdc * ref
        dx[0] = (vp - x[1]) / Lf
        dx[1] = (x[0] - x[1] / R) / Cf
    elif b.op == ops.OP_CARGA_RL_3F:
        R, L = p[0], p[1]
        va, vb, vc = s[in0[0]], s[in0[1]], s[in0[2]]
        vn = (va + vb + vc) / 3.0
        dx[0] = (va - vn - R * x[0]) / L
        dx[1] = (vb - vn - R * x[1]) / L
    elif b.op == ops.OP_BATERIA:
        E0, K, Q, A, B, R, tau = p[:7]
        eta_c = p[8]
        i = s[in0[0]]
        di = i / 3600.0
        if i < 0.0:
            di = i * eta_c / 3600.0
        if (x[0] >= 0.9 * Q and di > 0.0) or (x[0] <= 0.1 * Q and di < 0.0):
            di = 0.0
        dx[0] = di
        dx[1] = (i - x[1]) / tau
        if p[9] > 0.5:
            dx[2] = B * abs(i) * (-x[2] + A)
        else:
            dx[2] = 0.0
        if len(p) >= 14 and p[10] > 0.0:
            R_eff = R * (1.0 + p[12] * (x[3] - 25.0))
            dx[3] = (i * i * R_eff - (x[3] - p[13]) / p[11]) / p[10]
        else:
            dx[3] = 0.0
    elif b.op == ops.OP_BATERIA_ECM:
        Q_nom = p[0]
        R0, R1, C1, R2, C2 = p[4], p[5], p[6], p[7], p[8]
        Np = p[10]
        n_ocv = int(p[11])
        ocv_soc = p[12:12 + 2*n_ocv]
        ocv_idx = 12 + 2*n_ocv
        R_th_pack = p[ocv_idx+6]
        C_th_pack = p[ocv_idx+7]
        T_amb = p[ocv_idx+8]
        i = s[in0[0]]
        soc, V_rc1, V_rc2, T_cell = x[0], x[1], x[2], x[3]
        i_cell = i / Np
        eta_c = 0.99
        di_dt = -i_cell / 3600.0
        if i < 0.0:
            di_dt = -i_cell * eta_c / 3600.0
        dx[0] = di_dt
        if (soc <= 0.0 and dx[0] < 0.0) or (soc >= 1.0 and dx[0] > 0.0):
            dx[0] = 0.0
        ocv = np.interp(soc, p[12:12+2*int(p[11]):2], p[13:12+2*int(p[11]):2])
        dx[1] = (i_cell - V_rc1 / p[5]) / p[6]
        dx[2] = (i_cell - V_rc2 / p[7]) / p[8]
        R0_eff = p[4] * (1.0 + 0.003 * (T_cell - 25.0))
        P_loss = i * i * R0_eff * p[9] * Np
        dx[3] = (P_loss - (T_cell - p[ocv_idx+8]) / p[ocv_idx+6]) / p[ocv_idx+7]
    elif b.op == ops.OP_MASA_TERMICA:
        C = p[0]
        suma = sum(s[i] for i in in0)
        if p[2] > 0.0:
            suma -= (x[0] - p[1]) / p[2]
        dx[0] = suma / C
    elif b.op == ops.OP_VEHICULO:
        mass, Cd, A, Crr, rho, g = p[0], p[1], p[2], p[3], p[4], p[5]
        gr, r_w = p[6], p[7]
        eff_g, eff_r = p[8], p[9]
        T_mot = s[in0[1]]
        grade = s[in0[2]] if b.n_in >= 3 and in0[2] >= 0 else 0.0
        v = x[0]
        eff = eff_g if T_mot >= 0.0 else eff_r
        F_drive = T_mot * gr * eff / r_w
        F_aero = 0.5 * rho * Cd * A * v * abs(v)
        sign_v = 1.0 if v > 1e-3 else (-1.0 if v < -1e-3 else 0.0)
        F_rod = mass * g * (Crr * sign_v * np.cos(grade) + np.sin(grade))
        dx[0] = (F_drive - F_aero - F_rod) / mass
        if v <= 1e-9 and v >= -1e-9 and dx[0] < 0.0 and F_drive >= -1e-9 and abs(grade) < 1e-9:
            dx[0] = 0.0
    elif b.op == ops.OP_INDUCTOR:
        dx[0] = s[in0[0]] / p[0]
    elif b.op == ops.OP_EJE_FLEXIBLE:
        dx[0] = s[in0[0]]
        dx[1] = s[in0[1]]
    elif b.op == ops.OP_MNA:
        pass
    elif b.op == ops.OP_MUTUAL_INDUCTOR:
        L1, L2, M = p[0], p[1], p[2]
        v1 = s[b.in_idx[0]]
        v2 = s[b.in_idx[1]]
        det = L1 * L2 - M * M
        if abs(det) > 1e-12:
            dx[0] = (L2 * x[0] - M * x[1] + v1 * L2 - v2 * M) / det
            dx[1] = (-M * x[0] + L1 * x[1] + v2 * L1 - v1 * M) / det
        else:
            dx[0] = 0.0
            dx[1] = 0.0
def _maq_out(b, ctx, x):
    s = ctx.sig
    p = b.param
    out = b.out_idx
    in0 = b.in_idx
    if b.op == ops.OP_MAQ_PMAC:
        Ld, Lq, lam, P = p[1], p[2], p[3], p[4]
        ext = p[7] > 0.5
        iqs, ids = x[0], x[1]
        wm = s[in0[3]] if ext else x[2]
        the = (P / 2) * s[in0[4]] if ext else x[3]
        al, be = _inv_park(ids, iqs, the)
        s[out[0]], s[out[1]], s[out[2]] = _inv_clarke(al, be)
        s[out[3]], s[out[4]] = iqs, ids
        s[out[5]], s[out[6]] = wm, the / (P / 2)
        s[out[7]] = the
        lam_d, _ = _pmac_lut(p, ids)
        s[out[8]] = 1.5 * (P / 2) * (lam_d * iqs - Lq * iqs * ids)
    elif b.op == ops.OP_MAQ_INDUCCION:
        Li00, Li01, Li11, P = p[2], p[3], p[4], p[6]
        ext = p[9] > 0.5
        lqs, lds, lqr, ldr = x[:4]
        wm = s[in0[3]] if ext else x[4]
        thr = s[in0[4]] if ext else x[5]
        iqs = Li00 * lqs + Li01 * lqr
        ids = Li00 * lds + Li01 * ldr
        s[out[0]], s[out[1]], s[out[2]] = _inv_clarke(ids, iqs)
        s[out[3]], s[out[4]] = iqs, ids
        s[out[5]], s[out[6]] = wm, thr
        s[out[7]] = 0.0
        s[out[8]] = 1.5 * (P / 2) * (lds * iqs - lqs * ids)
    elif b.op == ops.OP_MAQ_SINCRONA:
        Liq, Lid, P = p[8:17], p[17:26], p[5]
        ext = len(p) > 26 and p[26] > 0.5
        lqs, lkq1, lkq2, lds, lfd, lkd = x[:6]
        wm = s[in0[4]] if ext else x[6]
        the = (P / 2) * s[in0[5]] if ext else x[7]
        iqs = Liq[0] * lqs + Liq[1] * lkq1 + Liq[2] * lkq2
        ids = Lid[0] * lds + Lid[1] * lfd + Lid[2] * lkd
        al, be = _inv_park(ids, iqs, the)
        s[out[0]], s[out[1]], s[out[2]] = _inv_clarke(al, be)
        s[out[3]], s[out[4]] = iqs, ids
        s[out[5]], s[out[6]] = wm, the / (P / 2)
        s[out[7]] = the
        s[out[8]] = 1.5 * (P / 2) * (lds * iqs - lqs * ids)
    elif b.op == ops.OP_MAQ_CC:
        ra, LAF = p[0], p[4]
        ext = p[7] > 0.5
        wm = s[in0[2]] if ext else x[2]
        s[out[0]], s[out[1]] = x[0], x[1]
        s[out[2]] = wm
        s[out[3]] = s[in0[3]] if ext else x[3]
        s[out[4]] = LAF * x[1] * x[0]
        e = LAF * x[1] * wm
        s[out[5]] = e
        s[out[6]] = e + ra * x[0]
    elif b.op == ops.OP_MAQ_DC_PM:
        ra, Kt = p[0], p[2]
        ext = p[5] > 0.5
        wm = s[in0[1]] if ext else x[1]
        s[out[0]] = x[0]
        s[out[1]] = wm
        s[out[2]] = s[in0[2]] if ext else x[2]
        s[out[3]] = Kt * x[0]
        e = Kt * wm
        s[out[4]] = e
        s[out[5]] = e + ra * x[0]
    elif b.op == ops.OP_EJE_MECANICO:
        s[out[0]] = x[0]
        s[out[1]] = x[1]
    elif b.op == ops.OP_INTEGRADOR:
        s[out[0]] = x[0]
    elif b.op == ops.OP_LIM_RAPIDEZ:
        s[out[0]] = x[0]
    elif b.op in (ops.OP_POT_BUCK, ops.OP_POT_BOOST, ops.OP_POT_BUCKBOOST):
        s[out[0]] = x[1]
        s[out[1]] = x[0]
    elif b.op == ops.OP_POT_RECT_3F:
        Rint = p[2]
        va, vb, vc = s[in0[0]], s[in0[1]], s[in0[2]]
        vrec = max(va, vb, vc) - min(va, vb, vc)
        ich = max((vrec - x[0]) / Rint, 0.0)
        s[out[0]] = x[0]
        s[out[1]] = ich
    elif b.op == ops.OP_POT_INV_3F:
        for k in range(3):
            s[out[k]] = x[3 + k]
            s[out[3 + k]] = x[k]
    elif b.op == ops.OP_POT_INV_1F:
        s[out[0]] = x[1]
        s[out[1]] = x[0]
    elif b.op == ops.OP_CARGA_RL_3F:
        s[out[0]] = x[0]
        s[out[1]] = x[1]
        s[out[2]] = -(x[0] + x[1])
    elif b.op == ops.OP_BATERIA:
        E0, K, Q, A, B, Vcap = p[0], p[1], p[2], p[3], p[4], p[7]
        R = p[5]
        if len(p) >= 14 and p[10] > 0.0:
            R = R * (1.0 + p[12] * (x[3] - 25.0))
        qmax = 0.9 * Q
        it = min(max(x[0], 0.0), qmax)
        den1 = Q - it
        den2 = max(it - 0.1 * Q, 1e-9 * Q)
        ifi = x[1]
        exp_h = x[2] if p[9] > 0.5 else A * np.exp(-B * it)
        if ifi >= 0.0:
            E = E0 - K * Q * (it + ifi) / den1 + exp_h
        else:
            E = E0 - K * Q * ifi / den2 - K * Q * it / den1 + exp_h
            E = min(E, Vcap)
        s[out[0]] = E - R * s[in0[0]]
        s[out[1]] = 1.0 - it / Q
        s[out[2]] = x[3]
    elif b.op == ops.OP_BATERIA_ECM:
        Q_nom, V_nom, V_min, V_max = p[0], p[1], p[2], p[3]
        R0, R1, C1, R2, C2 = p[4], p[5], p[6], p[7], p[8]
        Ns, Np = p[9], p[10]
        n_ocv = int(p[11])
        ocv_soc = p[12:12 + 2*n_ocv]
        ocv_idx = 12 + 2*n_ocv
        I_chg_cont, I_dch_cont = p[ocv_idx], p[ocv_idx+1]
        T_min_chg, T_max_chg = p[ocv_idx+2], p[ocv_idx+3]
        T_min_dch, T_max_dch = p[ocv_idx+4], p[ocv_idx+5]
        R_th_pack, C_th_pack, T_amb = p[ocv_idx+6], p[ocv_idx+7], p[ocv_idx+8]
        i = s[in0[0]]
        soc, V_rc1, V_rc2, T_cell = x[0], x[1], x[2], x[3]
        i_cell = i / Np
        ocv = np.interp(soc, ocv_soc[::2], ocv_soc[1::2])
        R0_eff = R0 * (1.0 + 0.003 * (T_cell - 25.0))
        V_term_cell = ocv - V_rc1 - V_rc2 - R0_eff * i_cell
        V_term_pack = V_term_cell * p[9]
        V_min_pack = V_min * Ns
        V_max_pack = V_max * Ns
        V_term_pack = np.clip(V_term_pack, V_min_pack, V_max_pack)
        R0_eff = R0 * (1.0 + 0.003 * (T_cell - 25.0))
        P_loss = i * i * R0_eff * Ns * Np
        I_chg_lim = p[ocv_idx] * Np
        I_dch_lim = p[ocv_idx+1] * Np
        T_min_chg = p[ocv_idx+2]
        T_max_chg = p[ocv_idx+3]
        T_min_dch = p[ocv_idx+4]
        T_max_dch = p[ocv_idx+5]
        if T_cell < T_min_chg or T_cell > T_max_chg:
            I_chg_lim = 0.0
        if T_cell < T_min_dch or T_cell > T_max_dch:
            I_dch_lim = 0.0
        s[out[0]] = V_term_pack
        s[out[1]] = soc
        s[out[2]] = T_cell
        s[out[3]] = P_loss
        s[out[4]] = I_chg_lim
        s[out[5]] = I_dch_lim
    elif b.op == ops.OP_MASA_TERMICA:
        s[out[0]] = x[0]
    elif b.op == ops.OP_VEHICULO:
        mass, Cd, A, Crr, rho, g = p[0], p[1], p[2], p[3], p[4], p[5]
        gr, r_w = p[6], p[7]
        eff = p[8]
        v = x[0]
        omega_m = s[in0[0]]
        grade = s[in0[2]] if b.n_in >= 3 and in0[2] >= 0 else 0.0
        F_aero = 0.5 * rho * Cd * A * v * v
        F_rod = mass * g * (Crr * np.cos(grade) + np.sin(grade))
        s[out[0]] = (F_aero + F_rod) * r_w / (gr * eff)
        s[out[1]] = v
        s[out[2]] = v * 3.6
        s[out[3]] = omega_m
        v_ref = omega_m * r_w / gr
        s[out[4]] = (0.5 * rho * Cd * A * v_ref * v_ref + F_rod) \
                    * r_w / (gr * eff)
        s[out[5]] = grade
        s[out[6]] = v * gr / r_w
    elif b.op == ops.OP_INDUCTOR:
        s[out[0]] = x[0]
    elif b.op == ops.OP_CAPACITOR:
        s[out[0]] = x[0]
    elif b.op == ops.OP_EJE_FLEXIBLE:
        K, B = p[0], p[1]
        s[out[0]] = K * (x[0] - x[1]) + B * (s[in0[0]] - s[in0[1]])
    elif b.op == ops.OP_VCVS:
        gain = b.param[0] if b.param else 1.0
        vp = s[b.in_idx[0]]
        vn = s[b.in_idx[1]]
        s[b.out_idx[0]] = gain * (vp - vn)
    elif b.op == ops.OP_VCCS:
        gm = b.param[0] if b.param else 1.0
        vp = s[b.in_idx[0]]
        vn = s[b.in_idx[1]]
        s[b.out_idx[0]] = gm * (vp - vn)
    elif b.op == ops.OP_MUTUAL_INDUCTOR:
        s[b.out_idx[0]] = x[0]
        s[b.out_idx[1]] = x[1]
    elif b.op == ops.OP_MNA:
        pass
def simular(bloques, dt, t_fin, rec_idx, metodo=0, max_iter=50, tol=1e-9,
            w_opt=1.0, orden_estatico=None):
    n_steps = int(round(t_fin / dt)) + 1
    n_sig = max(
        (max(b.out_idx) + 1 if b.out_idx else 0) for b in bloques
    ) if bloques else 0
    n_sig = max(n_sig, max(rec_idx) + 1 if rec_idx else 0)
    ctx = _Ctx(bloques, dt, metodo, max_iter, tol, w_opt, n_sig)
    est = ctx.est
    idx = {id(b): i for i, b in enumerate(bloques)}
    estaticos = [b for b in bloques if b.op in ops.ES_ESTATICO]
    if orden_estatico is not None:
        estaticos = [bloques[k] for k in orden_estatico]
    dinam = [b for b in bloques if b.op in ops.ES_DINAMICO]
    for b in bloques:
        b.in_idx = list(b.in_idx)
    _fuentes(bloques, ctx)
    for i, b in enumerate(bloques):
        if b.op == ops.OP_RELAY:
            _relay_update(b, ctx, est[i])
        elif b.op == ops.OP_DIODO:
            _diodo_update(b, ctx, est[i])
        elif b.op == ops.OP_RETENEDOR:
            ctx.sig[b.out_idx[0]] = est[i][0]
        elif b.op == ops.OP_MAQ_ESTADOS:
            ctx.sig[b.out_idx[0]] = est[i][0]
        elif b.op in (ops.OP_MAQ_PMAC, ops.OP_MAQ_INDUCCION,
                      ops.OP_MAQ_SINCRONA, ops.OP_MAQ_CC, ops.OP_MAQ_DC_PM,
                      ops.OP_EJE_MECANICO, ops.OP_INTEGRADOR,
                      ops.OP_POT_BUCK, ops.OP_POT_BOOST, ops.OP_POT_BUCKBOOST,
                      ops.OP_POT_RECT_3F, ops.OP_POT_INV_3F, ops.OP_POT_INV_1F,
                      ops.OP_CARGA_RL_3F, ops.OP_BATERIA, ops.OP_BATERIA_ECM,
                      ops.OP_LIM_RAPIDEZ,
                      ops.OP_MASA_TERMICA, ops.OP_EJE_FLEXIBLE,
                      ops.OP_VEHICULO,
                      ops.OP_INDUCTOR, ops.OP_CAPACITOR,
                      ops.OP_MUTUAL_INDUCTOR):
            _maq_out(b, ctx, est[i])
        elif b.op == ops.OP_PLL:
            Kp, Ki, wff = b.param[0], b.param[1], b.param[2]
            al, be = _clarke(ctx.sig[b.in_idx[0]], ctx.sig[b.in_idx[1]],
                             ctx.sig[b.in_idx[2]])
            d, q = _park(al, be, est[i][0])
            ctx.sig[b.out_idx[0]] = wff + Kp * q + Ki * est[i][1]
            ctx.sig[b.out_idx[1]] = est[i][0]
    rec = np.zeros((len(rec_idx), n_steps))
    maxd = [0.0]
    for _ in range(ctx.max_iter):
        maxd[0] = 0.0
        for b in estaticos:
            _eval_estatico(b, ctx, maxd)
        if maxd[0] < ctx.tol:
            break
    for r, sig_i in enumerate(rec_idx):
        rec[r, 0] = ctx.sig[sig_i]
    for s in range(1, n_steps):
        ctx.t += dt
        _fuentes(bloques, ctx)
        maxd = [0.0]
        for _ in range(ctx.max_iter):
            maxd[0] = 0.0
            for b in estaticos:
                _eval_estatico(b, ctx, maxd)
            if maxd[0] < ctx.tol:
                break
        if maxd[0] >= ctx.tol:
            raise RuntimeError(
                "El lazo algebraico no convergio "
                f"(t = {ctx.t:.6g} s, max_iter = {ctx.max_iter}, "
                f"tol = {ctx.tol:g}). Revisa el modelo o ajusta "
                "Modelo(max_iter=..., tol=..., w_opt=...)."
            )
        for i, b in enumerate(bloques):
            if b.op not in ops.ES_DINAMICO:
                continue
            x = est[idx[id(b)]]
            if b.op == ops.OP_TF:
                n = int(b.param[0])
                bd = b.param[1:1 + n + 1]
                ad = b.param[1 + n + 1:]
                u = ctx.sig[b.in_idx[0]]
                y = bd[0] * u
                for k in range(1, n + 1):
                    y += bd[k] * x[k - 1] - ad[k - 1] * x[n + k - 1]
                for k in range(n - 1, 0, -1):
                    x[k] = x[k - 1]
                x[0] = u
                for k in range(2 * n - 1, n, -1):
                    x[k] = x[k - 1]
                x[n] = y
                ctx.sig[b.out_idx[0]] = y
            elif b.op == ops.OP_PID:
                Kp, Ki, Kd, Tf, umin, umax = b.param
                e = ctx.sig[b.in_idx[0]]
                ud = (Kd * (e - x[1]) + Tf * x[2]) / (Tf + ctx.dt)
                u = Kp * e + Ki * x[0] + ud
                if not (u > umax and e > 0) and not (u < umin and e < 0):
                    x[0] += ctx.dt * e
                x[1], x[2] = e, ud
                u = max(umin, min(umax, u))
                ctx.sig[b.out_idx[0]] = u
            elif b.op == ops.OP_RELAY:
                _relay_update(b, ctx, x)
            elif b.op == ops.OP_DIODO:
                _diodo_update(b, ctx, x)
            elif b.op == ops.OP_RETENEDOR:
                _retenedor_update(b, ctx, x)
            elif b.op == ops.OP_MAQ_ESTADOS:
                _maq_estados_update(b, ctx, x)
            elif b.op in (ops.OP_MAQ_PMAC, ops.OP_MAQ_INDUCCION,
                          ops.OP_MAQ_SINCRONA, ops.OP_MAQ_CC, ops.OP_MAQ_DC_PM,
                          ops.OP_EJE_MECANICO, ops.OP_INTEGRADOR,
                          ops.OP_POT_BUCK, ops.OP_POT_BOOST, ops.OP_POT_BUCKBOOST,
                          ops.OP_POT_RECT_3F, ops.OP_POT_INV_3F, ops.OP_POT_INV_1F,
                          ops.OP_CARGA_RL_3F, ops.OP_BATERIA, ops.OP_BATERIA_ECM,
                          ops.OP_LIM_RAPIDEZ,
                          ops.OP_MASA_TERMICA, ops.OP_EJE_FLEXIBLE,
                          ops.OP_VEHICULO,
                          ops.OP_INDUCTOR, ops.OP_CAPACITOR):
                dx = np.zeros(b.n_state)
                _maq_deriv(b, ctx, x, dx)
                if ctx.metodo == 0:
                    x += ctx.dt * dx
                else:
                    x2 = x + 0.5 * ctx.dt * dx
                    k2 = np.zeros(b.n_state)
                    _maq_deriv(b, ctx, x2, k2)
                    x3 = x + 0.5 * ctx.dt * k2
                    k3 = np.zeros(b.n_state)
                    _maq_deriv(b, ctx, x3, k3)
                    x4 = x + ctx.dt * k3
                    k4 = np.zeros(b.n_state)
                    _maq_deriv(b, ctx, x4, k4)
                    x += (ctx.dt / 6) * (dx + 2 * k2 + 2 * k3 + k4)
                _maq_out(b, ctx, x)
            elif b.op == ops.OP_PLL:
                _pll_update(b, ctx, x)
        for _ in range(ctx.max_iter):
            maxd[0] = 0.0
            for b in estaticos:
                _eval_estatico(b, ctx, maxd)
            if maxd[0] < ctx.tol:
                break
        if maxd[0] >= ctx.tol:
            raise RuntimeError(
                "El lazo algebraico no convergio "
                f"(t = {ctx.t:.6g} s, max_iter = {ctx.max_iter}, "
                f"tol = {ctx.tol:g}). Revisa el modelo o ajusta "
                "Modelo(max_iter=..., tol=..., w_opt=...)."
            )
        for r, sig_i in enumerate(rec_idx):
            rec[r, s] = ctx.sig[sig_i]
    return rec
