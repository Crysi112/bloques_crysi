import numpy as np

from bloques_crysi import (
    Modelo, FuenteConstante, Ganancia, InvClarke, InvPark, Sensor,
    MaquinaImanesPermanentes, Puerto, ScopeTiempoReal,
    Demultiplexor, Multiplexor, Multiplicador, MasaTermica,Scope,
)

S2 = np.sqrt(2.0)

def _arma_pmsm(nombre, rs, Ld, Lq, lam_m, P, J, vq, tl=0.0, dt=1e-5):
    m = Modelo(dt=dt, metodo="rk4")
    vd = m.add(FuenteConstante("vd", 0.0))
    vq_s = m.add(FuenteConstante("vq", vq))
    tl_s = m.add(FuenteConstante("tl", tl))
    maq = m.add(MaquinaImanesPermanentes(
        "pm", rs=rs, Ld=Ld, Lq=Lq, lam_m=lam_m, P=P, J=J, Bm=0.0))
    m.conectar(tl_s.salida, maq.T_L)
    th = Sensor("th", maq, "sal", 7, 1, canales=["th"])
    ipk = m.add(InvPark("ipk"))
    icl = m.add(InvClarke("icl"))
    m.conectar(vd.salida, Puerto(ipk, "ent", 0, 1))
    m.conectar(vq_s.salida, Puerto(ipk, "ent", 1, 1))
    m.conectar(th, Puerto(ipk, "ent", 2, 1))
    m.conectar(ipk.salida, icl.entrada)
    m.conectar(icl.salida, maq.terminales)
    return m, maq, th

def main():
    print("PMAC 4.5-1 del libro (Krause, cap. 4.5)")
    print("=" * 62)
    rs = 3.4
    lam = 60.0 / 2.0 / np.sqrt(3.0) / (2.0 * (1000.0 * 2 * np.pi / 60.0))

    vq = S2 * 11.25                 
    print(f"lambda_m = {lam:.7f} V*s  (libro: 0.0827 V*s, redondeado)")
    print(f"comprobacion 4.5-2: v_ln pico a 1000 rpm = "
          f"{lam * 2.0 * (1000.0 * 2 * np.pi / 60.0):.4f} V "
          f"= 60 V pk-pk entre lineas")

    m, maq, th = _arma_pmsm("pm451", rs=rs, Ld=12.1e-3, Lq=12.1e-3,
                            lam_m=lam, P=4, J=5e-4, vq=vq)

    sc = m.add(Scope("scope", anchos=[1, 3, 1],
                               guiones=["wm (rad/s)", "ia (A)",
                                        "T devanado (C)"]))
    m.conectar(maq.sensorVelocidad(), sc.canales[0])
    m.conectar(maq.sensor3I(), sc.canales[1])

    demux = m.add(Demultiplexor("demux_i", n_canales=3))
    m.conectar(maq.sensor3I(), demux.entrada)
    mux_ia = m.add(Multiplexor("mux_ia", n_canales=2))
    m.conectar(demux.salidas[0], mux_ia.entradas[0])
    m.conectar(demux.salidas[0], mux_ia.entradas[1])
    mult = m.add(Multiplicador("mult"))
    m.conectar(mux_ia.salida, mult.entrada)
    perd = m.add(Ganancia("perd", rs))
    m.conectar(mult.salida, perd.entrada)

    masa = m.add(MasaTermica("masa", C_th=400.0, T_inicial=25.0,
                             T_amb=25.0, R_amb=0.3))
    m.conectar(perd.salida, masa.entradas[0])
    m.conectar(masa.salida, sc.canales[2])
    res = m.run(0.6, registrar=[sc, maq.sensorVelocidad(), maq.sensor3I(),
                                masa])

    we = float(res["wm"][-1]) * 2.0  
    print(f"sin carga  : we = {we:6.1f} rad/s electricos "
          f"({res['wm'][-1]:5.1f} rad/s mecanicos, "
          f"{res['wm'][-1]*60/(2*np.pi):5.0f} rpm)")
    print(f"  libro    : ~200 rad/s electricos (eje de la Fig. 4.5-3) "
          f"= 955 rpm")
    print(f"  analitico: vq/lambda_m = {vq/lam:6.1f} rad/s")

    m2, maq2, th2 = _arma_pmsm("pm451b", rs=3.4, Ld=12.1e-3, Lq=12.1e-3,
                               lam_m=lam, P=4, J=5e-4, vq=vq, tl=0.1)
    iqs_s = Sensor("iqs", maq2, "sal", 3, 1)
    res2 = m2.run(3.0, registrar=[maq2.sensorVelocidad(), maq2.sensor3I(),
                                  maq2.sensorPar(), iqs_s])
    we2 = float(res2["wm"][-1]) * 2.0
    iqs = 0.1 / ((3.0 / 2.0) * 2.0 * lam)       

    a2 = -12.1e-3 * 12.1e-3 * iqs / 3.4
    b2 = -lam
    c2 = vq - 3.4 * iqs
    we_eq = (-b2 - np.sqrt(b2 * b2 - 4 * a2 * c2)) / (2 * a2)
    print(f"con 0.1 N.m: we = {we2:6.1f} rad/s,  "
          f"Te = {res2['Te'][-1]:5.3f} N.m,  "
          f"iqs = {res2['iqs'][-1]:5.3f} A")
    print(f"  equilibrio qd completo: we = {we_eq:6.1f} rad/s")
    print(f"  libro (modelo ideal, ids = 0): we = "
          f"{(vq - 3.4*iqs)/lam:6.1f} rad/s,  iqs = {iqs:.3f} A")
    print()
    print("PMAC 14.3 del libro (Krause, cap. 14.3)")
    print("-" * 62)
    m3 = Modelo(dt=1e-5, metodo="rk4")
    vd3 = m3.add(FuenteConstante("vd3", 0.0))
    vq3 = m3.add(FuenteConstante("vq3", 120.0 / np.sqrt(3.0)))
    g = m3.add(Ganancia("g", 0.005))
    maq3 = m3.add(MaquinaImanesPermanentes(
        "pm143", rs=2.98, Ld=11.4e-3, Lq=11.4e-3, lam_m=0.156, P=4,
        J=0.005))
    th3 = Sensor("th3", maq3, "sal", 7, 1, canales=["th3"])
    ipk3 = m3.add(InvPark("ipk3"))
    icl3 = m3.add(InvClarke("icl3"))
    m3.conectar(vd3.salida, Puerto(ipk3, "ent", 0, 1))
    m3.conectar(vq3.salida, Puerto(ipk3, "ent", 1, 1))
    m3.conectar(th3, Puerto(ipk3, "ent", 2, 1))
    m3.conectar(ipk3.salida, icl3.entrada)
    m3.conectar(icl3.salida, maq3.terminales)
    wm3 = maq3.sensorVelocidad()
    m3.conectar(wm3, g.entrada)
    m3.conectar(g.salida, maq3.T_L)
    sc3 = m3.add(Scope("scope3", anchos=[1, 1, 3],
                                 guiones=["wm (rad/s)", "Te (N.m)",
                                          "ia", "ib", "ic"]))
    m3.conectar(wm3, sc3.canales[0])
    m3.conectar(maq3.sensorPar(), sc3.canales[1])
    m3.conectar(maq3.sensor3I(), sc3.canales[2])
    res3 = m3.run(2.5, registrar=[sc3, wm3, maq3.sensor3I(),
                                  maq3.sensorPar()])
    w = float(res3["wm"][-1])
    print(f"con TL=0.005*wm: wm = {w:6.1f} rad/s mecanicos, "
          f"Te = {res3['Te'][-1]:5.3f} N.m,  TL = {0.005*w:5.3f} N.m")
    _rs, _Lq, _Ld, _lam, _vq = 2.98, 11.4e-3, 11.4e-3, 0.156, 120.0 / np.sqrt(3.0)
    w_eq = None
    for wc in np.linspace(50.0, 400.0, 7001):
        wec = 2.0 * wc
        iqs_c = (_vq - wec * _lam) / (_rs + wec * wec * _Ld * _Lq / _rs)
        if 3.0 * _lam * iqs_c <= 0.005 * wc:
            w_eq = wc
            break
    print(f"  equilibrio qd completo: wm = {w_eq:6.1f} rad/s,  "
          f"iqs = {0.005*w_eq/(3.0*0.156):6.3f} A "
          f"(simulacion reproducida por la maquina)")
    iqs_eq = (_vq - 2.0 * w_eq * _lam) / (_rs + (2.0 * w_eq)**2 * _Ld * _Lq / _rs)
    ids_eq = 2.0 * w_eq * _Lq * iqs_eq / _rs
    print(f"  exactos             : wm = {w_eq:.4f} rad/s, "
          f"we = {2.0*w_eq:.4f} rad/s, iqs = {iqs_eq:.4f} A, "
          f"ids = {ids_eq:.4f} A")
    print(f"  libro (modelo ideal, ids = 0): wm = "
          f"{120.0/np.sqrt(3.0)/(0.156*2.0 + 2.98*0.005/(3.0*0.156)):6.1f} "
          f"rad/s,  iqs = "
          f"{0.005*120.0/np.sqrt(3.0)/(0.156*2.0 + 2.98*0.005/(3.0*0.156))/(3.0*0.156):4.2f} "
          f"A  (TL = 1 N.m a 200 rad/s)")
    print()
    print("(cerrar las figuras de los Scopes para terminar)")

main()
