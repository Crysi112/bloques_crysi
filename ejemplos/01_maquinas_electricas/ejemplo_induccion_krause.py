import numpy as np

from bloques_crysi import (
    Modelo, FuenteTrifasica, FuenteConstante, MaquinaInduccion, Puerto,
    Scope, Demultiplexor, Multiplexor, Multiplicador, Ganancia, MasaTermica,
)

WB = 2 * np.pi * 60.0

def _tb(hp):
    return hp * 745.699872 / (2 * np.pi * 60.0 / 2.0)

def _armar(dt, hp, TL, bloqueado=False):
    if hp == 3:
        Vll, rs, rr, Lm, Ll, J = 220.0, 0.435, 0.816, 26.13 / WB, \
            0.754 / WB, 0.089

    elif hp == 50:
        Vll, rs, rr, Lm, Ll, J = 460.0, 0.087, 0.228, 13.08 / WB, \
            0.302 / WB, 1.662
    else:
        raise ValueError(hp)

    m = Modelo(dt=dt, metodo="rk4")
    red = m.add(FuenteTrifasica("red", amplitud=Vll * np.sqrt(2.0 / 3.0),
                                frecuencia=60.0, fase=np.pi / 2.0))
    maq = m.add(MaquinaInduccion(
        "mi", rs=rs, rr=rr, Lm=Lm, Lls=Ll, Llr=Ll, P=4, J=J, Bm=0.0,
        mecanica_interna=not bloqueado))
    m.conectar(red.salida, maq.terminales)
    if bloqueado:
        w0 = m.add(FuenteConstante("w0", 0.0))
        th0 = m.add(FuenteConstante("th0", 0.0))
        m.conectar(w0.salida, Puerto(maq, "ent", 3, 1))
        m.conectar(th0.salida, Puerto(maq, "ent", 4, 1))
    else:
        tl = m.add(FuenteConstante("tl", 0.0 if TL is None else TL))
        m.conectar(tl.salida, maq.T_L)
    return m, maq

def _rms_i(res, ultimos=2000):
    ia = np.asarray(res["I"])[:, 0][-ultimos:]
    return float(np.sqrt(np.mean(ia * ia)))

def _corre_3hp():
    m, maq = _armar(1e-4, 3, None)
    demux = m.add(Demultiplexor("demux_i", n_canales=3))
    m.conectar(maq.sensor3I(), demux.entrada)
    mux_ia = m.add(Multiplexor("mux_ia", n_canales=2))
    m.conectar(demux.salidas[0], mux_ia.entradas[0])
    m.conectar(demux.salidas[0], mux_ia.entradas[1])
    mult = m.add(Multiplicador("mult"))
    m.conectar(mux_ia.salida, mult.entrada)
    perd = m.add(Ganancia("perd", 3.0 * 0.435))
    m.conectar(mult.salida, perd.entrada)
    masa = m.add(MasaTermica("masa", C_th=5000.0, T_inicial=25.0,
                             T_amb=25.0, R_amb=0.08))
    m.conectar(perd.salida, masa.entradas[0])

    sc = m.add(Scope("scope", anchos=[1, 1, 3, 1],
                     guiones=["wm (rad/s)", "Te (N.m)",
                              "ia", "ib", "ic", "T motor (C)"]))
    m.conectar(maq.sensorVelocidad(), sc.canales[0])
    m.conectar(maq.sensorPar(), sc.canales[1])
    m.conectar(maq.sensor3I(), sc.canales[2])
    m.conectar(masa.salida, sc.canales[3])
    res = m.run(2.0, registrar=[sc, maq.sensorVelocidad(),
                                maq.sensorPar(), maq.sensor3I(), masa])
    Te = np.asarray(res["Te"]).ravel()
    ws = 2 * np.pi * 60.0 / 2.0
    k = int(round(0.05 / m.dt))
    return ws, Te, res, k, masa

m, maq = _armar(1e-4, 3, None, bloqueado=True)
res = m.run(0.6, registrar=[maq.sensorPar(), maq.sensor3I()])
k1, k2 = int(0.45 / m.dt), int(0.60 / m.dt)
Te_stall = np.asarray(res["Te"]).ravel()[k1:k2].mean()
ia = np.asarray(res["I"])[:, 0][k1:k2]
print("Motor de induccion 3 hp (Tabla 6.10-1 de Krause)")
print("=" * 62)
print(f"par de arranque (rotor bloqueado): {Te_stall:7.2f} N.m  "
      f"(regimen estacionario exacto: 52.9704; libro, Ej. 6B: 51.90 "
      f"N.m, de su formula aproximada 6.9-19)")
print(f"corriente de arranque            : {np.sqrt(np.mean(ia*ia)):7.2f} A rms "
      f"(modelo qd completo: 65.74; libro, 6B-2: 64.8 A rms = 64.8260, "
      f"sin rama magnetizante)")

ws, Te, res, k, masa = _corre_3hp()
print(f"velocidad a los 2 s       : {res['wm'][-1]:7.1f} rad/s "
      f"({100*res['wm'][-1]/ws:.1f} % del sincronismo {ws:.1f} rad/s)")
print(f"pico de corriente de fase : {np.abs(res['I']).max():7.2f} A "
      f"(eje de la Fig. 6.10-5: 83.2 A; depende de la fase de conexion)")
print(f"pico de par en el arranque: {np.abs(Te[:k]).max():7.2f} N.m "
      f"(eje de la Fig. 6.10-5: 118.7 N.m = 10*TB)")
print(f"monitoreo termico (aux)   : T motor a los 2 s = "
      f"{res['masa'][-1]:5.1f} C  (Q ~ 3*rs*ia^2, MasaTermica)")

print("Punto nominal a plena carga (TL = TB exacto)")
print("-" * 62)
for hp, Vll, rpm_lib, A_lib in (
        (3, 220.0, 1710.0, 5.8),
        (50, 460.0, 1705.0, 46.8)):
    TB = _tb(hp)
    m, maq = _armar(1e-4, hp, TB)
    ws = 2 * np.pi * 60.0 / 2.0
    maq.estados_iniciales[4] = ws
    res = m.run(0.8 if hp == 3 else 1.5,
                registrar=[maq.sensorVelocidad(), maq.sensor3I(),
                           maq.sensorPar()])
    wm = float(res["wm"][-1])
    slip = 100.0 * (ws - wm) / ws
    print(f"  {hp:3d} hp: TB = {TB:8.4f} N.m (libro: {TB:.1f}), "
          f"w = {wm:7.2f} rad/s ({wm*60/(2*np.pi):5.0f} rpm), "
          f"deslizamiento {slip:5.2f} %  (libro: {rpm_lib:.0f} rpm, "
          f"{100*(1800-rpm_lib)/1800:.2f} %),  "
          f"Ia rms = {_rms_i(res):5.2f} A  (IB del libro: {A_lib} A; "
          f"la IB es la corriente BASE, no la del circuito)")
