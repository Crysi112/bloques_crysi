import numpy as np

from bloques_crysi import (
    Modelo, FuenteConstante, MaquinaDCImanesPermanentes, Scope,
    Multiplexor, Multiplicador, Ganancia, MasaTermica,
)

OZIN_NM = 0.00706155
MICRO_OZIN_S2 = 0.00706155

m = Modelo(dt=1e-5, metodo="rk4")
v = m.add(FuenteConstante("v", 12.0))
tl = m.add(FuenteConstante("tl", 0.0))
maq = m.add(MaquinaDCImanesPermanentes(
    "dc", r_a=3.8, L_a=1.14e-3, Kt=0.031, J=1.41e-5, Bm=2.82e-4))
m.conectar(v.salida, maq.entrada)
m.conectar(tl.salida, maq.T_L)
sc = m.add(Scope("scope", anchos=[1, 1],
                 guiones=["wm (rad/s)", "ia (A)"]))
m.conectar(maq.sensorVelocidad(), sc.canales[0])
m.conectar(maq.sensorCorriente(), sc.canales[1])
res = m.run(0.25, registrar=[sc, maq.sensorVelocidad(),
                             maq.sensorCorriente()])
wm = float(res["wm"][-1])
ia = float(res["ia"][-1])
wm_an = 0.031 * 12.0 / (0.031**2 + 2.82e-4 * 3.8)
ia_an = 2.82e-4 * wm_an / 0.031
print("Motor de CC fraccional de 12 V (Krause, cap. 10.11)")
print("=" * 62)
print(f"sin carga: wm = {wm:7.1f} rad/s ({wm*60/(2*np.pi):5.0f} rpm),  "
      f"ia = {ia:5.2f} A")
print(f"  analitico: {wm_an:7.1f} rad/s (1748 rpm),  ia = {ia_an:.2f} A")
print(f"  libro    : 183.0 rad/s (1748 rpm),  1.67 A (redondeado)")

kt = 2.0 * OZIN_NM
j = 150e-6 * MICRO_OZIN_S2
bm = 6.04e-6
print("Motor DC PM de juguete de 6 V (Krause, Ejemplo 10A)")
print("-" * 62)
print(f"kT = {kt:.7f} N*m/A  (2 oz*in/A; el libro redondea a 1.41e-2), "
      f"J = {j:.4e} kg*m2 (libro: 1.06e-6),  Bm = {bm:.3e} N*m*s/rad")
m = Modelo(dt=1e-5, metodo="rk4")
v = m.add(FuenteConstante("v", 6.0))
tl = m.add(FuenteConstante("tl", 0.0))
toy = m.add(MaquinaDCImanesPermanentes(
    "toy", r_a=7.0, L_a=0.12, Kt=kt, J=j, Bm=bm))
m.conectar(v.salida, toy.entrada)
m.conectar(tl.salida, toy.T_L)
mux_ia = m.add(Multiplexor("mux_ia", n_canales=2))
m.conectar(toy.sensorCorriente(), mux_ia.entradas[0])
m.conectar(toy.sensorCorriente(), mux_ia.entradas[1])
mult = m.add(Multiplicador("mult"))
m.conectar(mux_ia.salida, mult.entrada)
perd = m.add(Ganancia("perd", 7.0))
m.conectar(mult.salida, perd.entrada)
masa = m.add(MasaTermica("masa", C_th=60.0, T_inicial=25.0,
                         T_amb=25.0, R_amb=2.0))
m.conectar(perd.salida, masa.entradas[0])
sc = m.add(Scope("scope2", anchos=[1, 1, 1],
                 guiones=["wm (rad/s)", "ia (A)", "T devanado (C)"]))
m.conectar(toy.sensorVelocidad(), sc.canales[0])
m.conectar(toy.sensorCorriente(), sc.canales[1])
m.conectar(masa.salida, sc.canales[2])
res = m.run(0.25, registrar=[sc, toy.sensorVelocidad(),
                             toy.sensorCorriente(), masa])
wm = float(res["wm"][-1])
ia = float(res["ia"][-1])
wm_an = 6.0 / (kt + 7.0 * bm / kt)
ia_an = bm * wm_an / kt
print(f"sin carga: wm = {wm:7.1f} rad/s ({wm*60/(2*np.pi):5.0f} rpm), "
      f"ia = {ia:5.2f} A")
print(f"  analitico con Bm : {wm_an:7.1f} rad/s "
      f"({wm_an*60/(2*np.pi):5.0f} rpm),  ia = {ia_an:.2f} A")
print(f"  libro (10A-3)    : 351.1 rad/s = 3353 rpm,  ia = 0.15 A "
      f"(con kv redondeado a 1.41e-2)")
print(f"  limite ideal Bm=0: V/kv = {6.0/kt:.1f} rad/s "
      f"({6.0/kt*60/(2*np.pi):.0f} rpm)")
print(f"monitoreo termico (aux) : T devanado = "
      f"{res['masa'][-1]:.3f} C  (Q ~ ia^2*r_a ~ "
      f"{float(res['masa'][-1] - 25.0) * 60.0 / 0.25:.2f} W medios)")
