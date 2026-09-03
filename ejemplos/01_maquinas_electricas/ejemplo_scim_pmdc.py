import numpy as np

from bloques_crysi import (
    Modelo, FuenteTrifasica, FuenteConstante,
    MaquinaInduccion, MaquinaCorrienteContinua, Ganancia, Scope,
    Multiplexor, Multiplicador, MasaTermica,
)

R_LOAD = 1.8

m = Modelo(dt=1e-4)
red = m.add(FuenteTrifasica("red", amplitud=169.83, frecuencia=60.0))
scim = m.add(MaquinaInduccion(
    "scim", rs=12.45, rr=9.58, Lm=0.415, Lls=0.040, Llr=0.040,
    P=4, J=0.0022, mecanica_interna=True))
tl = m.add(FuenteConstante("tl", 0.0))
m.conectar(red.salida, scim.terminales)
m.conectar(tl.salida, scim.T_L)
res = m.run(t_fin=2.0, registrar=[scim.sensorVelocidad(),
                                  scim.sensorPar(), scim.sensor3I()])
wm = float(res["wm"][-1])
print("Motor SCIM solo (Lab-Volt 8221, 4 polos, 60 Hz)")
print("=" * 62)
print(f"velocidad sincrona : 188.50 rad/s (1800 RPM)")
print(f"velocidad final    : {wm:8.2f} rad/s  "
      f"({wm*60/(2*np.pi):6.0f} RPM, "
      f"slip {(188.5-wm)/188.5*100:.2f} %)")

m = Modelo(dt=1e-5)
red = m.add(FuenteTrifasica("red", amplitud=169.83, frecuencia=60.0))
scim = m.add(MaquinaInduccion(
    "scim", rs=12.45, rr=9.58, Lm=0.415, Lls=0.040, Llr=0.040,
    P=4, J=0.0022, mecanica_interna=False))
pmdc = m.add(MaquinaCorrienteContinua(
    "pmdc", r_a=0.295, L_a=436e-6, r_f=1.0, L_f=1e-3, L_AF=0.1098,
    J=0.001, mecanica_interna=False))
vf = m.add(FuenteConstante("vf", 1.0))
carga = m.add(Ganancia("carga", -R_LOAD))
m.conectar(red.salida, scim.terminales)
m.conectar(vf.salida, pmdc.campo)
m.conectar(pmdc.sensorCorriente(), carga.entrada)
m.conectar(carga.salida, pmdc.entrada)
eje = m.acoplar_maquinas(scim, pmdc, J_eq=0.003, Bm_eq=0.0)

mux_ia = m.add(Multiplexor("mux_ia", n_canales=2))
m.conectar(pmdc.sensorCorriente(), mux_ia.entradas[0])
m.conectar(pmdc.sensorCorriente(), mux_ia.entradas[1])
mult = m.add(Multiplicador("mult"))
m.conectar(mux_ia.salida, mult.entrada)
perd = m.add(Ganancia("perd", 0.295))
m.conectar(mult.salida, perd.entrada)
masa = m.add(MasaTermica("masa", C_th=1500.0, T_inicial=25.0,
                         T_amb=25.0, R_amb=0.35))
m.conectar(perd.salida, masa.entradas[0])

sc = m.add(Scope("banco", anchos=[1, 1, 1, 1],
                 guiones=["velocidad (rad/s)", "par (N.m)",
                          "ia generador (A)", "T devanado (C)"]))
m.conectar(scim.sensorVelocidad(), sc.canales[0])
m.conectar(scim.sensorPar(), sc.canales[1])
m.conectar(pmdc.sensorCorriente(), sc.canales[2])
m.conectar(masa.salida, sc.canales[3])

res = m.run(t_fin=2.0, registrar=[
    sc, scim.sensorVelocidad(), scim.sensorPar(), scim.sensor3I(),
    pmdc.sensorCorriente(), pmdc.sensorEa(), pmdc.sensorVoltajeTerminal(),
    masa,
])
t = res.t
k = t > 1.5
wm = np.mean(res["wm"][k])
ia_r = float(np.mean(res["ia"][k]))
tdev = float(res["masa"][-1])
print("Banco SCIM (motor) + PMDC (generador)")
print("=" * 62)
print(f"velocidad sincrona : 188.50 rad/s (1800 RPM, 60 Hz, 4 polos)")
print(f"velocidad regimen  : {wm:8.2f} rad/s  ({wm*60/(2*np.pi):6.0f} RPM, "
      f"slip {(188.5-wm)/188.5*100:.2f} %)")
print(f"par del motor      : {np.mean(res['Te'][k]):8.3f} N.m")
print(f"corriente de linea : {np.sqrt(np.mean(res['I'][k, 0]**2)):8.3f} A rms")
print(f"generador: ia      : {ia_r:8.3f} A  "
      f"Ea = {np.mean(res['Ea'][k]):6.2f} V  "
      f"V_t = {np.mean(res['V_t'][k]):6.2f} V")
print(f"perdidas de cobre  : {ia_r**2*0.295:8.3f} W  (Q = ia^2*r_a)")
print(f"T devanado final   : {tdev:8.1f} C  (MasaTermica, "
      f"delta estacionario {ia_r**2*0.295*0.35:.1f} C)")
print(f"placa del 8221     : 175 W, 208 V, 1.2 A, 1670 RPM a plena carga")
