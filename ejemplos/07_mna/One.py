from bloques_crysi import (Modelo, FuenteTrifasica, FuenteConstante,
                           MaquinaInduccion, MaquinaDCImanesPermanentes,
                           Ganancia, Scope, MedidorPotencia, Puerto)
import numpy as np

R_LOAD = 4.5517

m = Modelo(dt=dt, metodo="rk4")

red = m.add(FuenteTrifasica("red", amplitud=169.83, frecuencia=60.0))

scim = m.add(MaquinaInduccion(
    "scim", rs=12.45, rr=9.58, Lm=0.415, Lls=0.040, Llr=0.040,
    P=4, J=0.0022, mecanica_interna=False))

pmdc = m.add(MaquinaDCImanesPermanentes(
    "pmdc", r_a=0.295, L_a=436e-6, Kt=0.1098,
    J=0.001, mecanica_interna=False))

carga = m.add(Ganancia("carga", -R_LOAD))

m.conectar(red.salida, scim.terminales)
m.conectar(pmdc.sensorCorriente(), carga.entrada)
m.conectar(carga.salida, pmdc.entrada)
eje = m.acoplar_maquinas(scim, pmdc, J_eq=0.003, Bm_eq=0.0)

med_scim = m.add(MedidorPotencia("med_scim"))
med_pmdc = m.add(MedidorPotencia("med_pmdc", fases=1))
m.conectar(red.salida, med_scim.entrada)    

m.conectar(scim.sensor3I(), med_scim.corrientes)
m.conectar(scim.sensorPar(), Puerto(med_scim, "ent", 6, 1))
m.conectar(scim.sensorVelocidad(), Puerto(med_scim, "ent", 7, 1))
m.conectar(pmdc.sensorVoltajeTerminal(), med_pmdc.entrada)
m.conectar(pmdc.sensorCorriente(), med_pmdc.corrientes)
m.conectar(pmdc.sensorPar(), Puerto(med_pmdc, "ent", 2, 1))
m.conectar(pmdc.sensorVelocidad(), Puerto(med_pmdc, "ent", 3, 1))

sc_banco = m.add(Scope("Banco", anchos=[1, 1, 1, 3, 2],
                       guiones=["wm (rad/s)", "Te (N.m)", "ia (A)",
                                "P SCIM (W)", "Q SCIM (var)", "P_m SCIM (W)",
                                "P_e PMDC (W)", "P_m PMDC (W)"]))

m.conectar(scim.sensorVelocidad(), sc_banco.canales[0])
m.conectar(scim.sensorPar(), sc_banco.canales[1])
m.conectar(pmdc.sensorCorriente(), sc_banco.canales[2])
m.conectar(med_scim.salida, sc_banco.canales[3])
m.conectar(med_pmdc.salida, sc_banco.canales[4])

res = m.run(t_fin=0.5, registrar=[
    sc_banco, scim.sensorVelocidad(), scim.sensorPar(), scim.sensor3I(),
    pmdc.sensorCorriente(), pmdc.sensorEa(), pmdc.sensorVoltajeTerminal(),
    med_scim.salida, med_pmdc.salida,
])

k = res.t > 0.4
wm = np.mean(res["wm"][k])
print("Banco SCIM (motor) + PMDC (generador), R_LOAD =", R_LOAD, "ohm")
print("velocidad sincrona : 188.50 rad/s (1800 RPM, 60 Hz, 4 polos)")
print(f"velocidad regimen  : {wm:8.2f} rad/s  ({wm*60/(2*np.pi):6.0f} RPM, "
      f"slip {(188.5-wm)/188.5*100:.2f} %)")
print(f"par del motor      : {np.mean(res['Te'][k]):8.3f} N.m")
print(f"corriente de linea : {np.sqrt(np.mean(res['I'][k, 0]**2)):8.3f} A rms")
print(f"generador: ia      : {np.mean(res['ia'][k]):8.3f} A  "
      f"Ea = {np.mean(res['Ea'][k]):6.2f} V  "
      f"V_t = {np.mean(res['V_t'][k]):6.2f} V")
print(f"SCIM (motor)      : P_e = {np.mean(res['med_scim'][k, 0]):7.1f} W  "
      f"P_m = {np.mean(res['med_scim'][k, 1]):7.1f} W")
print(f"PMDC (generador)  : P_e = {np.mean(res['med_pmdc'][k, 0]):7.1f} W  "
      f"P_m = {np.mean(res['med_pmdc'][k, 1]):7.1f} W  (P_e < 0: entrega a la carga)")
print("placa del 8221     : 175 W, 208 V, 1.2 A, 1670 RPM a plena carga")
