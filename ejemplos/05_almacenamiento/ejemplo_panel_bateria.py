def _ruta(rel):
    for p in (os.path.join(os.getcwd(), 'ejemplos', rel),
              os.path.join(os.getcwd(), rel)):
        if os.path.exists(p):
            return p
    return os.path.join(os.getcwd(), rel)

import os

from bloques_crysi import (
    Modelo, PanelSolar, Bateria, FuenteCSV, FuenteConstante, FuenteEscalon,
    Suma, Puerto, Scope,
)

DT = 1.0
T_END = 2 * 3600.0
T_CEL = 35.0
V_DESC = 36.0
T_CARGA = 3600.0
CSV_G = _ruta("data/irradiancia.csv")

def _corre():
    m = Modelo(dt=DT)
    pan = m.add(PanelSolar("pan"))
    g = m.add(FuenteCSV("g", CSV_G))
    t = m.add(FuenteConstante("t", T_CEL))
    v = m.add(FuenteEscalon("v", V_DESC, t_paso=T_CARGA,
                            valor_inicial=30.0))
    m.conectar(g.salida, Puerto(pan, "ent", 0, 1))
    m.conectar(t.salida, Puerto(pan, "ent", 1, 1))
    m.conectar(v.salida, Puerto(pan, "ent", 2, 1))

    bat = m.add(Bateria("bat", tipo="plomo_acido", Vnom=24.0,
                        SOCinit=60.0, eta_c=0.95))
    carg = m.add(Suma("carg", signos=[1.0, -1.0]))
    cons = m.add(FuenteConstante("cons", 4.0))
    m.conectar(pan.salida, Puerto(carg, "ent", 0, 1))
    m.conectar(cons.salida, Puerto(carg, "ent", 1, 1))
    m.conectar(carg.salida, bat.entrada)

    sc = m.add(Scope("scope", anchos=[1, 1, 1, 1, 1],
                     guiones=["I panel (A)", "Vbat (V)", "SOC",
                              "I bat (A)", "G (W/m2)"]))
    m.conectar(pan.salida, sc.canales[0])
    m.conectar(bat.sensorVoltaje(), sc.canales[1])
    m.conectar(bat.sensorSOC(), sc.canales[2])
    m.conectar(carg.salida, sc.canales[3])
    m.conectar(g.salida, sc.canales[4])

    res = m.run(t_fin=T_END, registrar=[sc])
    return res

res = _corre()
ip = res["scope"][:, 0]
vbat = res["scope"][:, 1]
soc = res["scope"][:, 2]
print("Panel solar con bateria (plomo-acido 24 V), perfil G desde CSV")
print("=" * 62)
print(f"perfil G          : 800 -> nube 250 (t=0.75-1.25 h) -> 800 W/m2")
print(f"panel a G=800     : I(V=36) = {ip[0]:.2f} A  (~{ip[0]*36.0:.0f} W)")
print(f"en la nube (1 h)  : G = 400 W/m2 -> I panel = "
      f"{ip[int(3600/DT)]:.2f} A  (punto minimo G=250 a t=0.75 h)")
print(f"consumo           : 4 A -> I bat inicial = {ip[0]-4.0:+.2f} A")
print(f"SOC 0 h           : {soc[0]*100:.1f} %")
print(f"SOC 1 h (nube)    : {soc[int(T_CARGA/DT)]*100:.1f} %")
print(f"SOC 2 h (V=36)    : {soc[-1]*100:.1f} %")
print(f"Vbat final        : {vbat[-1]:.2f} V")

import matplotlib.pyplot as plt
res = _corre()
ip = res["scope"][:, 0]
vbat = res["scope"][:, 1]
soc = res["scope"][:, 2]
t = res.t
plt.figure(figsize=(10, 4))
plt.plot(t / 3600, ip, label="I panel (A)")
plt.plot(t / 3600, vbat, label="Vbat (V)")
plt.plot(t / 3600, soc * 100, label="SOC (%)")
plt.xlabel("tiempo [h]")
plt.grid(True, alpha=0.3)
plt.legend()
plt.title("Panel solar 260 Wp con bateria: nube en el perfil de irradiancia")

pan = PanelSolar("pan")
pan.graficar_curvas(Gs=(200, 500, 800, 1000), T=T_CEL)
plt.show()
