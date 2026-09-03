from bloques_crysi import (
    Modelo, Bateria, FuenteConstante, FuenteEscalon, Suma, Puerto, Scope,
    Multiplexor, Multiplicador, Ganancia, MasaTermica,
)

DT = 1.0
T_END = 4 * 3600.0
I_DESC = 20.0
I_CARG = -20.0
T_CARGA = 2 * 3600.0
T_AMB = 25.0
C_TH = 1.2e4
R_AMB = 0.5

def _corre():
    m = Modelo(dt=DT)
    bat = m.add(Bateria("bat", tipo="plomo_acido"))
    i = m.add(Suma("i", signos=[1.0, -1.0]))
    desc = m.add(FuenteConstante("desc", I_DESC))
    carg = m.add(FuenteEscalon("carg", I_DESC - I_CARG, t_paso=T_CARGA))
    m.conectar(desc.salida, Puerto(i, "ent", 0, 1))
    m.conectar(carg.salida, Puerto(i, "ent", 1, 1))
    m.conectar(i.salida, bat.entrada)

    mux_i = m.add(Multiplexor("mux_i", n_canales=2))
    m.conectar(i.salida, mux_i.entradas[0])
    m.conectar(i.salida, mux_i.entradas[1])
    mult = m.add(Multiplicador("mult"))
    m.conectar(mux_i.salida, mult.entrada)
    perd = m.add(Ganancia("perd", bat.param[5]))
    m.conectar(mult.salida, perd.entrada)
    masa = m.add(MasaTermica("masa", C_th=C_TH, T_inicial=T_AMB,
                             T_amb=T_AMB, R_amb=R_AMB))
    m.conectar(perd.salida, masa.entradas[0])

    sc = m.add(Scope("scope", anchos=[1, 1, 1, 1],
                     guiones=["Vbat (V)", "SOC", "I (A)", "T celda (C)"]))
    m.conectar(bat.sensorVoltaje(), sc.canales[0])
    m.conectar(bat.sensorSOC(), sc.canales[1])
    m.conectar(i.salida, sc.canales[2])
    m.conectar(masa.salida, sc.canales[3])

    res = m.run(t_fin=T_END, registrar=[sc, masa])
    return res, bat

res, bat = _corre()
vbat = res["scope"][:, 0]
soc = res["scope"][:, 1]
tcel = res["masa"]
print("Bateria plomo-acido (Tremblay-Dessaint, tipo Simulink)")
print("=" * 62)
print(f"parametros        : E0={bat.param[0]:.4f} V, K={bat.param[1]:.5f},")
print(f"                    Q={bat.param[2]:g} Ah, A={bat.param[3]:.3f} V,")
print(f"                    B={bat.param[4]:.4f}, R={bat.param[5]:.5f} ohm,")
print(f"                    Vcap={bat.param[7]:.3f} V")
print(f"descarga 20 A     : V desde {vbat[0]:.2f} V, SOC {soc[0]*100:.1f} %")
print(f"  a las 2 h       : V={vbat[int(2*3600/DT)]:.2f} V, "
      f"SOC={soc[int(2*3600/DT)]*100:.1f} %")
print(f"carga -20 A       : V final={vbat[-1]:.2f} V, SOC final="
      f"{soc[-1]*100:.1f} %")
print(f"termica           : perdidas I^2*R = {I_DESC**2*bat.param[5]:.2f} W")
print(f"                    T celda {tcel[0]:.1f} -> {tcel[-1]:.1f} C  "
      f"(delta estacionario teorico {I_DESC**2*bat.param[5]*R_AMB:.2f} C)")

import matplotlib.pyplot as plt
res, bat = _corre()
vbat = res["scope"][:, 0]
soc = res["scope"][:, 1]
tcel = res["masa"]
t = res.t
plt.figure(figsize=(10, 4))
plt.plot(t / 3600, vbat, label="Vbat")
plt.plot(t / 3600, soc * 100, label="SOC (%)")
plt.plot(t / 3600, tcel, label="T celda (C)")
plt.xlabel("tiempo [h]")
plt.grid(True, alpha=0.3)
plt.legend()
plt.title("Bateria plomo-acido 12 V / 100 Ah: descarga, carga y termica")

bat.graficar_curvas(corriente=20.0)
plt.show()
