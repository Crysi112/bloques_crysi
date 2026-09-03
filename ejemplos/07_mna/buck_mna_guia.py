import numpy as np
import matplotlib.pyplot as plt
from bloques_crysi import Modelo, PulsoRectangular, FuenteConstante
from bloques_crysi.mna import SubredMNA, Nodo, VSource, Switch, Diodo, Inductor, Capacitor, Resistor

dt_sim = 1e-6
m = Modelo(dt=dt_sim)

gnd = Nodo("0")
n_in = Nodo("in")
n_mid = Nodo("mid")
n_out = Nodo("out")

componentes = [
    VSource("Vin", n_in, gnd, idx_u=0),
    Switch("SW1", n_in, n_mid, idx_ctrl=1, Ron=1e-3, Roff=1e6),
    Diodo("D1", gnd, n_mid, Ron=1e-3, Roff=1e6, Vf=0.7),
    Inductor("L1", n_mid, n_out, L=100e-6),
    Capacitor("C1", n_out, gnd, C=10e-6),
    Resistor("Rload", n_out, gnd, R=5.0)
]

mediciones_v = [(n_out, gnd)]
mediciones_i = []

circuito = m.add(SubredMNA("Buck_MNA", nodos=[gnd, n_in, n_mid, n_out], componentes=componentes, dt=dt_sim, mediciones_v=mediciones_v, mediciones_i=mediciones_i, precomputar=True))

v_dc = m.add(FuenteConstante("Vdc_val", 24.0))
pwm = m.add(PulsoRectangular("PWM", amplitud=1.0, periodo=10e-6, duty=0.5))

m.conectar(v_dc.salida, circuito.entrada[0])
m.conectar(pwm.salida, circuito.entrada[1])

res = m.run(t_fin=0.002, registrar=[circuito.salida[0], pwm.salida])
vout = res["Buck_MNA"] if "Buck_MNA" in res else res["Buck_MNA[0]"]
print(f"Vout final: {float(vout[-1]):.3f} V (ideal 12V)")
print(f"Vout medio ultimos 0.5ms: {float(np.mean(vout[-500:])):.3f} V")
print(f"PWM samples: {res['PWM'][:5].tolist()}")

plt.figure(figsize=(10,4))
plt.plot(res.t*1e3, vout, label="Vout (MNA precomputado)")
plt.title("Buck MNA Guia - 24V 50% -> ~12V")
plt.xlabel("Tiempo [ms]"); plt.ylabel("Tension [V]")
plt.grid(); plt.legend()
plt.tight_layout()
plt.savefig("ejemplos/buck_mna_guia.png", dpi=150)
print("Guardado ejemplos/buck_mna_guia.png")
