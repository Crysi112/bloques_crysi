import numpy as np
import matplotlib.pyplot as plt
from bloques_crysi import Modelo, FuenteConstante, ConvertidorBuck
from bloques_crysi.bloques import PulsoRectangular
from bloques_crysi.mna import SubredMNA, Nodo, Resistor, Capacitor, Inductor, VSource, Switch, Diodo

Vdc = 12.0
L = 1e-3
C = 100e-6
R = 10.0
fsw = 5000.0
duty = 0.5
dt = 1e-6
t_fin = 0.01

m1 = Modelo(dt=dt)
with m1:
    vdc = m1.add(FuenteConstante("vdc", Vdc))
    pwm = m1.add(PulsoRectangular("pwm", amplitud=1.0, periodo=1/fsw, duty=duty))
    buck = m1.add(ConvertidorBuck("buck", L=L, C=C, R=R))
    m1.conectar(vdc.salida, buck.entrada)
    m1.conectar(pwm.salida, buck.d)

res1 = m1.run(t_fin, registrar=["buck"])
vout_mono = res1["buck"][:,0]
t = res1.t

gnd = Nodo("GND")
n1 = Nodo("N1")
n2 = Nodo("N2")
n3 = Nodo("N3")

vs = VSource("Vs", n1, gnd, idx_u=0)
sw = Switch("S1", n1, n2, idx_ctrl=1, Ron=1e-3, Roff=1e6)
l = Inductor("L", n2, n3, L)
cap = Capacitor("C", n3, gnd, C)
r = Resistor("R", n3, gnd, R)
diodo = Diodo("D1", gnd, n2, Ron=1e-3, Roff=1e6, Vf=0.0)

red = SubredMNA("buck_mna", nodos=[gnd, n1, n2, n3],
                componentes=[vs, sw, l, cap, r, diodo],
                dt=dt,
                mediciones_v=[(n3, gnd)],
                mediciones_i=[vs])

m2 = Modelo(dt=dt)
with m2:
    vdc2 = m2.add(FuenteConstante("vdc", Vdc))
    pwm2 = m2.add(PulsoRectangular("pwm", amplitud=1.0, periodo=1/fsw, duty=duty))
    bloque = m2.add(red)
    m2.conectar(vdc2.salida, bloque.entrada[0:1])
    m2.conectar(pwm2.salida, bloque.entrada[1:2])

res2 = m2.run(t_fin, registrar=["buck_mna"])
vout_mna = np.asarray(res2["buck_mna"])[:,0]
i_vs_mna = np.asarray(res2["buck_mna"])[:,1]

print(f"Vout monolitico final: {vout_mono[-1]:.3f} V (esperado ~{Vdc*duty:.1f} V)")
print(f"Vout MNA final:        {vout_mna[-1]:.3f} V")
print(f"I_Vs MNA pico:         {np.max(np.abs(i_vs_mna)):.3f} A")
print(f"Error relativo Vout:   {abs(vout_mna[-1]-vout_mono[-1])/vout_mono[-1]*100:.2f}%")

plt.figure(figsize=(10,6))
plt.subplot(2,1,1)
plt.plot(t*1e3, vout_mono, label="Monolitico (promediado)")
plt.plot(t*1e3, vout_mna, '--', label="SubredMNA (conmutado)")
plt.axhline(Vdc*duty, color='k', linestyle=':', label=f"Vdc*duty={Vdc*duty:.1f}V")
plt.ylabel("Vout [V]")
plt.legend()
plt.title("Buck 12V -> ~6V : Monolitico vs SubredMNA (MNA+Dommel+LU)")

plt.subplot(2,1,2)
plt.plot(t*1e3, i_vs_mna, label="I_Vs (MNA)")
plt.ylabel("I [A]")
plt.xlabel("t [ms]")
plt.legend()
plt.tight_layout()
plt.savefig("demo_buck_mna.png", dpi=150)
print("Grafico guardado en demo_buck_mna.png")
