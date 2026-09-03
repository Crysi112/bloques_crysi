import numpy as np
from bloques_crysi import Modelo
from bloques_crysi.bloques import FuenteConstante, PulsoRectangular
from bloques_crysi.mna import SubredMNA, Nodo, Resistor, Capacitor, Inductor, VSource, Switch, Diodo

Vdc=12.0; L=1e-3; C=100e-6; R=10.0; fsw=5000; duty=0.5; dt=1e-6; t_fin=0.005

gnd=Nodo("GND"); n1=Nodo("N1"); n2=Nodo("N2"); n3=Nodo("N3")
vs=VSource("Vs", n1, gnd, idx_u=0)
sw=Switch("S1", n1, n2, idx_ctrl=1, Ron=1e-3, Roff=1e6)
l=Inductor("L", n2, n3, L)
cap=Capacitor("C", n3, gnd, C)
r=Resistor("R", n3, gnd, R)
diodo=Diodo("D1", gnd, n2, Ron=1e-3, Roff=1e6, Vf=0.0)

red=SubredMNA("buck", nodos=[gnd,n1,n2,n3], componentes=[vs,sw,l,cap,r,diodo], dt=dt,
              mediciones_v=[(n1,gnd),(n2,gnd),(n3,gnd)], mediciones_i=[vs])
print("n_x", red.n_x, "n_u", red.n_u, "param len", len(red.param))
print("param header", red.param[:15])

m=Modelo(dt=dt)
with m:
    vdc=m.add(FuenteConstante("vdc", Vdc))
    pwm=m.add(PulsoRectangular("pwm", amplitud=1.0, periodo=1/fsw, duty=duty))
    b=m.add(red)
    m.conectar(vdc.salida, b.entrada[0:1])
    m.conectar(pwm.salida, b.entrada[1:2])

res=m.run(t_fin, registrar=["buck","pwm","vdc"])
print(res["buck"].shape)
v1=res["buck"][:,0]; v2=res["buck"][:,1]; v3=res["buck"][:,2]; i_vs=res["buck"][:,3] if res["buck"].shape[1]>3 else None
pwm_sig=np.asarray(res["pwm"]).ravel() if res["pwm"].ndim==1 else res["pwm"][:,0]
print("v1 mean", np.mean(v1[-1000:]), "v2 mean", np.mean(v2[-1000:]), "v3 mean", np.mean(v3[-1000:]))
print("v1 max", np.max(v1), "v2 max", np.max(v2), "v3 max", np.max(v3))
for k in range(0, 10):
    print(f"t={res.t[k]*1e6:.0f}us pwm={pwm_sig[k]:.0f} v1={v1[k]:.2f} v2={v2[k]:.2f} v3={v3[k]:.4f}")
idx_on = np.where(pwm_sig>0.5)[0]
idx_off = np.where(pwm_sig<=0.5)[0]
print("when pwm=1, v2 close to v1?", np.mean(v2[idx_on[-100:]]), "v1", np.mean(v1[idx_on[-100:]]))
print("when pwm=0, v2 should be ~0 (diode clamp)", np.mean(v2[idx_off[-100:]]))
