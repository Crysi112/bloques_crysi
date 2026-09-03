import numpy as np

from bloques_crysi import (
    Modelo, FuenteTrifasica, FuenteConstante, MaquinaInduccion,
    FiltroPasoBajo, Scope,
)

WB = 2 * np.pi * 60.0

m = Modelo(dt=1e-4, metodo="rk4")

red = m.add(FuenteTrifasica("red", amplitud=220.0 * np.sqrt(2.0 / 3.0),
                            frecuencia=60.0))
tl = m.add(FuenteConstante("tl", 0.0))
maq = m.add(MaquinaInduccion(
    "mi", rs=0.435, rr=0.816, Lm=26.13 / WB, Lls=0.754 / WB,
    Llr=0.754 / WB, P=4, J=0.089, Bm=0.0))
m.conectar(red.salida, maq.terminales)
m.conectar(tl.salida, maq.T_L)

te_f = m.add(FiltroPasoBajo("te_f", fc=10.0, orden=1))

m.conectar(maq.sensorPar(), te_f.entrada)
sc1 = m.add(Scope("scope1", anchos=[1, 1, 1],
                  guiones=["wm (rad/s)", "Te (N.m)",
                           "Te filtrado (N.m)"]))

sc2 = m.add(Scope("scope2", anchos=[3],
                  guiones=["ia", "ib", "ic"]))
m.conectar(maq.sensorVelocidad(), sc1.canales[0])
m.conectar(maq.sensorPar(), sc1.canales[1])
m.conectar(te_f.salida, sc1.canales[2])
m.conectar(maq.sensor3I(), sc2.canales[0])

res = m.run(2.0, registrar=[sc1, sc2, maq.sensorVelocidad(),
                            maq.sensorPar(), maq.sensor3I(), te_f])
ws = 2 * np.pi * 60.0 / 2.0

print(f"velocidad final       : {res['wm'][-1]:7.1f} rad/s "
      f"({100*res['wm'][-1]/ws:.1f} % del sincronismo {ws:.1f} rad/s)")
print(f"pico de corriente     : {np.abs(res['I']).max():7.2f} A")
