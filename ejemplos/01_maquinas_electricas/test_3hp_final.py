import numpy as np
import warnings
warnings.filterwarnings("ignore", message="numba cannot be imported")
from bloques_crysi import (
    Modelo, FuenteTrifasica, MaquinaInduccion, MasaTermica,
    FuenteEscalon,
)
WB = 2 * np.pi * 60.0
m = Modelo(dt=1e-7, metodo="rk4")

red = m.add(FuenteTrifasica("red", amplitud=220 * np.sqrt(2/3), frecuencia=60, fase=np.pi/3))
maq = m.add(MaquinaInduccion("mi", rs=0.435, rr=0.816, Lm=26.13/WB,
                             Lls=0.754/WB, Llr=0.754/WB, P=4, J=0.089, Bm=0.0))
m.conectar(red.salida, maq.terminales)
masa = m.add(MasaTermica("masa", 5000, maq.sensorPerdidasEstator(),
                         T_inicial=25.0, T_amb=25.0, R_amb=0.08))
m.scope("", maq.sensorVelocidad(), maq.sensorPar(),
        guiones=[r"$\omega_m\ [\mathrm{rad/s}]$", r"$T_e\ [\mathrm{N\cdot m}]$"],
        xy_mode=(0, 1), mostrar=True)
m.scope("Transitorio_Corrientes",
        maq.sensor3I()[0], maq.sensor3I()[1], maq.sensor3I()[2],
        guiones=[r"$i_{as}\ [\mathrm{A}]$", r"$i_{bs}\ [\mathrm{A}]$", r"$i_{cs}\ [\mathrm{A}]$"],
        mostrar=True)
m.scope("Corrientes_Rotor",
        maq.sensorCorrienteRotor(),
        guiones=[r"$i'_{ar}\ [\mathrm{A}]$", r"$i'_{br}\ [\mathrm{A}]$", r"$i'_{cr}\ [\mathrm{A}]$"],
        mostrar=True)
m.scope("Dinamica_Completa",
        maq.sensorVelocidad(), maq.sensorPar(),
        maq.sensorPerdidasEstator(), masa.salida,
        guiones=[r"$\omega_m\ [\mathrm{rad/s}]$",
                 r"$T_e\ [\mathrm{N\cdot m}]$",
                 r"$P_{\mathrm{Cu}}\ [\mathrm{W}]$",
                 r"$T\ [^\circ\mathrm{C}]$"],
        mostrar=True)
res = m.run(0.5, registrar=[maq.sensorVelocidad(), maq.sensorPar(),
                            maq.sensor3I(), masa.salida])
