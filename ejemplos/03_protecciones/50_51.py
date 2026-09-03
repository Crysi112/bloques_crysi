import math
import numpy as np
import matplotlib.pyplot as plt
from bloques_crysi.red import RedOpenDSS

red = RedOpenDSS(nombre="CasoIntegral_M5", v_slack_kv_ll=13.8, f_hz=60.0)
backend = red.compilar()

bus_fuente = backend.dss.text("? Vsource.source.bus1").strip()
if not bus_fuente: 
    bus_fuente = "sourcebus"

backend.dss.text("Edit Vsource.source basekv=13.8 pu=1.0 R1=0.0541 X1=0.5414 R0=0.1353 X0=1.3535")
backend.dss.text("New Linecode.FEEDER_M5 nphases=3 units=mi rmatrix=[0.225 | 0.075 0.225 | 0.075 0.075 0.225] xmatrix=[0.525 | 0.175 0.525 | 0.175 0.175 0.525]")
backend.dss.text(f"New Line.L_F_X bus1={bus_fuente} bus2=BusX length=1 units=mi linecode=FEEDER_M5 phases=3")
backend.dss.text("New Load.C_X bus1=BusX phases=3 kv=13.8 kw=2550 kvar=1580.34")

backend.dss.text("Set VoltageBases=[13.8]")
backend.dss.text("CalcVoltageBases")
backend.dss.text("Solve Mode=Snap")
backend.dss.text("Solve Mode=FaultStudy")

ruta = backend.dss.text("Export FaultStudy").strip()
fila = next((l for l in open(ruta) if "busx" in l.lower()), None)
if fila:
    datos = fila.split(',')
    i3ph, islg = float(datos[1]), float(datos[2])
else:
    i3ph, islg = 0.0, 0.0

I_load = 3e6 / (math.sqrt(3) * 13800)
ct = next(v for v in [50, 75, 100, 150, 200, 300, 400, 600, 800, 1200] if v >= I_load * 1.25)
tap = next(t for t in [4, 5, 6, 8, 10, 12] if t >= (I_load / (ct / 5)) * 1.5)
pickup_51 = tap * (ct / 5)
tms = 0.5 * ((i3ph / pickup_51)**0.02 - 1) / 0.14 if i3ph > pickup_51 else 0.0
pickup_50 = i3ph * 0.6

print(f"I_carga: {I_load:.1f}A | Falla 3φ: {i3ph:.0f}A | Falla SLG: {islg:.0f}A")
print(f"Relé 50 -> CT: {ct}:5 | Pickup Instantáneo: {pickup_50:.1f}A")
print(f"I_carga: {I_load:.1f}A | Falla 3φ: {i3ph:.0f}A | Falla SLG: {islg:.0f}A")
print(f"Relé 51 -> CT: {ct}:5 | TAP: {tap}A | Pickup: {pickup_51:.0f}A | TMS: {tms:.3f}")

i_51 = np.linspace(pickup_51, pickup_50, 500)
m = i_51 / pickup_51
t_51 = tms * (0.14 / (m**0.02 - 1))

plt.figure(figsize=(9, 6))
plt.loglog(i_51, t_51, label='Curva Compuesta 50/51 (Inversa + Truncamiento Instantáneo)', color='blue', linewidth=2)
plt.vlines(x=pickup_50, ymin=0.01, ymax=tms * (0.14 / ((pickup_50 / pickup_51)**0.02 - 1)), color='red', linestyle='-', linewidth=2, label=f'Pickup 50 ({pickup_50:.1f} A)')
plt.hlines(y=0.01, xmin=pickup_50, xmax=i3ph * 1.2, color='red', linestyle='-', linewidth=2)

plt.axvline(x=pickup_51, color='orange', linestyle=':', label=f'Pickup 51 ({pickup_51:.1f} A)')
plt.axvline(x=i3ph, color='purple', linestyle='-.', label=f'Falla 3φ ({i3ph:.0f} A)')

plt.xlabel('Corriente Primaria (A)', fontsize=11, fontweight='bold')
plt.ylabel('Tiempo de Disparo (s)', fontsize=11, fontweight='bold')
plt.title('Curva Característica TCC Coordinada (50/51)', fontsize=12, fontweight='bold')
plt.grid(True, which="both", ls="--", alpha=0.6)
plt.legend(loc='upper right', frameon=True)
plt.xlim(100, i3ph * 1.3)
plt.ylim(0.005, 100)
plt.tight_layout()
plt.show()
