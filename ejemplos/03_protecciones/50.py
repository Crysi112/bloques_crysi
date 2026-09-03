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
pickup_50 = i3ph * 0.6

print(f"I_carga: {I_load:.1f}A | Falla 3φ: {i3ph:.0f}A | Falla SLG: {islg:.0f}A")
print(f"Relé 50 -> CT: {ct}:5 | Pickup Instantáneo: {pickup_50:.1f}A")

i_vals = np.linspace(100, i3ph * 1.2, 1000)
t_vals = np.where(i_vals >= pickup_50, 0.0, 10.0)

plt.figure(figsize=(9, 6))
plt.plot(i_vals, t_vals, label='50 - Instantáneo', color='red', linewidth=2)
plt.axvline(x=pickup_50, color='red', linestyle='--', label=f'Pickup 50 ({pickup_50:.1f} A)')
plt.xlabel('Corriente Primaria (A)', fontsize=11, fontweight='bold')
plt.ylabel('Tiempo de Disparo (s)', fontsize=11, fontweight='bold')
plt.title('Curva TCC - Función 50 (Instantáneo)', fontsize=12, fontweight='bold')
plt.grid(True, which="both", ls="--", alpha=0.6)
plt.legend(loc='upper right', frameon=True)
plt.xlim(100, i3ph * 1.3)
plt.ylim(-0.1, 2)
plt.tight_layout()
plt.show()
