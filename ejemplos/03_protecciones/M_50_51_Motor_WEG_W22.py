import math
import numpy as np
import matplotlib.pyplot as plt
from bloques_crysi.red import RedOpenDSS

v_ll = 0.220
red = RedOpenDSS(nombre="Motor_Protection", v_slack_kv_ll=v_ll, f_hz=60.0)
backend = red.compilar()

bus_fuente = backend.dss.text("? Vsource.source.bus1").strip()
if not bus_fuente:
    bus_fuente = "sourcebus"

R1_src, X1_src = 0.0008543, 0.0042714
backend.dss.text(
    f"Edit Vsource.source basekv={v_ll} pu=1.0 "
    f"R1={R1_src} X1={X1_src} R0={R1_src} X0={X1_src}"
)


backend.dss.text(
    "New Linecode.CABLE_MTR nphases=3 units=km "
    "rmatrix=[2.56 | 0.0 2.56 | 0.0 0.0 2.56] "
    "xmatrix=[0.171 | 0.0 0.171 | 0.0 0.0 0.171]"
)

L_km = 0.0196
backend.dss.text(f"New Line.L_Motor bus1={bus_fuente} bus2=BusMotor length={L_km} phases=3 linecode=CABLE_MTR")

hp = 7.5
eff_motor = 0.917
pf_motor = 0.81

kw_motor = (hp * 0.746) / eff_motor
kvar_motor = kw_motor * math.tan(math.acos(pf_motor))
backend.dss.text(f"New Load.MotorLoad bus1=BusMotor phases=3 kv={v_ll} kw={kw_motor:.3f} kvar={kvar_motor:.3f} model=1")

backend.dss.text(f"Set VoltageBases=[{v_ll}]")
backend.dss.text("CalcVoltageBases")
backend.dss.text("Solve")

V_LN = (v_ll * 1000) / math.sqrt(3)

Z1_src = complex(R1_src, X1_src)
Z0_src = complex(R1_src, X1_src)  

r_line = 2.56 * L_km
x_line = 0.171 * L_km

rm_line = 0.0 * L_km
xm_line = 0.0 * L_km

Z1_line = complex(r_line, x_line)
Z0_line = complex(r_line + 2 * rm_line, x_line + 2 * xm_line)

Z1_total = Z1_src + Z1_line
Z0_total = Z0_src + Z0_line

Icc3_motor = V_LN / abs(Z1_total)
Icc1_motor = (3 * V_LN) / abs(2 * Z1_total + Z0_total)

fla = 19.46
lra = 138.19
fuse_std = 35.0

pickup_51 = fla * 1.25
pickup_50 = lra * 1.7
tms_51 = 0.3

K = 250.0 

i_51 = np.linspace(pickup_51, pickup_50, 500)
m = i_51 / pickup_51
t_51 = tms_51 * (0.14 / (m**0.02 - 1))

i_fuse = np.linspace(fuse_std, Icc3_motor * 1.5, 500)
t_fuse = K * (fuse_std / i_fuse)**2

S_trafo = 500_000.0
I_fla_trafo = S_trafo / (math.sqrt(3) * (v_ll * 1000))
Z_pu_trafo = 0.045
I_max_trafo = I_fla_trafo / Z_pu_trafo
K_trafo = 1250.0  

i_pu_trafo_valido = np.linspace(2.0, min(40.0, I_max_trafo / I_fla_trafo), 500)
i_trafo_valido = i_pu_trafo_valido * I_fla_trafo
t_trafo_valido = K_trafo / i_pu_trafo_valido**2

I_pu_check = Icc3_motor / I_fla_trafo
t_check = K_trafo / I_pu_check**2
extrapola = I_pu_check < 2.0  
area_mm2 = 8.37   
K_cable = 115.0   
constante_dano_cable = (K_cable * area_mm2)**2
i_cable = np.linspace(pickup_51, Icc3_motor * 2.0, 500)
t_cable = constante_dano_cable / (i_cable**2)

print(f"FLA: {fla} A | LRA: {lra} A | Fusible: {fuse_std} A")
print(f"Cable: 8 AWG (Seccion={area_mm2} mm2, Constante K={K_cable})")
print(f"Pickup 51: {pickup_51:.2f} A | Pickup 50: {pickup_50:.1f} A")
print(f"Cortocircuito Trifasico (Motor): {Icc3_motor:.1f} A")
print(f"Cortocircuito Monofasico (Motor): {Icc1_motor:.1f} A")


plt.figure(figsize=(10, 7))

plt.loglog(i_51, t_51, label='Curva Compuesta 50/51', color='blue', linewidth=2)
plt.vlines(x=pickup_50, ymin=0.01,
           ymax=tms_51 * (0.14 / ((pickup_50 / pickup_51)**0.02 - 1)),
           color='red', linestyle='-', linewidth=2, label=f'Pickup 50 ({pickup_50:.1f} A)')
plt.hlines(y=0.01, xmin=pickup_50, xmax=Icc3_motor * 1.5, color='red', linestyle='-', linewidth=2)
plt.loglog(i_fuse, t_fuse, label=f'Fusible ({fuse_std:.0f} A)', color='green', linestyle='-.', linewidth=2)

plt.loglog(i_trafo_valido, t_trafo_valido,
           label='Curva Dano Trafo (Cat. I, C57.109, 2x-40x)', color='darkred', linestyle='-', linewidth=2)

plt.loglog(i_cable, t_cable, label='Dano Conductor (8 AWG, Cu/PVC)', color='magenta', linestyle=':', linewidth=2.5)

plt.axvline(x=pickup_51, color='orange', linestyle=':', label=f'Pickup 51 ({pickup_51:.1f} A)')
plt.axvline(x=lra, color='purple', linestyle='--', label=f'LRA Arranque ({lra:.1f} A)')
plt.axvline(x=Icc1_motor, color='brown', linestyle='-.', linewidth=1.5, label=f'Icc1 Monofasico ({Icc1_motor:.1f} A)')
plt.axvline(x=Icc3_motor, color='black', linestyle='-.', linewidth=1.5, label=f'Icc3 Trifasico ({Icc3_motor:.1f} A)')

plt.xlabel('Corriente Primaria (A)', fontsize=11, fontweight='bold')
plt.ylabel('Tiempo de Disparo (s)', fontsize=11, fontweight='bold')
plt.title('Coordinacion TCC: Motor 7.5 HP Trifasico (220V) - WEG 00718ET3E213T-W22', fontsize=12, fontweight='bold')
plt.grid(True, which="both", ls="--", alpha=0.6)
plt.legend(loc='upper right', frameon=True, fontsize=9)

plt.xlim(10, 50000)
plt.ylim(0.005, 10000)
plt.tight_layout()
plt.show()