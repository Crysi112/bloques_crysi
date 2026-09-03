import csv
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

backend.dss.text("Solve mode=FaultStudy")
ruta_faultstudy = backend.dss.text("Export Faultstudy").strip()
icc3_opendss = None
icc1_opendss = None
try:
    with open(ruta_faultstudy, newline="") as f:
        lector = csv.DictReader(f)
        lector.fieldnames = [c.strip() for c in lector.fieldnames]
        print("Columnas exportadas por OpenDSS:", lector.fieldnames)
        for fila in lector:
            fila_limpia = {k.strip(): v.strip() for k, v in fila.items()}
            if fila_limpia["Bus"].upper() == "BUSMOTOR":
                icc3_opendss = float(fila_limpia["3-Phase"])
                icc1_opendss = float(fila_limpia["1-Phase"])
                break
except (FileNotFoundError, KeyError) as e:
    print(f"No se pudo leer/parsear el export de FaultStudy ({e}). "
          f"Revisa el nombre de columnas impreso arriba.")

if icc3_opendss is not None:
    diff_3f = (icc3_opendss - Icc3_motor) / Icc3_motor * 100
    print(f"Icc3 manual:   {Icc3_motor:.1f} A")
    print(f"Icc3 OpenDSS:  {icc3_opendss:.1f} A  (dif. {diff_3f:+.2f}%)")
if icc1_opendss is not None:
    diff_1f = (icc1_opendss - Icc1_motor) / Icc1_motor * 100
    print(f"Icc1 manual:   {Icc1_motor:.1f} A")
    print(f"Icc1 OpenDSS:  {icc1_opendss:.1f} A  (dif. {diff_1f:+.2f}%)")


fla = 19.46
lra = 138.19
fuse_std = 35.0

pickup_51 = fla * 1.25
pickup_50 = lra * 1.7

A, B, p = 28.2, 0.1217, 2.0
td_51 = 3.0

K = 250.0  
i_51 = np.linspace(pickup_51 * 1.001, pickup_50, 500)
m = i_51 / pickup_51
t_51 = td_51 * (A / (m**p - 1) + B)

m_pickup50 = pickup_50 / pickup_51
t_en_pickup50 = td_51 * (A / (m_pickup50**p - 1) + B)

i_fuse = np.linspace(fuse_std, Icc3_motor * 1.5, 500)
t_fuse = K * (fuse_std / i_fuse) ** 2

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
constante_daño_cable = (K_cable * area_mm2)**2
i_cable = np.linspace(pickup_51, Icc3_motor * 2.0, 500)
t_cable = constante_daño_cable / (i_cable**2)

ct_prim, ct_sec = 75.0, 5.0
rtc = ct_prim / ct_sec
ct_clase_alf, ct_sn_va = 20.0, 15.0

rho_cu75, area_ct_mm2, long_ct_m = 0.0214, 5.26, 20.0
R_cable_ct = 2.0 * long_ct_m * rho_cu75 / area_ct_mm2
S_cable_ct = ct_sec**2 * R_cable_ct
S_rele_ct = 0.10            
R_contactos = 0.05          
S_contactos = ct_sec**2 * R_contactos
S_total_ct = S_cable_ct + S_rele_ct + S_contactos

n_req = Icc3_motor / ct_prim
Rct_est = 0.05              
S_int_est = ct_sec**2 * Rct_est
alf_real = ct_clase_alf * (ct_sn_va + S_int_est) / (S_total_ct + S_int_est)
Rct_max = ((ct_clase_alf * ct_sn_va - n_req * S_total_ct) / (n_req - ct_clase_alf)) / ct_sec**2
Vs_req = (Icc3_motor / rtc) * (R_cable_ct + 0.004 + R_contactos + Rct_est)

fla_sec = fla / rtc
pk51_sec = pickup_51 / rtc
pk50_sec = pickup_50 / rtc
pk50_xfla = pickup_50 / fla
lra_sec = lra / rtc
icc_sec = Icc3_motor / rtc

i_alf_real = alf_real * ct_prim          
i_ad_sel = 110.0 * rtc                   
i_ct_therm = np.linspace(max(i_alf_real, pickup_50), Icc3_motor * 2.0, 300)
t_ct_therm = (4500.0**2 * 1.0) / i_ct_therm**2   

checks_ct = [
    ("51P 1.622A en [0.50, 10.00]A sec.", 0.50 <= pk51_sec <= 10.00),
    ("50P 12.07xFLA en [0.10, 20.00]xFLA", 0.10 <= pk50_xfla <= 20.00),
    ("TD=3 en [0.50, 15.00] (US)", 0.50 <= td_51 <= 15.00),
    ("LRA sec 9.21A < 15A continuo", lra_sec < 15.0),
    ("Icc sec 164.1A < 500A/1s", icc_sec < 500.0),
    ("Ith CT 4500A > Icc 2461.7A", 60.0 * ct_prim >= Icc3_motor),
    ("ALF real 48.7 >= requerido 32.82", alf_real >= n_req),
    ("Burden 5.42VA < 15VA nominal", S_total_ct < ct_sn_va),
]

print(f"FLA: {fla} A | LRA: {lra} A | Fusible: {fuse_std} A")
print(f"Cable: 8 AWG (Seccion={area_mm2} mm2, Constante K={K_cable})")
print(f"Pickup 51: {pickup_51:.2f} A | Pickup 50: {pickup_50:.1f} A")
print(f"Cortocircuito Trifasico (Motor, manual): {Icc3_motor:.1f} A")
print(f"Cortocircuito Monofasico (Motor, manual): {Icc1_motor:.1f} A")
print(f"CT MBS SASK 31.6 75/5A 5P20 15VA (RTC={rtc:.0f})")
print(f"Burden: cable {S_cable_ct:.2f} + rele {S_rele_ct:.2f} + contactos {S_contactos:.2f} = {S_total_ct:.2f} VA ({S_total_ct/ct_sn_va*100:.1f}% de 15VA)")
print(f"ALF requerido: {n_req:.2f} | ALF real (Rct={Rct_est} ohm): {alf_real:.1f} | Rct max admisible: {Rct_max:.3f} ohm")
print(f"Vs requerida: {Vs_req:.1f} V | Secundario: FLA {fla_sec:.3f}A, 51pk {pk51_sec:.3f}A, 50pk {pk50_sec:.2f}A ({pk50_xfla:.2f}xFLA), LRA {lra_sec:.2f}A, Icc {icc_sec:.1f}A")
for nombre, ok in checks_ct:
    print(f"[{'OK' if ok else 'FALLA'}] {nombre}")
if not all(ok for _, ok in checks_ct):
    raise SystemExit("Verificacion CT/SEL-710 fallida: revisar ajustes.")

plt.figure(figsize=(10, 7))

plt.loglog(i_51, t_51, label='Curva 51 (EI, IEEE C37.112, TD=3 provisional)', color='blue', linewidth=2)
plt.vlines(x=pickup_50, ymin=0.01, ymax=t_en_pickup50,
           color='red', linestyle='-', linewidth=2, label=f'Pickup 50 ({pickup_50:.1f} A)')
plt.hlines(y=0.01, xmin=pickup_50, xmax=Icc3_motor * 1.5, color='red', linestyle='-', linewidth=2)
plt.loglog(i_fuse, t_fuse, label=f'Fusible ({fuse_std:.0f} A)', color='green', linestyle='-.', linewidth=2)

plt.loglog(i_trafo_valido, t_trafo_valido,
           label='Curva Daño Trafo (Cat. I, C57.109, 2x-40x)', color='darkred', linestyle='-', linewidth=2)

plt.loglog(i_cable, t_cable, label='Daño Conductor (8 AWG, Cu/PVC)', color='magenta', linestyle=':', linewidth=2.5)

plt.loglog(i_ct_therm, t_ct_therm, label='Daño Térmico CT (Ith 4500A/1s)', color='teal', linestyle=':', linewidth=1.5)
plt.axvline(x=i_alf_real, color='teal', linestyle='--', linewidth=1.5, label=f'Límite lineal CT (ALF real {alf_real:.1f}, {i_alf_real:.0f} A)')
plt.axvline(x=i_ad_sel, color='gray', linestyle=':', linewidth=1.5, label=f'Saturación A/D SEL-710 ({i_ad_sel:.0f} A)')

plt.axvline(x=pickup_51, color='orange', linestyle=':', label=f'Pickup 51 ({pickup_51:.1f} A)')
plt.axvline(x=lra, color='purple', linestyle='--', label=f'LRA Arranque ({lra:.1f} A)')
plt.axvline(x=Icc1_motor, color='brown', linestyle='-.', linewidth=1.5, label=f'Icc1 Monofasico ({Icc1_motor:.1f} A)')
plt.axvline(x=Icc3_motor, color='black', linestyle='-.', linewidth=1.5, label=f'Icc3 Trifasico ({Icc3_motor:.1f} A)')

if icc3_opendss is not None:
    plt.axvline(x=icc3_opendss, color='black', linestyle=':', linewidth=1.2,
                label=f'Icc3 OpenDSS ({icc3_opendss:.1f} A)')

plt.xlabel('Corriente Primaria (A)', fontsize=11, fontweight='bold')
plt.ylabel('Tiempo de Disparo (s)', fontsize=11, fontweight='bold')
plt.title('Coordinacion TCC: Motor 7.5 HP Trifasico (220V) - WEG 00718ET3E213T-W22', fontsize=12, fontweight='bold')
plt.grid(True, which="both", ls="--", alpha=0.6)
plt.legend(loc='upper right', frameon=True, fontsize=9)

plt.xlim(10, 50000)
plt.ylim(0.005, 10000)
plt.tight_layout()
plt.show()