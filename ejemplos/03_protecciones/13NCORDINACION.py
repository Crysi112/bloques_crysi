"""Esquema de protecciones alimentador IEEE 13 nodos (ANSI/IEEE C37.2).

Datos + factory crear_proteccion + modulo bloques_crysi.tcc.
Sin funciones ni ciclos de calculo: todo lo reusable vive en la libreria.
"""

import cmath
import math
from pathlib import Path

import numpy as np
from bloques_crysi.red import RedOpenDSS
from bloques_crysi import (
    crear_proteccion, curva_51, curva_fusible, curva_itm, curva_dano,
    icc_faultstudy, verificar_cti, figura_tcc, tabla_ajustes, marcas_icc,
)
from bloques_crysi.tcc import K_FUSIBLE_K

BASE = Path(__file__).resolve().parent.parent.parent
IMG = BASE / "docs" / "img"
IMG.mkdir(parents=True, exist_ok=True)

NXFM1 = 4.16 / 0.48
Ibase_sub_mt = 5_000_000.0 / (math.sqrt(3) * 4160.0)
Ibase_xfm1_mt = 500_000.0 / (math.sqrt(3) * 4160.0)
Ibase_xfm1_bt = 500_000.0 / (math.sqrt(3) * 480.0)

IR_BT, ISD_BT, II_BT = 680.0, 2720.0, 8000.0
INRUSH_BT = 12.0 * Ibase_xfm1_bt
print(f"\nITM 634: Ii={II_BT:.0f} A vs inrush 12xIn={INRUSH_BT:.0f} A "
      f"(margen {II_BT / INRUSH_BT:.2f}x, debe ser >1.0)")

# ================================================================
# 1. RED IEEE 13 NODOS
# ================================================================
red = RedOpenDSS(nombre="IEEE_13_Completo", v_slack_kv_ll=115.0, f_hz=60.0)
backend = red.compilar()

bus_fuente = backend.dss.text("? Vsource.source.bus1").strip() or "sourcebus"
backend.script([
    "Edit Vsource.source basekv=115.0 pu=1.0 Isc3=10000 Isc1=10000",
    "New Transformer.Sub Phases=3 Windings=2 XHL=8",
    f"~ wdg=1 bus={bus_fuente} conn=Delta kv=115 kva=5000 %r=1",
    "~ wdg=2 bus=650 conn=Wye kv=4.16 kva=5000 %r=1",
    "New Transformer.XFM1 Phases=3 Windings=2 XHL=2 wdg=1 bus=633 conn=Wye kv=4.16 kva=500 wdg=2 bus=634 conn=Wye kv=0.48 kva=500",
    "New Linecode.601 nphases=3 units=mi rmatrix=[0.3465 | 0.1560 0.3375 | 0.1580 0.1535 0.3414] xmatrix=[1.0179 | 0.5017 1.0478 | 0.4236 0.3849 1.0348]",
    "New Linecode.602 nphases=3 units=mi rmatrix=[0.7526 | 0.1580 0.7475 | 0.1560 0.1535 0.7436] xmatrix=[1.1814 | 0.4236 1.1983 | 0.5017 0.3849 1.2112]",
    "New Linecode.603 nphases=2 units=mi rmatrix=[1.3294 | 0.2066 1.3238] xmatrix=[1.3471 | 0.4591 1.3569]",
    "New Linecode.604 nphases=2 units=mi rmatrix=[1.3294 | 0.2066 1.3238] xmatrix=[1.3471 | 0.4591 1.3569]",
    "New Linecode.605 nphases=1 units=mi rmatrix=[1.3294] xmatrix=[1.3471]",
    "New Linecode.606 nphases=3 units=mi rmatrix=[0.1600 | 0.0500 0.1600 | 0.0500 0.0500 0.1600] xmatrix=[0.1200 | 0.0400 0.1200 | 0.0400 0.0400 0.1200]",
    "New Linecode.607 nphases=1 units=mi rmatrix=[0.3500] xmatrix=[0.1000]",
    "New Line.L_650_632 bus1=650 bus2=632 length=2000 units=ft phases=3 linecode=601",
    "New Line.L_632_671 bus1=632 bus2=671 length=2000 units=ft phases=3 linecode=601",
    "New Line.L_632_633 bus1=632 bus2=633 length=500 units=ft phases=3 linecode=602",
    "New Line.L_632_645 bus1=632.2.3 bus2=645.2.3 length=500 units=ft phases=2 linecode=603",
    "New Line.L_645_646 bus1=645.2.3 bus2=646.2.3 length=300 units=ft phases=2 linecode=603",
    "New Line.L_671_684 bus1=671.1.3 bus2=684.1.3 length=300 units=ft phases=2 linecode=604",
    "New Line.L_684_611 bus1=684.3 bus2=611.3 length=300 units=ft phases=1 linecode=605",
    "New Line.L_684_652 bus1=684.1 bus2=652.1 length=800 units=ft phases=1 linecode=607",
    "New Line.L_671_692 bus1=671 bus2=692 length=100 units=ft phases=3 linecode=601",
    "New Line.L_692_675 bus1=692 bus2=675 length=500 units=ft phases=3 linecode=606",
    "New Load.Load_671 bus1=671 phases=3 kv=4.16 kw=1155 kvar=660 model=1 conn=Delta",
    "New Load.Load_634 bus1=634 phases=3 kv=0.48 kw=400 kvar=290 model=1 conn=Wye",
    "New Load.L_645B bus1=645.2 phases=1 kv=2.4 kw=170 kvar=125 model=1",
    "New Load.L_646BC bus1=646.2.3 phases=2 kv=4.16 kw=230 kvar=132 model=1 conn=Delta",
    "New Load.L_675A bus1=675.1 phases=1 kv=2.4 kw=485 kvar=190 model=1",
    "New Load.L_675B bus1=675.2 phases=1 kv=2.4 kw=68 kvar=60 model=1",
    "New Load.L_675C bus1=675.3 phases=1 kv=2.4 kw=290 kvar=212 model=1",
    "New Load.L_692 bus1=692 phases=3 kv=4.16 kw=170 kvar=151 model=1 conn=Delta",
    "New Capacitor.CAP675 bus1=675 phases=3 kv=4.16 kvar=600 conn=Wye",
    "New Load.L_611 bus1=611.3 phases=1 kv=2.4 kw=170 kvar=80 model=1",
    "New Capacitor.CAP611 bus1=611.3 phases=1 kv=2.4 kvar=100",
    "New Load.L_652 bus1=652.1 phases=1 kv=2.4 kw=128 kvar=86 model=1",
    "Set VoltageBases=[115, 4.16, 0.48]", "CalcVoltageBases", "Solve",
])

# ================================================================
# 2. CORTOCIRCUITO POR BARRA
# ================================================================
BUSES = ["650", "632", "671", "633", "634", "645", "646", "684", "611", "652", "692", "675"]
print("Icc 3F/1F por barra (A):")
icc3, icc1 = icc_faultstudy(backend, BUSES)

# ================================================================
# 3. RELES ANSI (FACTORY)
# ================================================================
r87_sub = crear_proteccion("87T", "87_Sub", I_min=0.2, slope1=0.25, slope2=0.60, I_break=2.0)
r51_AT = crear_proteccion("51P", "51_AT115", I_pickup=35.0, dial_tds=4.5, tipo_curva="IEEE_VERY_INVERSE")
r51NT = crear_proteccion("51NT", "51NT_Sub", I_pickup=140.0, dial_tds=5.0, tipo_curva="IEC_STANDARD_INVERSE")
r51_cab = crear_proteccion("51P", "51_Cab650", I_pickup=735.0, dial_tds=3.0, tipo_curva="IEEE_EXTREMELY_INVERSE")
r51N_cab = crear_proteccion("51N", "51N_Cab650", I_pickup=150.0, dial_tds=3.0, tipo_curva="IEEE_EXTREMELY_INVERSE")
r51_tr = crear_proteccion("51P", "51_Troncal", I_pickup=580.0, dial_tds=2.5, tipo_curva="IEEE_EXTREMELY_INVERSE")
r51N_tr = crear_proteccion("51N", "51N_Troncal", I_pickup=95.0, dial_tds=2.0, tipo_curva="IEC_STANDARD_INVERSE")
r51_671 = crear_proteccion("51P", "51_Acom671", I_pickup=230.0, dial_tds=1.5, tipo_curva="IEC_STANDARD_INVERSE")
r46_671 = crear_proteccion("46", "46_Acom671", I2_pickup=27.7, t_retardo=3.0)
r51B_645 = crear_proteccion("51P", "51B_Ramal645", I_pickup=180.0, dial_tds=2.0, tipo_curva="IEEE_MODERATELY_INVERSE")
r51C_645 = crear_proteccion("51P", "51C_Ramal645", I_pickup=85.0, dial_tds=2.0, tipo_curva="IEEE_MODERATELY_INVERSE")
r51N_645 = crear_proteccion("51N", "51N_Ramal645", I_pickup=145.0, dial_tds=2.0, tipo_curva="IEEE_MODERATELY_INVERSE")
r79_645 = crear_proteccion("79", "79_Ramal645", tiempos_muertos=(2.0, 15.0))
r79_tr = crear_proteccion("79", "79_Troncal")
r52_650 = crear_proteccion("52", "52_Cab650", t_apertura_mecanica=0.05)
r86_sub = crear_proteccion("86", "86_Sub")
r51_692 = crear_proteccion("51P", "51_Subt692", I_pickup=260.0, dial_tds=1.5, tipo_curva="IEEE_VERY_INVERSE")
r51N_692 = crear_proteccion("51N", "51N_Subt692", I_pickup=160.0, dial_tds=1.5, tipo_curva="IEEE_VERY_INVERSE")
r79_692_bloq = crear_proteccion("79", "79_Subt692_Bloqueado", tiempos_muertos=())
r51A_684 = crear_proteccion("51P", "51A_Cab684", I_pickup=80.0, dial_tds=1.5, tipo_curva="IEC_STANDARD_INVERSE")
r51C_684 = crear_proteccion("51P", "51C_Cab684", I_pickup=90.0, dial_tds=1.5, tipo_curva="IEC_STANDARD_INVERSE")
r51N_684 = crear_proteccion("51N", "51N_Cab684", I_pickup=90.0, dial_tds=1.5, tipo_curva="IEC_STANDARD_INVERSE")
r59N_675 = crear_proteccion("59N", "59N_Cap675", V_pickup=0.05, t_retardo=0.2)
r62_675 = crear_proteccion("62", "62_DescargaCap", t_retardo=300.0)
r50_LTC = crear_proteccion("50", "50_Lockout_Regulador", I_pickup=1400.0)

print("\nZona subterranea 692/675/652: 79 BLOQUEADO (sin recierre en XLPE)")
print("\nAjustes instanciados:")
tabla_ajustes((r87_sub, r51_AT, r51NT, r51_cab, r51N_cab, r51_tr, r51N_tr,
               r51_671, r46_671, r51B_645, r51C_645, r51N_645,
               r51_692, r51N_692, r51A_684, r51C_684, r51N_684,
               r59N_675, r62_675, r50_LTC))

ok_pas, id_pas, ir_pas = r87_sub.evaluar_disparo(2.51 + 0j, -2.51 + 0j)
ok_int, id_int, ir_int = r87_sub.evaluar_disparo(2.51 + 0j, 2.51 + 0j)
print(f"\n87T pasante: trip={ok_pas} (Idiff={id_pas:.2f}, Irest={ir_pas:.2f})")
print(f"87T interna: trip={ok_int} (Idiff={id_int:.2f}, Irest={ir_int:.2f})")

# ================================================================
# 4. SELECTIVIDAD (CTI)
# ================================================================
print("\nParejas de coordinacion (aguas abajo -> aguas arriba):")
verificar_cti([
    ("ITM-BT I/S @634", float(curva_itm(icc3["634"], ir=IR_BT, isd=ISD_BT, ii=II_BT)),
    
     "Fusible 100T @633", float(curva_fusible(icc3["634"] / NXFM1, 100.0)), 0.20),
    ("Fusible 100T @633", float(curva_fusible(icc3["633"], 100.0)),
     "51 Troncal @632", curva_51(r51_tr, icc3["633"]), 0.20),
    ("51 Acometida @671", curva_51(r51_671, icc3["671"]),
     "51 Troncal @632", curva_51(r51_tr, icc3["671"]), 0.30),
    ("51 Troncal @671", curva_51(r51_tr, icc3["671"]),
     "51 Cabecera @650", curva_51(r51_cab, icc3["671"]), 0.30),
    ("51B Ramal @645", curva_51(r51B_645, icc3["645"]),
     "51 Cabecera @650", curva_51(r51_cab, icc3["645"]), 0.30),
    ("51 Cabecera @650", curva_51(r51_cab, icc3["650"]),
     "51 AT115 (referido)", curva_51(r51_AT, icc3["650"] / (115.0 / 4.16)), 0.20),
    ("51 Subt692 @675", curva_51(r51_692, icc3["675"]),
     "51 Troncal @632", curva_51(r51_tr, icc3["675"]), 0.30),
    ("Fusible 100T @646", float(curva_fusible(icc3["646"], 100.0)),
     "51B Ramal @645", curva_51(r51B_645, icc3["646"]), 0.30),
    ("CLF 100A @652 (aprox. T)", float(curva_fusible(icc3["652"], 100.0)),
     "51A Cab684", curva_51(r51A_684, icc3["652"]), 0.20),
    ("Fusible 65K @611", float(curva_fusible(icc3["611"], 65.0, K=K_FUSIBLE_K, n=2.5, umbral_pu=1.35)),
     "51C Cab684", curva_51(r51C_684, icc3["611"]), 0.20),
    ("Fusible 125T lateral @684", float(curva_fusible(icc3["611"], 125.0)),
     "51C Cab684", curva_51(r51C_684, icc3["611"]), 0.20),
])

# ================================================================
# 5. COMPLEMENTARIAS 21/27/59/81/49/63
# ================================================================
ANG_LINEA = math.radians(75.0)
r21_Z1 = crear_proteccion("21_Z1", "21_L115_Z1", z_alcance_ohm=cmath.rect(0.64, ANG_LINEA), t_retardo=0.0)
r21_Z2 = crear_proteccion("21_Z2", "21_L115_Z2", z_alcance_ohm=cmath.rect(0.96, ANG_LINEA), t_retardo=0.4)
r59_N1 = crear_proteccion("59", "59_Sub_N1", V_pickup=1.10, t_retardo=10.0)
r59_N2 = crear_proteccion("59", "59_Sub_N2", V_pickup=1.20, t_retardo=0.16)
r27_N1 = crear_proteccion("27", "27_Mot671_N1", V_pickup=0.90, t_retardo=5.0)
r27_N2 = crear_proteccion("27", "27_Aisl_N2", V_pickup=0.80, t_retardo=1.0)
r81O = crear_proteccion("81O", "81O_Sub", f_pickup=60.5, t_retardo=2.0)
r81U_1 = crear_proteccion("81U", "81U_Sub_N1", f_pickup=59.5, t_retardo=10.0)
r81U_2 = crear_proteccion("81U", "81U_Ind671_N2", f_pickup=59.0, t_retardo=0.2)
r81U_3 = crear_proteccion("81U", "81U_Total_N3", f_pickup=58.5, t_retardo=0.1)
r49_sub = crear_proteccion("49", "49_TrafoSub", I_nominal=693.9)
r63_buch = crear_proteccion("63", "63_Buchholz", t_retardo=0.0)
r63_spr = crear_proteccion("63", "63_PresionSubita", t_retardo=0.0)

print("\nComplementarias (V/f/Z/termicas):")
tabla_ajustes((r21_Z1, r21_Z2, r59_N1, r59_N2, r27_N1, r27_N2,
               r81O, r81U_1, r81U_2, r81U_3, r49_sub, r63_buch, r63_spr))
print("  Placa termica: Top-Oil 90/105 C | Hot-Spot 110/120 C (IEEE C57.91)")

t1, _, _, ez1 = r21_Z1.paso(cmath.rect(0.40, ANG_LINEA), 1.0 + 0j, 0.02)
print(f"\n21 Z1: falla 50% -> en_zona={ez1:.0f} trip={t1:.0f} (instantaneo)")
r21_Z1b = crear_proteccion("21_Z1", "21_L115_Z1b", z_alcance_ohm=cmath.rect(0.64, ANG_LINEA), t_retardo=0.0)
t1b, _, _, ez1b = r21_Z1b.paso(cmath.rect(0.80, ANG_LINEA), 1.0 + 0j, 0.02)
t2, _, _, ez2 = r21_Z2.paso(cmath.rect(0.80, ANG_LINEA), 1.0 + 0j, 0.4)
print(f"21 Z1/Z2: falla 80% -> Z1 en_zona={ez1b:.0f} trip={t1b:.0f} | Z2 en_zona={ez2:.0f} trip={t2:.0f} (0.4 s)")
r21_Z1c = crear_proteccion("21_Z1", "21_L115_Z1c", z_alcance_ohm=cmath.rect(0.64, ANG_LINEA), t_retardo=0.0)
t1c, _, _, ez1c = r21_Z1c.paso(cmath.rect(2.00, ANG_LINEA), 1.0 + 0j, 0.02)
print(f"21 Z1: falla externa -> en_zona={ez1c:.0f} trip={t1c:.0f} (restringe)")
print(f"27 N2 (0.75 pu, 1 s): trip={r27_N2.paso(0.75, 1.0)[0]:.0f}")
print(f"59 N2 (1.25 pu, 0.16 s): trip={r59_N2.paso(1.25, 0.16)[0]:.0f}")
print(f"81U N2 (58.8 Hz, 0.2 s): trip={r81U_2.paso(58.8, 0.2)[0]:.0f}")
print(f"81O (60.6 Hz, 2 s): trip={r81O.paso(60.6, 2.0)[0]:.0f}")
print(f"63 Buchholz (senal=1): trip={r63_buch.paso(1.0, 0.01)[0]:.0f} (al 86 Lockout)")

# ================================================================
# 6. TCC INTERACTIVA UNIFICADA (TODO REFERIDO A 4.16 KV)
# ================================================================
i_plot = np.logspace(1, 4.5, 1000)
figura_tcc(
    curvas=[
        ("51 Cabecera 650 (735A)", i_plot, curva_51(r51_cab, i_plot), "black", "solid", 2.5, True),
        ("51 Troncal 632-671 (580A)", i_plot, curva_51(r51_tr, i_plot), "blue", "solid", 2.5, True),
        ("51NT Neutro Subestacion (140A)", i_plot, curva_51(r51NT, i_plot), "brown", "dash", 2.0, True),
        ("51 Acometida Ind. 671 (230A)", i_plot, curva_51(r51_671, i_plot), "red", "solid", 2.5, True),
        ("51B Ramal Bifasico 645 (180A)", i_plot, curva_51(r51B_645, i_plot), "green", "solid", 2.5, True),
        ("51 Subterraneo 692 (260A)", i_plot, curva_51(r51_692, i_plot), "purple", "solid", 2.5, "legendonly"),
        ("Fusible 100T (Trafo 633)", i_plot, curva_fusible(i_plot, 100.0), "orange", "dashdot", 2.5, True),
        ("Fusible 125T (Ramal 645 y Lateral 684)", i_plot, curva_fusible(i_plot, 125.0), "magenta", "dashdot", 2.5, True),
        ("Fusible 100T @646 (Carga B-C)", i_plot, curva_fusible(i_plot, 100.0), "cyan", "dashdot", 2.5, "legendonly"),
        ("Fusible 65K @611 (Cap. 100 kVAr)", i_plot, curva_fusible(i_plot, 65.0, K=K_FUSIBLE_K, n=2.5, umbral_pu=1.35), "olive", "dashdot", 2.5, "legendonly"),
        ("ITM 634 (BT ref. a 4.16kV)", i_plot, curva_itm(i_plot * NXFM1, ir=IR_BT, isd=ISD_BT, ii=II_BT), "darkred", "solid", 2.5, "legendonly"),
        ("Dano Trafo Principal 5MVA", i_plot, curva_dano(i_plot, Ibase_sub_mt), "darkgray", "solid", 2.5, "legendonly"),
        ("Dano Trafo XFM-1 500kVA", i_plot, curva_dano(i_plot, Ibase_xfm1_mt), "darkgoldenrod", "solid", 2.5, "legendonly"),
    ],
    marcas=marcas_icc(icc3, {"650": "black", "671": "blue", "645": "green", "633": "orange",
                             "634": "darkred", "675": "purple", "692": "teal",
                             "646": "cyan", "611": "olive", "652": "gray"}, {"634": NXFM1}),
    inrush=("Inrush XFM-1 (0.1s)", 833.0, 0.1),
    titulo="Coordinacion TCC Interactiva - IEEE 13 Nodos",
    subtitulo="Todas las curvas referidas a base Primaria 4.16 kV",
    archivo=IMG / "tcc_interactiva_ieee13.html",
)
