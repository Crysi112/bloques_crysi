"""Esquema de protecciones alimentador IEEE 13 nodos (ANSI/IEEE C37.2).

Secciones 1-4 y 6-7 segun memoria de calculo: 87T/50/51/51NT en
subestacion, cabecera 650, troncal 632-671, ramal 632-645 (51B/51C/51N+79),
XFM-1 (fusible 100T + ITM LSI+G en BT) y acometida industrial 671 (51+46).
"""

import csv
import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from bloques_crysi.red import RedOpenDSS
from bloques_crysi import crear_proteccion

INF = float("inf")
BASE = Path(__file__).resolve().parent.parent.parent
IMG = BASE / "docs" / "img"
IMG.mkdir(parents=True, exist_ok=True)

# ================================================================
# 1. RED IEEE 13 NODOS (MODELO PARCIAL: 650/632/671/633/634/645)
# ================================================================
v_ll_kv = 4.16
red = RedOpenDSS(nombre="IEEE_13_Completo", v_slack_kv_ll=115.0, f_hz=60.0)
backend = red.compilar()

bus_fuente = backend.dss.text("? Vsource.source.bus1").strip() or "sourcebus"

backend.dss.text("Edit Vsource.source basekv=115.0 pu=1.0 Isc3=10000 Isc1=10000")

backend.dss.text("New Transformer.Sub Phases=3 Windings=2 XHL=8")
backend.dss.text(f"~ wdg=1 bus={bus_fuente} conn=Delta kv=115 kva=5000 %r=1")
backend.dss.text("~ wdg=2 bus=650 conn=Wye kv=4.16 kva=5000 %r=1")

backend.dss.text("New Linecode.601 nphases=3 units=mi rmatrix=[0.3465 | 0.1560 0.3375 | 0.1580 0.1535 0.3414] xmatrix=[1.0179 | 0.5017 1.0478 | 0.4236 0.3849 1.0348]")
backend.dss.text("New Linecode.602 nphases=3 units=mi rmatrix=[0.7526 | 0.1580 0.7475 | 0.1560 0.1535 0.7436] xmatrix=[1.1814 | 0.4236 1.1983 | 0.5017 0.3849 1.2112]")
backend.dss.text("New Linecode.603 nphases=2 units=mi rmatrix=[1.3294 | 0.2066 1.3238] xmatrix=[1.3471 | 0.4591 1.3569]")
backend.dss.text("New Linecode.604 nphases=2 units=mi rmatrix=[1.3294 | 0.2066 1.3238] xmatrix=[1.3471 | 0.4591 1.3569]")
backend.dss.text("New Linecode.605 nphases=1 units=mi rmatrix=[1.3294] xmatrix=[1.3471]")
backend.dss.text("New Linecode.606 nphases=3 units=mi rmatrix=[0.1600 | 0.0500 0.1600 | 0.0500 0.0500 0.1600] xmatrix=[0.1200 | 0.0400 0.1200 | 0.0400 0.0400 0.1200]")
backend.dss.text("New Linecode.607 nphases=1 units=mi rmatrix=[0.3500] xmatrix=[0.1000]")

backend.dss.text("New Line.L_650_632 bus1=650 bus2=632 length=2000 units=ft phases=3 linecode=601")
backend.dss.text("New Line.L_632_671 bus1=632 bus2=671 length=2000 units=ft phases=3 linecode=601")
backend.dss.text("New Line.L_632_633 bus1=632 bus2=633 length=500 units=ft phases=3 linecode=602")
backend.dss.text("New Line.L_632_645 bus1=632.2.3 bus2=645.2.3 length=500 units=ft phases=2 linecode=603")
backend.dss.text("New Line.L_645_646 bus1=645.2.3 bus2=646.2.3 length=300 units=ft phases=2 linecode=603")
backend.dss.text("New Line.L_671_684 bus1=671.1.3 bus2=684.1.3 length=300 units=ft phases=2 linecode=604")
backend.dss.text("New Line.L_684_611 bus1=684.3 bus2=611.3 length=300 units=ft phases=1 linecode=605")
backend.dss.text("New Line.L_684_652 bus1=684.1 bus2=652.1 length=800 units=ft phases=1 linecode=607")
backend.dss.text("New Line.L_671_692 bus1=671 bus2=692 length=100 units=ft phases=3 linecode=601")
backend.dss.text("New Line.L_692_675 bus1=692 bus2=675 length=500 units=ft phases=3 linecode=606")

backend.dss.text("New Transformer.XFM1 Phases=3 Windings=2 XHL=2 wdg=1 bus=633 conn=Wye kv=4.16 kva=500 wdg=2 bus=634 conn=Wye kv=0.48 kva=500")
backend.dss.text("New Load.Load_671 bus1=671 phases=3 kv=4.16 kw=1155 kvar=660 model=1 conn=Delta")
backend.dss.text("New Load.Load_634 bus1=634 phases=3 kv=0.48 kw=400 kvar=290 model=1 conn=Wye")
backend.dss.text("New Load.L_645B bus1=645.2 phases=1 kv=2.4 kw=170 kvar=125 model=1")
backend.dss.text("New Load.L_646BC bus1=646.2.3 phases=2 kv=4.16 kw=230 kvar=132 model=1 conn=Delta")
backend.dss.text("New Load.L_675A bus1=675.1 phases=1 kv=2.4 kw=485 kvar=190 model=1")
backend.dss.text("New Load.L_675B bus1=675.2 phases=1 kv=2.4 kw=68 kvar=60 model=1")
backend.dss.text("New Load.L_675C bus1=675.3 phases=1 kv=2.4 kw=290 kvar=212 model=1")
backend.dss.text("New Load.L_692 bus1=692 phases=3 kv=4.16 kw=170 kvar=151 model=1 conn=Delta")
backend.dss.text("New Capacitor.CAP675 bus1=675 phases=3 kv=4.16 kvar=600 conn=Wye")
backend.dss.text("New Load.L_611 bus1=611.3 phases=1 kv=2.4 kw=170 kvar=80 model=1")
backend.dss.text("New Capacitor.CAP611 bus1=611.3 phases=1 kv=2.4 kvar=100")
backend.dss.text("New Load.L_652 bus1=652.1 phases=1 kv=2.4 kw=128 kvar=86 model=1")

backend.dss.text("Set VoltageBases=[115, 4.16, 0.48]")
backend.dss.text("CalcVoltageBases")
backend.dss.text("Solve")

# ================================================================
# 2. CORTOCIRCUITO (FAULTSTUDY) EN BARRAS MODELADAS
# ================================================================
backend.dss.text("Solve mode=FaultStudy")
ruta_fs = backend.dss.text("Export Faultstudy").strip()
icc3, icc1 = {}, {}
try:
    with open(ruta_fs, newline="") as f:
        lector = csv.DictReader(f)
        lector.fieldnames = [c.strip() for c in lector.fieldnames]
        for fila in lector:
            fila = {k.strip(): v.strip() for k, v in fila.items()}
            try:
                icc3[fila["Bus"].upper()] = float(fila["3-Phase"])
                icc1[fila["Bus"].upper()] = float(fila["1-Phase"])
            except (KeyError, ValueError):
                pass
except OSError as e:
    print(f"Sin export FaultStudy ({e}); se continua sin Icc de OpenDSS.")

BUSES = ["650", "632", "671", "633", "634", "645", "646", "684", "611", "652", "692", "675"]
print("Icc 3F/1F por barra (A):")
for b in BUSES:
    print(f"  {b}: 3F={icc3.get(b, float('nan')):.1f}  1F={icc1.get(b, float('nan')):.1f}")

# ================================================================
# 3. RELES ANSI (FACTORY) - AJUSTES DE LA MEMORIA DE CALCULO
# ================================================================
r87_sub = crear_proteccion("87T", "87_Sub", I_min=0.2, slope1=0.25, slope2=0.60, I_break=2.0)
r51_AT = crear_proteccion("51P", "51_AT115", I_pickup=35.0, dial_tds=4.5,
                          tipo_curva="IEEE_VERY_INVERSE")
r51NT = crear_proteccion("51NT", "51NT_Sub", I_pickup=140.0, dial_tds=5.0,
                         tipo_curva="IEC_STANDARD_INVERSE")
r51_cab = crear_proteccion("51P", "51_Cab650", I_pickup=735.0, dial_tds=3.0,
                           tipo_curva="IEEE_EXTREMELY_INVERSE")
r51N_cab = crear_proteccion("51N", "51N_Cab650", I_pickup=150.0, dial_tds=3.0,
                            tipo_curva="IEEE_EXTREMELY_INVERSE")
r51_tr = crear_proteccion("51P", "51_Troncal", I_pickup=580.0, dial_tds=2.5,
                          tipo_curva="IEEE_EXTREMELY_INVERSE")
r51N_tr = crear_proteccion("51N", "51N_Troncal", I_pickup=95.0, dial_tds=2.0,
                           tipo_curva="IEC_STANDARD_INVERSE")
r51_671 = crear_proteccion("51P", "51_Acom671", I_pickup=230.0, dial_tds=1.5,
                           tipo_curva="IEC_STANDARD_INVERSE")
r46_671 = crear_proteccion("46", "46_Acom671", I2_pickup=27.7, t_retardo=3.0)
r51B_645 = crear_proteccion("51P", "51B_Ramal645", I_pickup=180.0, dial_tds=2.0,
                            tipo_curva="IEEE_MODERATELY_INVERSE")
r51C_645 = crear_proteccion("51P", "51C_Ramal645", I_pickup=85.0, dial_tds=2.0,
                            tipo_curva="IEEE_MODERATELY_INVERSE")
r51N_645 = crear_proteccion("51N", "51N_Ramal645", I_pickup=145.0, dial_tds=2.0,
                            tipo_curva="IEEE_MODERATELY_INVERSE")
r79_645 = crear_proteccion("79", "79_Ramal645", tiempos_muertos=(2.0, 15.0))
r79_tr = crear_proteccion("79", "79_Troncal")
r52_650 = crear_proteccion("52", "52_Cab650", t_apertura_mecanica=0.05)
r86_sub = crear_proteccion("86", "86_Sub")

# Transicion aereo-subterranea 692 -> 675 (desensibilizado por desbalance)
r51_692 = crear_proteccion("51P", "51_Subt692", I_pickup=260.0, dial_tds=1.5,
                           tipo_curva="IEEE_VERY_INVERSE")
r51N_692 = crear_proteccion("51N", "51N_Subt692", I_pickup=160.0, dial_tds=1.5,
                            tipo_curva="IEEE_VERY_INVERSE")
# Sin recierre en zona subterranea (XLPE no admite 79)
r79_692_bloq = crear_proteccion("79", "79_Subt692_Bloqueado", tiempos_muertos=())

# Cabecera derivacion bifasica 671 -> 684 (Sec. 6.1)
r51A_684 = crear_proteccion("51P", "51A_Cab684", I_pickup=80.0, dial_tds=1.5,
                            tipo_curva="IEC_STANDARD_INVERSE")
r51C_684 = crear_proteccion("51P", "51C_Cab684", I_pickup=90.0, dial_tds=1.5,
                            tipo_curva="IEC_STANDARD_INVERSE")
r51N_684 = crear_proteccion("51N", "51N_Cab684", I_pickup=90.0, dial_tds=1.5,
                            tipo_curva="IEC_STANDARD_INVERSE")

# Bancos de capacitores (IEEE 1036 / NEC 460)
r59N_675 = crear_proteccion("59N", "59N_Cap675", V_pickup=0.05, t_retardo=0.2)
r62_675 = crear_proteccion("62", "62_DescargaCap", t_retardo=300.0)

# Bloqueo LTC del regulador 650-632 durante falla (IEEE C57.15)
r50_LTC = crear_proteccion("50", "50_Lockout_Regulador", I_pickup=1400.0)

print("\nZona subterranea 692/675/652: 79 BLOQUEADO (sin recierre en XLPE)")

print("\nAjustes instanciados:")
for r in (r87_sub, r51_AT, r51NT, r51_cab, r51N_cab, r51_tr, r51N_tr,
          r51_671, r46_671, r51B_645, r51C_645, r51N_645,
          r51_692, r51N_692, r51A_684, r51C_684, r51N_684,
          r59N_675, r62_675, r50_LTC):
    print(f"  {r.codigo_ansi:6s} {r.nombre:20s} param={r.param}")

# Demo 87T: condicion pasante (no opera) vs falla interna (opera)
ok_pas, id_pas, ir_pas = r87_sub.evaluar_disparo(2.51 + 0j, -2.51 + 0j)
ok_int, id_int, ir_int = r87_sub.evaluar_disparo(2.51 + 0j, 2.51 + 0j)
print(f"\n87T pasante: trip={ok_pas} (Idiff={id_pas:.2f}, Irest={ir_pas:.2f})")
print(f"87T interna: trip={ok_int} (Idiff={id_int:.2f}, Irest={ir_int:.2f})")

# ================================================================
# 4. CURVAS AUXILIARES (FUSIBLES T, ITM BT, DANO TRAFOS C57.109)
# ================================================================
K_T = 0.35 * (833.0 / 100.0) ** 2  # fusible 100T: 0.35 s a 833 A (inrush XFM-1)
K_K = 0.15 * (500.0 / 65.0) ** 2   # fusible 65K: aprox. 0.054 s a 500 A

def t_fuseT(i_rat, i):
    i = np.asarray(i, dtype=float)
    return np.where(i > i_rat, K_T * (i_rat / i) ** 2, INF)

def t_fuseK(i_rat, i):
    i = np.asarray(i, dtype=float)
    return np.where(i > (i_rat * 1.35), K_K * (i_rat / i) ** 2.5, INF)

IR_BT, ISD_BT, II_BT = 680.0, 2720.0, 6800.0  # ITM 800 A marco, LSI+G

def t_itm(i):
    i = np.asarray(i, dtype=float)
    t = np.full_like(i, INF)
    m_l = (i >= IR_BT) & (i < ISD_BT)
    m_s = (i >= ISD_BT) & (i < II_BT)
    m_i = i >= II_BT
    t[m_l] = 10.0 * (6.0 * IR_BT / i[m_l]) ** 2
    t[m_s] = 0.2
    t[m_i] = 0.02
    return t

NXFM1 = 4.16 / 0.48  # relacion XFM-1 para referir corrientes

def t_dano_trafo(i, i_base, K=1250.0):
    i = np.asarray(i, dtype=float)
    return np.where(i > 0, K / (i / i_base) ** 2, INF)

Ibase_sub_mt = 5_000_000.0 / (math.sqrt(3) * 4160.0)   # 693.9 A
Ibase_xfm1_bt = 500_000.0 / (math.sqrt(3) * 480.0)     # 601.4 A
Ibase_xfm1_mt = 500_000.0 / (math.sqrt(3) * 4160.0)    # 69.4 A

def t51(rele, i):
    i = np.atleast_1d(np.asarray(i, dtype=float))
    return np.array([rele.calcular_tiempo_trip(float(v)) for v in i])

def t51s(rele, i):
    return float(t51(rele, i)[0])

# ================================================================
# 5. VERIFICACION DE SELECTIVIDAD (CTI)
# ================================================================
def fmt(t):
    return f"{t:.3f} s" if np.isfinite(t) else "no arranca"

print("\nParejas de coordinacion (aguas abajo -> aguas arriba):")
parejas = [
    ("ITM-BT I/S @634", lambda: float(t_itm(icc3["634"])),
     "Fusible 100T @633", lambda: float(t_fuseT(100.0, icc3["634"] / NXFM1)), 0.20),
    ("Fusible 100T @633", lambda: float(t_fuseT(100.0, icc3["633"])),
     "51 Troncal @632", lambda: t51s(r51_tr, icc3["633"]), 0.20),
    ("51 Acometida @671", lambda: t51s(r51_671, icc3["671"]),
     "51 Troncal @632", lambda: t51s(r51_tr, icc3["671"]), 0.30),
    ("51 Troncal @671", lambda: t51s(r51_tr, icc3["671"]),
     "51 Cabecera @650", lambda: t51s(r51_cab, icc3["671"]), 0.30),
    ("51B Ramal @645", lambda: t51s(r51B_645, icc3["645"]),
     "51 Cabecera @650", lambda: t51s(r51_cab, icc3["645"]), 0.30),
    ("51 Cabecera @650", lambda: t51s(r51_cab, icc3["650"]),
     "51 AT115 (referido)", lambda: t51s(r51_AT, icc3["650"] / (115.0 / 4.16)), 0.20),
    ("51 Subt692 @675", lambda: t51s(r51_692, icc3["675"]),
     "51 Troncal @632", lambda: t51s(r51_tr, icc3["675"]), 0.30),
    ("Fusible 100T @646", lambda: float(t_fuseT(100.0, icc3["646"])),
     "51B Ramal @645", lambda: t51s(r51B_645, icc3["646"]), 0.30),
    ("CLF 100A @652 (aprox. T)", lambda: float(t_fuseT(100.0, icc3["652"])),
     "51A Cab684", lambda: t51s(r51A_684, icc3["652"]), 0.20),
    ("Fusible 65K @611", lambda: float(t_fuseK(65.0, icc3["611"])),
     "51C Cab684", lambda: t51s(r51C_684, icc3["611"]), 0.20),
]
for dw, f_dw, up, f_up, cti_min in parejas:
    try:
        t_dw, t_up = f_dw(), f_up()
        cti = t_up - t_dw
        ok = bool(np.isfinite(cti) and cti >= cti_min)
        print(f"  [{'OK' if ok else 'REVISAR'}] {dw} ({fmt(t_dw)}) < {up} ({fmt(t_up)}) "
              f"CTI={fmt(cti)} (min {cti_min:.2f} s)")
    except KeyError as e:
        print(f"  [SIN DATO] {dw} vs {up}: falta Icc en barra {e}")

# ================================================================
# 6. TCC MT 4.16 KV
# ================================================================
i_mt = np.logspace(math.log10(30.0), math.log10(30000.0), 500)

plt.figure(figsize=(10, 7))
plt.loglog(i_mt, t51(r51_cab, i_mt), label="51 Cabecera 650 (735 A, EI, TDS 3.0)", color="black", linewidth=2)
plt.loglog(i_mt, t51(r51_tr, i_mt), label="51 Troncal 632-671 (580 A, EI, TDS 2.5)", color="blue", linewidth=2)
plt.loglog(i_mt, t51(r51_671, i_mt), label="51 Acometida 671 (230 A, SI, TDS 1.5)", color="red", linewidth=2)
plt.loglog(i_mt, t51(r51B_645, i_mt), label="51B Ramal 645 (180 A, MI, TDS 2.0)", color="green", linewidth=2)
plt.loglog(i_mt, t51(r51NT, i_mt), label="51NT Subestacion (140 A, SI, TDS 5.0)", color="brown", linestyle="--", linewidth=1.5)
plt.loglog(i_mt, t_fuseT(100.0, i_mt), label="Fusible 100T @633", color="orange", linestyle="-.", linewidth=2)
plt.loglog(i_mt, t_fuseT(125.0, i_mt), label="Fusible 125T @645", color="magenta", linestyle="-.", linewidth=1.5)
plt.loglog(i_mt, t_fuseT(100.0, i_mt), label="Fusible 100T @646 (Carga B-C)", color="cyan", linestyle="-.", linewidth=1.5)
plt.loglog(i_mt, t_fuseK(65.0, i_mt), label="Fusible 65K @611 (Cap. 100 kVAr)", color="olive", linestyle="-.", linewidth=1.5)
plt.loglog(i_mt, t_dano_trafo(i_mt, Ibase_sub_mt), label="Dano trafo 5 MVA (C57.109 Cat. II)", color="darkred", linewidth=1.5)
plt.plot(833.0, 0.1, "ro", label="Inrush XFM-1 (833 A, 0.1 s)")

for b, c in (("650", "black"), ("671", "blue"), ("645", "green"), ("633", "orange"),
             ("675", "purple"), ("646", "cyan"), ("611", "olive"), ("652", "gray")):
    if b in icc3:
        plt.axvline(x=icc3[b], color=c, linestyle=":", linewidth=1.2, label=f"Icc3@{b} ({icc3[b]:.0f} A)")

plt.xlabel("Corriente MT (A)", fontsize=11, fontweight="bold")
plt.ylabel("Tiempo (s)", fontsize=11, fontweight="bold")
plt.title("IEEE 13 nodos: coordinacion MT 4.16 kV", fontsize=12, fontweight="bold")
plt.grid(True, which="both", ls="--", alpha=0.6)
plt.legend(loc="upper right", frameon=True, fontsize=8)
plt.xlim(30, 30000)
plt.ylim(0.01, 2000)
plt.tight_layout()
plt.savefig(IMG / "tcc_13nodos_mt.png", dpi=150)
plt.close()
print(f"Figura MT: {IMG / 'tcc_13nodos_mt.png'}")

# ================================================================
# 7. TCC BT 480 V (XFM-1 + ITM)
# ================================================================
i_bt = np.logspace(math.log10(300.0), math.log10(60000.0), 500)

plt.figure(figsize=(10, 7))
plt.loglog(i_bt, t_itm(i_bt), label="ITM 634 (L 680 A / S 2720 A 0.2 s / I 6800 A)", color="red", linewidth=2)
plt.loglog(i_bt, t_fuseT(100.0, i_bt / NXFM1), label="Fusible 100T referido a BT", color="orange", linestyle="-.", linewidth=2)
plt.loglog(i_bt, t_dano_trafo(i_bt, Ibase_xfm1_bt), label="Dano XFM-1 500 kVA (C57.109 Cat. I)", color="darkred", linewidth=1.5)
if "634" in icc3:
    plt.axvline(x=icc3["634"], color="black", linestyle=":", linewidth=1.5, label=f"Icc3@634 ({icc3['634']:.0f} A)")
plt.axhline(y=0.2, color="gray", linestyle=":", linewidth=1.0, label="Banda S/G 0.2 s")

plt.xlabel("Corriente BT (A)", fontsize=11, fontweight="bold")
plt.ylabel("Tiempo (s)", fontsize=11, fontweight="bold")
plt.title("IEEE 13 nodos: coordinacion BT 480 V (nodo 634)", fontsize=12, fontweight="bold")
plt.grid(True, which="both", ls="--", alpha=0.6)
plt.legend(loc="upper right", frameon=True, fontsize=9)
plt.xlim(300, 60000)
plt.ylim(0.01, 200)
plt.tight_layout()
plt.savefig(IMG / "tcc_13nodos_bt.png", dpi=150)
plt.close()
print(f"Figura BT: {IMG / 'tcc_13nodos_bt.png'}")
