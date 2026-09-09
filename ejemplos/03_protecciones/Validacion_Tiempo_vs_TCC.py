"""Validacion temporal: bloques 50/51/86/52 en el motor vs TCC estatica.

Parte A: tres 51 (m = 2, 5, 10) y un 50 sobre fuentes constantes; el tiempo
         de disparo simulado debe igualar calcular_tiempo_trip.
Parte B: cadena 50 -> 86 -> 52 ante escalon de corriente; el 52 solo abre
         en cruce por cero tras su retardo mecanico.
"""

import numpy as np

from bloques_crysi import (
    Modelo, FuenteConstante, FuenteEscalon, FuenteSeno, crear_proteccion,
)

PICKUP_51 = 19.46 * 1.25
PICKUP_50 = 138.19 * 1.7
I_FALLA = 5 * PICKUP_51

# ================================================================
# PARTE A: 51 DINAMICO VS CURVA TCC + 50 EXACTO
# ================================================================
modelo_a = Modelo(dt=1e-3)
with modelo_a:
    src = modelo_a.add(FuenteConstante("i_falla", I_FALLA))
    src50 = modelo_a.add(FuenteConstante("i_50", 300.0))
    reles = {m: modelo_a.add(crear_proteccion(
        "51P", f"51_m{m:g}", I_pickup=I_FALLA / m,
        dial_tds=3.0, tipo_curva="IEEE_EXTREMELY_INVERSE")) for m in (2.0, 5.0, 10.0)}
    r50 = modelo_a.add(crear_proteccion("50P", "50_Test", I_pickup=PICKUP_50, t_retardo=0.01))
    for rele in list(reles.values()):
        modelo_a.conectar(src.salida, rele.entrada)
    modelo_a.conectar(src50.salida, r50.entrada)
    scope_a = modelo_a.scope("TCC", *[reles[m].salida[0] for m in (2.0, 5.0, 10.0)], r50.salida[0])

res_a = modelo_a.run(t_fin=29.0)
print("Parte A: Rele51 en motor vs calcular_tiempo_trip (EI, TD=3)")
print(f"  {'m (I/Ipu)':>9s} {'predicho [s]':>12s} {'simulado [s]':>12s} {'err [%]':>8s}")
for k, m in enumerate((2.0, 5.0, 10.0)):
    t_pred = reles[m].calcular_tiempo_trip(I_FALLA)
    t_sim = res_a.t[np.argmax(res_a["TCC"][:, k] > 0.5)]
    err = abs(t_sim - t_pred) / t_pred * 100.0
    print(f"  {m:>9.1f} {t_pred:>12.3f} {t_sim:>12.3f} {err:>7.3f}")
    assert err < 0.5, f"51 no reproduce la TCC a m={m}"

t50 = res_a.t[np.argmax(res_a["TCC"][:, 3] > 0.5)]
print(f"Parte A: Rele50 a 300 A disparo en {t50:.4f} s (ajuste 0.0100 s)")
assert abs(t50 - 0.01) <= 1e-3, "50 no respeta el retardo definido"

# ================================================================
# PARTE B: CADENA 50 -> 86 -> 52 (MALLA ALINEADA A 60 HZ)
# ================================================================
DT = 1.0 / 6000.0
modelo_b = Modelo(dt=DT)
with modelo_b:
    step = modelo_b.add(FuenteEscalon("falla", valor_final=300.0, t_paso=0.05, valor_inicial=0.0))
    seno = modelo_b.add(FuenteSeno("ired", amplitud=300.0 * 2.0 ** 0.5, frecuencia=60.0))
    cero = modelo_b.add(FuenteConstante("cero", 0.0))
    r50b = modelo_b.add(crear_proteccion("50P", "50_Cad", I_pickup=PICKUP_50, t_retardo=0.01))
    r86b = modelo_b.add(crear_proteccion("86", "86_Cad"))
    r52b = modelo_b.add(crear_proteccion("52", "52_Cad", t_apertura_mecanica=0.04))
    modelo_b.conectar(step.salida, r50b.entrada)
    modelo_b.conectar(r50b.salida[0], r86b.entrada_trip)
    modelo_b.conectar(cero.salida, r86b.entrada_reset)
    modelo_b.conectar(r86b.salida[0], r52b.entrada_trip)
    modelo_b.conectar(seno.salida, r52b.entrada_i)
    scope_b = modelo_b.scope("Cad", r50b.salida[0], r86b.salida[0], r52b.salida[0])

res_b = modelo_b.run(t_fin=0.15)
d = res_b["Cad"]
t50b = res_b.t[np.argmax(d[:, 0] > 0.5)]
t86b = res_b.t[np.argmax(d[:, 1] > 0.5)]
t52b = res_b.t[np.where((d[:, 2] < 0.5) & (res_b.t > 0.09))[0][0]]
print(f"Parte B: 50 disparo en {t50b:.4f} s, 86 enclavo en {t86b:.4f} s (esperado 0.0600 s)")
print(f"Parte B: 52 abrio en {t52b:.4f} s (esperado 0.1000 s, cruce por cero)")
assert abs(t50b - 0.06) <= 2 * DT and abs(t86b - 0.06) <= 2 * DT
assert abs(t52b - 0.10) <= 2 * DT
print("OK: bloques 50/51/86/52 reproducen la TCC en el dominio del tiempo")
