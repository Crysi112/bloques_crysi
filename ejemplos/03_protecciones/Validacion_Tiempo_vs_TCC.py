"""Validacion temporal: bloques 50/51/86/52 vs prediccion TCC estatica.

Parte A: el integrador de disco de Rele51 (paso a paso) debe disparar en el
         tiempo cerrado calcular_tiempo_trip para varios multiplos de pickup.
Parte B: cadena completa 50 -> 86 -> 52 ante falla senoidal; el 52 solo abre
         en cruce por cero tras su retardo mecanico.
"""

from bloques_crysi import crear_proteccion

PICKUP_51 = 19.46 * 1.25
PICKUP_50 = 138.19 * 1.7
TD_51 = 3.0

# ================================================================
# PARTE A: RELE 51 DINAMICO VS CURVA TCC
# ================================================================
print("Parte A: Rele51 paso a paso vs calcular_tiempo_trip (EI, TD=3)")
print(f"  {'m (I/Ipu)':>9s} {'predicho [s]':>12s} {'simulado [s]':>12s} {'err [%]':>8s}")
for m in (2.0, 5.0, 10.0):
    rele = crear_proteccion("51P", f"51_m{m:g}", I_pickup=PICKUP_51,
                            dial_tds=TD_51, tipo_curva="IEEE_EXTREMELY_INVERSE")
    i_falla = m * PICKUP_51
    t_pred = rele.calcular_tiempo_trip(i_falla)
    dt = 1e-3
    t, trip = 0.0, 0.0
    while trip < 0.5:
        trip, _, _ = rele.paso(i_falla, dt)
        t += dt
    err = abs(t - t_pred) / t_pred * 100.0
    print(f"  {m:>9.1f} {t_pred:>12.3f} {t:>12.3f} {err:>7.3f}")
    assert err < 0.5, f"51 no reproduce la TCC a m={m}"

r50 = crear_proteccion("50P", "50_Test", I_pickup=PICKUP_50, t_retardo=0.01)
dt = 1e-4
t, trip = 0.0, 0.0
while trip < 0.5:
    trip, _, _ = r50.paso(300.0, dt)
    t += dt
print(f"Parte A: Rele50 a 300 A disparo en {t:.4f} s (ajuste 0.0100 s)")
assert abs(t - 0.01) <= dt, "50 no respeta el retardo definido"

# ================================================================
# PARTE B: CADENA 50 -> 86 -> 52 CON FALLA SENOIDAL
# ================================================================
# Malla alineada a 60 Hz (dt = 1/6000) para que los cruces por cero
# caigan exactamente sobre muestras de simulacion.
DT = 1.0 / 6000.0
F = 60.0
T0 = 0.05
I_RMS_FALLA = 300.0
I_PICO = I_RMS_FALLA * (2.0 ** 0.5)

import math

r50b = crear_proteccion("50P", "50_Cadena", I_pickup=PICKUP_50, t_retardo=0.01)
r86b = crear_proteccion("86", "86_Cadena")
r52b = crear_proteccion("52", "52_Cadena", t_apertura_mecanica=0.04)

t = 0.0
t50 = t86 = t52 = None
while t <= 0.15:
    i_rms = I_RMS_FALLA if t >= T0 else 0.0          # estimador RMS ideal
    i_inst = I_PICO * math.sin(2.0 * math.pi * F * (t - T0)) if t >= T0 else 0.0
    trip50, _, _ = r50b.paso(i_rms, DT)
    if trip50 > 0.5 and t50 is None:
        t50 = t
    bloq, _ = r86b.paso(trip50, 0.0)
    if bloq > 0.5 and t86 is None:
        t86 = t
    estado, _, _ = r52b.paso(bloq, i_inst, DT)
    if estado < 0.5 and t52 is None:
        t52 = t
    t += DT

print(f"Parte B: 50 disparo en {t50:.4f} s (esperado 0.0600 s)")
print(f"Parte B: 86 enclavo en {t86:.4f} s (esperado 0.0600 s)")
print(f"Parte B: 52 abrio en {t52:.4f} s (esperado 0.1000 s, cruce por cero)")
assert t50 is not None and abs(t50 - 0.06) <= DT
assert t86 is not None and abs(t86 - 0.06) <= DT
assert t52 is not None and abs(t52 - 0.10) <= DT
assert r86b.estados_iniciales[0] == 1.0, "86 debe quedar enclavado"
assert r52b.estados_iniciales[0] == 0.0, "52 debe quedar abierto"
print("OK: bloques 50/51/86/52 reproducen la TCC en el dominio del tiempo")
