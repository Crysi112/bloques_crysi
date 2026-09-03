"""Refinamiento genuino de eventos (zero-crossing) en iterar.

Item 1: el paso se parte de verdad en t_cruce (minipaso, conmutacion
en t_cruce con dt=0, segundo minipaso).
Item 2: la foto guardada incluye senales y relojes t_next_update, asi
que las sondas y sub-pasos no corrompen los bloques multirate.
"""

import numpy as np
import pytest

from bloques_crysi import (
    FuenteRampa, Integrador, Modelo, PID, Relay,
)


def _corre(m, t_fin, eventos, registrar=None):
    """Todos los chunks de iterar concatenados (bordes compartidos)."""
    if registrar is None:
        registrar = [m.bloques[-1].salida]
    chunks = list(m.iterar(t_fin, chunk=t_fin, registrar=registrar,
                           eventos=eventos, profundidad=12))
    t = np.concatenate([chunks[0].t] + [c.t[1:] for c in chunks[1:]])
    cols = {k: np.concatenate([chunks[0][k]] + [c[k][1:] for c in chunks[1:]])
            for k in chunks[0]}
    return t, cols


def _arma_relay_rampa(m, t_inicio):
    """Rampa u = t - t_inicio -> Relay que conmuta en t = 5 + t_inicio."""
    rampa = m.add(FuenteRampa("r", pendiente=1.0, t_inicio=t_inicio))
    relay = m.add(Relay("sw", umbral_on=5.0, umbral_off=0.0))
    m.conectar(rampa.salida, relay.entrada)
    return rampa, relay


def test_evento_refina_el_cruce_en_el_tiempo_correcto():
    # dt=1: el cruce (t=5.7) cae dentro del paso [5, 6]; la muestra
    # refinada debe caer en t=5.7 +/- dt/2^12, sin duplicar la grilla
    m = Modelo(dt=1.0)
    _, relay = _arma_relay_rampa(m, 0.7)
    res = _corre(m, 10.0, [(relay.salida, 0.5)])
    ts = np.asarray(res[0])
    muestras = np.asarray(res[1]["sw"]).ravel()
    assert 5.0 in ts and 6.0 in ts
    extra = ts[(ts > 5.0) & (ts < 6.0)]
    assert len(extra) == 1
    assert abs(extra[0] - 5.7) < 1e-2
    assert len(ts) == len(np.unique(ts))          # sin tiempos repetidos
    assert muestras[5] == 0.0 and muestras[7] == 1.0  # grilla intacta


def test_evento_no_adelanta_relojes_multirate():
    # PID con Ts=2: la biseccion y los sub-pasos NO deben corromper
    # t_next_update (antes de la correccion perdia la actualizacion
    # de t=6 y el PID conmutaba una vez menos)
    m = Modelo(dt=1.0)
    _, relay = _arma_relay_rampa(m, 0.7)
    pid = m.add(PID("pid", Kp=1.0, Ki=0.1, Ts=2.0))
    m.conectar(relay.salida, pid.entrada)
    res = _corre(m, 10.0, [(relay.salida, 0.5)], registrar=[pid.salida])
    salidas = np.asarray(res[1]["pid"]).ravel()
    cambios = int(np.sum(np.abs(np.diff(salidas)) > 1e-12))
    assert cambios == 3, (salidas, cambios)


def test_evento_no_perturba_estados_de_grilla():
    # integrador lineal: los sub-pasos en t_cruce no alteran los
    # estados en los instantes de la grilla (coinciden con run())
    m = Modelo(dt=0.1)
    _, relay = _arma_relay_rampa(m, 0.65)
    integ = m.add(Integrador("x"))
    m.conectar(relay.salida, integ.entrada)
    res = _corre(m, 6.0, [(relay.salida, 0.5)], registrar=[integ.salida])
    m2 = Modelo(dt=0.1)
    _, relay2 = _arma_relay_rampa(m2, 0.65)
    integ2 = m2.add(Integrador("x"))
    m2.conectar(relay2.salida, integ2.entrada)
    res2 = m2.run(6.0, registrar=[integ2.salida])
    # comparar solo los instantes de la grilla (la muestra refinada va aparte)
    t = np.asarray(res[0])
    grilla = np.isin(np.round(t, 12), np.round(res2.t, 12))
    assert grilla.sum() == len(res2.t)
    np.testing.assert_allclose(np.asarray(res[1]["x"]).ravel()[grilla],
                               res2["x"].ravel(),
                               rtol=1e-12, atol=1e-12)


def test_evento_varias_conmutaciones():
    # subida y bajada: cada cruce de la salida del relay se refina
    # (el relay se apaga en u <= umbral_off = 0, o sea t = 10)
    m = Modelo(dt=1.0)
    rampa = m.add(FuenteRampa("r", pendiente=1.0, t_inicio=0.0))
    relay = m.add(Relay("sw", umbral_on=5.0, umbral_off=0.0))
    m.conectar(rampa.salida, relay.entrada)
    m2 = Modelo(dt=1.0)
    r2 = m2.add(FuenteRampa("r2", pendiente=-1.0, t_inicio=0.0, offset=10.0))
    relay2 = m2.add(Relay("sw2", umbral_on=5.0, umbral_off=0.0))
    m2.conectar(r2.salida, relay2.entrada)
    res = _corre(m, 10.0, [(relay.salida, 0.5)])
    res2 = _corre(m2, 15.0, [(relay2.salida, 0.5)])
    ts = np.asarray(res[0])
    ts2 = np.asarray(res2[0])
    assert len(ts) == 11 + 1                     # 11 de grilla + 1 cruce
    extra = ts[(ts > 4.0) & (ts < 5.0)]          # flip en u>=5 -> t=5.0
    assert len(extra) == 1 and abs(extra[0] - 5.0) < 1e-2
    assert len(ts2) == 16 + 1                    # bajada en t=10
    extra2 = ts2[(ts2 > 9.0) & (ts2 < 10.0)]
    assert len(extra2) == 1 and abs(extra2[0] - 10.0) < 1e-2