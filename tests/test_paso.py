"""API de simulacion paso a paso: iniciar/paso/iterar deben coincidir con run()."""

import numpy as np
import pytest

from bloques_crysi import (
    Modelo, FuenteConstante, FuenteSeno, FuenteTrifasica, Relay,
    MaquinaImanesPermanentes,
)


def _arma_simple(dt=1e-4):
    """PMAC alimentada: estados dinamicos + lazo estatico."""
    m = Modelo(dt=dt)
    src = m.add(FuenteTrifasica("red", 310.0, 50.0))
    tl = m.add(FuenteConstante("tl", 5.0))
    maq = m.add(MaquinaImanesPermanentes(
        "pmsm", rs=0.1, Ld=1e-3, Lq=1e-3, lam_m=0.1, P=6, J=0.01, Bm=0.001))
    m.conectar(src.salida, maq.terminales)
    m.conectar(tl.salida, maq.T_L)
    return m, maq


def test_iterar_concatena_igual_a_run():
    m, maq = _arma_simple()
    reg = [maq.sensorVelocidad(), maq.sensor3I()]
    res = m.run(0.2, registrar=reg)
    chunks = list(m.iterar(0.2, registrar=reg, chunk=0.05))
    assert 4 <= len(chunks) <= 5
    # los chunks comparten el instante del borde -> reconstruir la serie
    t = np.concatenate([chunks[0].t] + [c.t[1:] for c in chunks[1:]])
    wm = np.concatenate([chunks[0]["wm"]] + [c["wm"][1:] for c in chunks[1:]])
    I = np.concatenate([chunks[0]["I"]] + [c["I"][1:] for c in chunks[1:]])
    np.testing.assert_allclose(t, res.t)
    np.testing.assert_allclose(wm, res["wm"], rtol=1e-12, atol=0)
    np.testing.assert_allclose(I, res["I"], rtol=1e-12, atol=0)
    # memoria acotada por el chunk
    for c in chunks:
        assert len(c.t) <= 501


def test_iterar_sin_chunk_igual_a_run():
    m, maq = _arma_simple()
    reg = [maq.sensorVelocidad()]
    res = m.run(0.1, registrar=reg)
    (unico,) = list(m.iterar(0.1, registrar=reg))
    np.testing.assert_allclose(unico.t, res.t)
    np.testing.assert_allclose(unico["wm"], res["wm"])


def test_paso_por_paso_igual_a_run():
    m, maq = _arma_simple()
    reg = [maq.sensorVelocidad(), maq.sensor3I()]
    res = m.run(0.1, registrar=reg)
    vals0 = m.iniciar(registrar=reg)
    assert abs(vals0["wm"] - res["wm"][0]) < 1e-12
    for _ in range(int(round(0.1 / m.dt))):
        vals = m.paso()
    np.testing.assert_allclose(vals["wm"], res["wm"][-1], rtol=1e-12, atol=0)
    np.testing.assert_allclose(vals["I"], res["I"][-1], rtol=1e-12, atol=0)


def test_paso_requiere_iniciar():
    m = Modelo(dt=1e-3)
    m.add(FuenteConstante("f", 1.0))
    with pytest.raises(RuntimeError, match="iniciar"):
        m.paso()


def test_iterar_divergente_lanza():
    from .test_algebraico import _lazo_divergente
    m, ref, suma, k = _lazo_divergente()
    gen = m.iterar(0.02, registrar=[k], chunk=0.01)
    with pytest.raises(RuntimeError, match="no convergio"):
        next(gen)


def test_paso_divergente_lanza():
    from .test_algebraico import _lazo_divergente
    m, ref, suma, k = _lazo_divergente()
    m.iniciar(registrar=[k])
    with pytest.raises(RuntimeError, match="no convergio"):
        m.paso()


def _arma_relay(dt=1e-3):
    """Seno 1 Hz -> Relay con histeresis +/- 0.5."""
    m = Modelo(dt=dt)
    src = m.add(FuenteSeno("sn", amplitud=1.0, frecuencia=1.0))
    rly = m.add(Relay("r", umbral_on=0.5, umbral_off=-0.5))
    m.conectar(src.salida, rly.entrada)
    return m, src, rly


def _concatena(chunks):
    t = np.concatenate([chunks[0].t] + [c.t[1:] for c in chunks[1:]])
    cols = {k: np.concatenate([chunks[0][k]] + [c[k][1:] for c in chunks[1:]])
            for k in chunks[0]}
    return t, cols


def test_eventos_refinan_cruce_de_umbral():
    m, src, rly = _arma_relay()
    chunks = list(m.iterar(1.0, registrar=[src, rly], chunk=1.0,
                           eventos=[(src.salida, 0.5)], profundidad=12))
    t, cols = _concatena(chunks)
    m2, src2, rly2 = _arma_relay()
    res2 = m2.run(1.0, registrar=[src2, rly2])
    # cruces exactos de sin(2*pi*t) = 0.5: t = 1/12 y 5/12 s
    t_esperados = np.array([1 / 12, 5 / 12])
    for te in t_esperados:
        k = np.argmin(np.abs(t - te))
        assert abs(t[k] - te) < 1e-5          # resolución dt/2^12
        assert abs(cols["sn"][k] - 0.5) < 1e-4
    # las muestras de la grilla uniforme no cambian (estados restaurados)
    mask = np.isin(np.round(t, 12), np.round(res2.t, 12))
    assert mask.sum() == len(res2.t)
    assert mask.sum() + 2 == len(t)           # 2 muestras refinadas
    np.testing.assert_allclose(cols["sn"][mask], res2["sn"], rtol=1e-12, atol=0)
    np.testing.assert_allclose(cols["r"][mask], res2["r"], rtol=1e-12, atol=0)


def test_eventos_sin_cruce_no_insertan():
    m, src, rly = _arma_relay()
    chunks = list(m.iterar(1.0, registrar=[src], chunk=1.0,
                           eventos=[(src.salida, 0.99)], profundidad=12))
    t, _ = _concatena(chunks)
    # sin(2*pi*t) cruza 0.99 dos veces por periodo (sube y baja)
    assert len(t) == int(round(1.0 / m.dt)) + 1 + 2