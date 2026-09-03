"""Validación del PanelSolar: C == numpy, física y curvas IV."""

import numpy as np
import pytest

from bloques_crysi import Modelo, PanelSolar, FuenteConstante, Puerto
from .helpers import correr_c_y_numpy, flat_c


def _arma_panel(V=20.0, G=1000.0, T=25.0, **kw):
    m = Modelo(dt=1e-3)
    pan = m.add(PanelSolar("pan", **kw))
    g = m.add(FuenteConstante("g", G))
    t = m.add(FuenteConstante("t", T))
    v = m.add(FuenteConstante("v", V))
    m.conectar(g.salida, Puerto(pan, "ent", 0, 1))
    m.conectar(t.salida, Puerto(pan, "ent", 1, 1))
    m.conectar(v.salida, Puerto(pan, "ent", 2, 1))
    return m, pan


def _compara(m, t_fin, registrar, rtol=1e-9, atol=1e-9):
    res_c, arr_np, n_steps = correr_c_y_numpy(m, t_fin, registrar)
    arr_c = flat_c(res_c, registrar, m, n_steps)
    np.testing.assert_allclose(arr_c, arr_np, rtol=rtol, atol=atol)


# ---------------------------------------------------------- C == numpy

def test_panel_c_numpy():
    for V, G, T in ((0.0, 1000.0, 25.0), (20.0, 1000.0, 25.0),
                    (37.6, 1000.0, 25.0), (30.0, 500.0, 50.0),
                    (10.0, 100.0, 0.0)):
        m, pan = _arma_panel(V=V, G=G, T=T)
        _compara(m, 0.1, [pan.salida])


def test_panel_curva_con_el_simulador():
    # la curvaIV() (numpy puro) coincide con el bloque en el modelo
    m, pan = _arma_panel(V=20.0)
    Vs, Is, _ = pan.curvaIV()
    res = m.run(0.1, registrar=[pan.salida])
    assert abs(res["pan"][-1] - np.interp(20.0, Vs, Is)) < 1e-5


# ---------------------------------------------------------- física

def test_panel_cortocircuito_y_circuito_abierto():
    pan = PanelSolar("pan")
    # I(V=0) ~ Isc ; I(V=Voc) ~ 0
    Isc = _resolver_panel_sim(pan, 0.0)
    Ioc = _resolver_panel_sim(pan, 37.6)
    assert abs(Isc - 9.12) / 9.12 < 0.02, Isc
    assert abs(Ioc) < 1e-6, Ioc


def test_panel_curva_mpp():
    pan = PanelSolar("pan")
    Vs, Is, Ps = pan.curvaIV()
    k = np.argmax(Ps)
    assert Ps[k] > 0.95 * pan._datos["Vmp"] * pan._datos["Imp"], Ps[k]
    assert abs(Vs[k] - 29.9) < 1.5, Vs[k]
    assert abs(Is[k] - 8.63) < 0.3, Is[k]
    # potencia nula en corto y en abierto (evaluando exactamente V=0, V=Voc)
    assert abs(Ps[0]) < 1e-6
    assert abs(_resolver_panel_sim(pan, 37.6)) < 1e-6


def test_panel_crece_con_irradiancia():
    m, pan = _arma_panel(V=20.0)
    res = m.run(0.1, registrar=[pan.salida])
    i1 = res["pan"][-1]
    m2, pan2 = _arma_panel(V=20.0, G=500.0)
    res2 = m2.run(0.1, registrar=[pan2.salida])
    i2 = res2["pan"][-1]
    assert 0.48 < i2 / i1 < 0.52, (i1, i2)


def test_panel_temperatura():
    # a igual V y G, con ki > 0 la corriente crece con T (Iph sube)
    m, pan = _arma_panel(V=20.0, T=25.0)
    res = m.run(0.1, registrar=[pan.salida])
    i25 = res["pan"][-1]
    m2, pan2 = _arma_panel(V=20.0, T=60.0)
    res2 = m2.run(0.1, registrar=[pan2.salida])
    i60 = res2["pan"][-1]
    assert i60 > i25 + 0.05, (i25, i60)


def test_panel_sin_luz():
    m, pan = _arma_panel(V=0.0, G=0.0)
    res = m.run(0.1, registrar=[pan.salida])
    assert abs(res["pan"][-1]) < 1e-6


def test_panel_voc_derateo_termico():
    # a 50 C el Voc cae ~ -0.3%/C (25 C -> 37.6 V): la curva pasa por el
    # Voc deratado y I(V=Voc_STC) es fuertemente negativa
    pan = PanelSolar("pan")
    T = 50.0
    Voc_T = 37.6 * (1.0 - 0.003 * (T - 25.0))
    Ioc_T = _resolver_panel_sim(pan, Voc_T, T=T)
    Ioc_25 = _resolver_panel_sim(pan, 37.6, T=T)
    assert abs(Ioc_T) < 0.1, Ioc_T
    assert Ioc_25 < -5.0, Ioc_25  # 37.6 V ya no es circuito abierto a 50 C


def test_panel_voltaje_alto_finito():
    # un V externo absurdo no debe producir inf/NaN (clamp de exp en el
    # Newton-Raphson); solo un valor finito (aunque fuera de rango)
    pan = PanelSolar("pan")
    for V in (100.0, 1000.0, 1e6):
        I = _resolver_panel_sim(pan, V)
        assert np.isfinite(I), (V, I)


def test_panel_parametros_invalidos():
    with pytest.raises(ValueError):
        PanelSolar("p", Ns=0)
    with pytest.raises(ValueError):
        PanelSolar("p", Isc=-1.0)
    with pytest.raises(ValueError):
        PanelSolar("p", n=0.0)
    with pytest.raises(ValueError):
        PanelSolar("p", ki=-0.1)


def _resolver_panel_sim(pan, V, T=25.0):
    """Corriente del panel a V fijo via el simulador (C)."""
    m = Modelo(dt=1e-3)
    m.add(pan)
    g = m.add(FuenteConstante("g", 1000.0))
    t = m.add(FuenteConstante("t", T))
    v = m.add(FuenteConstante("v", V))
    m.conectar(g.salida, Puerto(pan, "ent", 0, 1))
    m.conectar(t.salida, Puerto(pan, "ent", 1, 1))
    m.conectar(v.salida, Puerto(pan, "ent", 2, 1))
    res = m.run(0.1, registrar=[pan.salida])
    return res["pan"][-1]