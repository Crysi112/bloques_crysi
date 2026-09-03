"""Validación de la Bateria Tremblay-Dessaint: C == numpy y física."""

import numpy as np
import pytest

from bloques_crysi import (
    Modelo, Bateria, FuenteConstante, FuenteEscalon, Suma, Puerto,
)
from .helpers import correr_c_y_numpy, flat_c


def _compara(modelo, t_fin, registrar, rtol=1e-9, atol=1e-10, metodo=0):
    modelo.metodo = metodo
    res_c, arr_np, n_steps = correr_c_y_numpy(modelo, t_fin, registrar)
    arr_c = flat_c(res_c, registrar, modelo, n_steps)
    assert arr_c.shape == arr_np.shape, (arr_c.shape, arr_np.shape)
    np.testing.assert_allclose(arr_c, arr_np, rtol=rtol, atol=atol)
    return arr_c, arr_np


def _arma_descarga(i=40.0, dt=1.0, t_fin=2000.0, **kw):
    m = Modelo(dt=dt)
    bat = m.add(Bateria("bat", **kw))
    src = m.add(FuenteConstante("i", i))
    m.conectar(src.salida, bat.entrada)
    return m, bat


def _arma_carga_descarga(t_paso=5000.0):
    m = Modelo(dt=1.0)
    bat = m.add(Bateria("bat"))
    s = m.add(Suma("sum", signos=[1.0, -1.0]))
    e = m.add(FuenteEscalon("e", 80.0, t_paso=t_paso))
    m.conectar(Puerto(bat, "sal", 0, 1), Puerto(s, "ent", 0, 1))
    m.conectar(e.salida, Puerto(s, "ent", 1, 1))
    m.conectar(s.salida, bat.entrada)
    return m, bat


# ---------------------------------------------------------- C == numpy

def test_bateria_c_numpy_euler():
    m, bat = _arma_descarga()
    _compara(m, 1000.0, [bat.salida])


def test_bateria_c_numpy_rk4():
    m, bat = _arma_descarga()
    _compara(m, 1000.0, [bat.salida], metodo=1)


def test_bateria_c_numpy_perfil_carga():
    m, bat = _arma_carga_descarga()
    _compara(m, 12000.0, [bat.salida])


# ---------------------------------------------------------- física

def test_bateria_v0_y_soc_inicial():
    # V0 = Vfull - R*i ; Vfull = Vnom*(1+f_vfull) = 13.8, R = rho*Vnom/Qrated
    m, bat = _arma_descarga(i=40.0)
    arr, _ = _compara(m, 4.0, [bat.salida], atol=1e-8)
    vbat, soc, T = arr
    Vfull = 13.8
    R = 0.0024
    assert abs(vbat[0] - (Vfull - R * 40.0)) < 1e-3
    assert abs(soc[0] - 1.0) < 1e-9
    assert soc[-1] < soc[0]          # descargando el SOC baja


def test_bateria_soc_integracion():
    # it = integral(i/3600) -> SOC = 1 - it/Q ; Q = 1.09*100 = 109
    m, bat = _arma_descarga(i=40.0)
    arr, _ = _compara(m, 2000.0, [bat.salida], atol=1e-8)
    soc = arr[1]
    esperado = 1.0 - 40.0 * 2000.0 / 3600.0 / 109.0
    assert abs(soc[-1] - esperado) < 1e-3


def test_bateria_carga_se_detiene_al_90():
    # 1000 s a 40 A (descarga) y luego 2000 s a -40 A (carga): el modelo
    # de Simulink frena la carga en it = 0.1*Q (SOC 90 %) y la tension
    # queda acotada por Vcap = 1.25*Vfull = 17.25
    m = Modelo(dt=1.0)
    bat = m.add(Bateria("bat"))
    s = m.add(Suma("sum", signos=[1.0, -1.0]))
    c = m.add(FuenteConstante("c", 40.0))
    e = m.add(FuenteEscalon("e", 80.0, t_paso=1000.0))
    m.conectar(c.salida, Puerto(s, "ent", 0, 1))
    m.conectar(e.salida, Puerto(s, "ent", 1, 1))
    m.conectar(s.salida, bat.entrada)
    arr, _ = _compara(m, 3000.0, [bat.salida], atol=1e-8)
    vbat, soc, _ = arr
    assert 0.895 < soc[-1] < 0.905
    assert np.isfinite(vbat).all()
    assert vbat.max() < 20.0
    assert vbat[-1] > 16.0           # subida de fin de carga


def test_bateria_fin_de_descarga():
    # descargando a 40 A, it llega a 0.9*Q (SOC 10 %) y se detiene,
    # con la caida abrupta de tension de fin de descarga
    m, bat = _arma_descarga(i=40.0)
    arr, _ = _compara(m, 9000.0, [bat.salida], atol=1e-8)
    vbat, soc, _ = arr
    assert 0.095 < soc[-1] < 0.105
    assert vbat[-1] < 12.0           # caida final por debajo de Vnom
    np.testing.assert_allclose(soc[-50:], soc[-1], atol=1e-4)  # ya no baja


def test_bateria_carga_acotada_por_vcap():
    # cargando desde SOCinit=95 % (it0 < 0.1*Q), la carga se frena y la
    # tension queda limitada por Vcap = 1.25*Vfull = 17.25
    m, bat = _arma_descarga(i=-40.0, SOCinit=95.0)
    arr, _ = _compara(m, 6000.0, [bat.salida], atol=1e-8)
    vbat, soc, _ = arr
    assert soc[-1] > 0.94            # no se descarga cargando
    assert vbat.max() <= 17.25 + 0.1
    assert np.isfinite(vbat).all()


# ---------------------------------------------------------- curvas

def test_curva_pasa_por_puntos_de_ajuste():
    bat = Bateria("bat")
    its, v = bat.curva()             # 1C = Qrated/3600
    i_c = 100.0 / 3600.0
    R = 0.0024
    Vfull = 13.8
    assert abs(v[0] - (Vfull - R * i_c)) < 5e-4
    idx = int(np.argmin(np.abs(its - 0.9 * bat.param[2])))
    assert abs(v[idx] - 12.0) < 1e-3     # punto nominal (Qnom, Vnom)


def test_curva_descarga_decreciente():
    bat = Bateria("tipo")
    its, v = bat.curva(n=500)
    d = np.diff(v)
    assert (d <= 1e-9).all()


def test_curva_carga_fin_limitado():
    bat = Bateria("tipo")
    its, v = bat.curva(carga=True, n=500)
    assert v.max() <= 17.25 + 1e-6    # Vcap = 1.25*Vfull
    assert np.isfinite(v).all()


# ---------------------------------------------------------- tipos

def test_tipos_y_alias():
    for nombre in ("plomo_acido", "litio", "niquel_cadmio",
                   "niquel_metal_hidruro"):
        bat = Bateria("b", tipo=nombre)
        E0, K, Q, A, B = bat.param[:5]
        assert K > 0 and A > 0 and B > 0
        assert Q >= 100.0
        assert E0 > 0
    for alias in ("lead-acid", "li", "NiMH", "plomo-ácido", "ni-cd"):
        assert Bateria("b", tipo=alias).tipo in (
            "plomo_acido", "litio", "niquel_cadmio", "niquel_metal_hidruro")


def test_tipo_invalido_y_socinit():
    with pytest.raises(ValueError):
        Bateria("b", tipo="zinc_aire")
    with pytest.raises(ValueError):
        Bateria("b", SOCinit=5.0)
    with pytest.raises(ValueError):
        Bateria("b", SOCinit=101.0)
    with pytest.raises(ValueError):
        Bateria("b", Vnom=-12.0)
    with pytest.raises(ValueError):
        Bateria("b", Qmax=50.0)


# ------------------------------------------- mejoras: eta_c e histeresis

def test_eta_c_c_numpy():
    m = Modelo(dt=1.0)
    bat = m.add(Bateria("bat", tipo="litio", SOCinit=50.0, eta_c=0.5))
    src = m.add(FuenteConstante("i", -10.0))   # carga
    m.conectar(src.salida, bat.entrada)
    _compara(m, 3600.0, [bat.salida], atol=1e-8)
    _compara(m, 3600.0, [bat.salida], atol=1e-8, metodo=1)


def test_eta_c_mitad_velocidad_de_carga():
    # con eta_c=0.5 el SOC sube a la mitad de velocidad con la misma I
    deltas = []
    for eta in (1.0, 0.5):
        m, bat = _arma_descarga(-10.0, dt=1.0, t_fin=3600.0,
                                tipo="litio", SOCinit=50.0, eta_c=eta)
        res = m.run(3600.0, registrar=[bat.salida])
        deltas.append(res["bat"][-1, 1] - 0.5)
    assert abs(deltas[1] / deltas[0] - 0.5) < 0.02, deltas


def test_histeresis_c_numpy():
    m = Modelo(dt=0.1)
    bat = m.add(Bateria("bat", tipo="plomo_acido", histeresis=True))
    src = m.add(FuenteConstante("i", 20.0))
    m.conectar(src.salida, bat.entrada)
    _compara(m, 7200.0, [bat.salida], atol=1e-7)
    _compara(m, 7200.0, [bat.salida], atol=1e-7, metodo=1)


def test_histeresis_subtensa_tras_descarga():
    # tras descarga prolongada exp_h -> A: Vbat > sin histeresis
    vs = []
    for his in (False, True):
        m = Modelo(dt=0.1)
        bat = m.add(Bateria("bat", tipo="plomo_acido", histeresis=his))
        src = m.add(FuenteConstante("i", 20.0))
        m.conectar(src.salida, bat.entrada)
        res = m.run(7200.0, registrar=[bat.salida])
        vs.append(res["bat"][-1, 0])
    assert vs[1] > vs[0] + 0.4, vs


def test_histeresis_subtensa_tras_carga():
    # en carga exp_h -> +A (termino siempre aditivo): Vbat > sin histeresis
    vs = []
    for his in (False, True):
        m = Modelo(dt=0.1)
        bat = m.add(Bateria("bat", tipo="plomo_acido", SOCinit=50.0,
                            histeresis=his))
        src = m.add(FuenteConstante("i", -20.0))
        m.conectar(src.salida, bat.entrada)
        res = m.run(3600.0, registrar=[bat.salida])
        vs.append(res["bat"][-1, 0])
    assert vs[1] > vs[0] + 0.4, vs


def test_eta_c_invalido():
    with pytest.raises(ValueError):
        Bateria("b", eta_c=0.0)
    with pytest.raises(ValueError):
        Bateria("b", eta_c=1.5)
