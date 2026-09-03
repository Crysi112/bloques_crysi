"""Bateria con modelo termico: estado T, calentamiento Joule y R(T)."""

import numpy as np
import pytest

from bloques_crysi import Bateria, FuenteConstante, Modelo, Puerto
from bloques_crysi.backend_numpy import simular
from .helpers import resolver_registro


def _arma(termica=True, i=40.0, C_th=30_000.0, R_th=0.3, alpha=0.0035,
          T_amb=25.0):
    m = Modelo(dt=1.0)
    bat = m.add(Bateria("bat", termica=termica, C_th=C_th, R_th=R_th,
                        alpha=alpha, T_amb=T_amb))
    f = m.add(FuenteConstante("f", i))
    m.conectar(f.salida, bat.entrada)
    return m, bat


def test_temperatura_arranca_en_ambiente():
    m, bat = _arma(T_amb=30.0)
    res = m.run(1.0, registrar=[bat.salida])
    assert abs(res["bat"][0, 2] - 30.0) < 1e-9


def test_calentamiento_joule_y_equilibrio():
    # T' = (I^2*R(T) - (T - T_amb)/R_th)/C_th ; en equilibrio:
    # T_inf = T_amb + I^2*R(T_inf)*R_th
    m, bat = _arma(i=40.0, C_th=10_000.0, R_th=1.0, alpha=0.0, T_amb=25.0)
    res = m.run(200_000.0, registrar=[bat.salida])
    T = res["bat"][:, 2]
    R = bat.param[5]
    T_inf = 25.0 + 40.0 ** 2 * R * 1.0
    assert abs(T[-1] - T_inf) / T_inf < 1e-3
    assert T[0] == 25.0
    # monotona creciente hacia el equilibrio
    assert np.all(np.diff(T) >= -1e-9)


def test_equilibrio_joule_referencia_euler():
    # integracion Euler a mano con los mismos parametros
    m, bat = _arma(i=25.0, C_th=20_000.0, R_th=0.5, alpha=0.0, T_amb=25.0)
    res = m.run(4000.0, registrar=[bat.salida])
    T = res["bat"][:, 2]
    R = bat.param[5]
    T_ref = 25.0 + 25.0 ** 2 * R * 0.5 * (1.0 - np.exp(-4000.0 / 10_000.0))
    assert abs(T[-1] - T_ref) < 1e-4


def test_R_depende_de_temperatura():
    # alpha > 0: a temperatura alta R crece -> Vbat cae mas (misma I)
    m, bat = _arma(i=40.0, C_th=20_000.0, R_th=2.0, alpha=0.004, T_amb=25.0)
    res = m.run(200_000.0, registrar=[bat.salida])
    vbat = res["bat"][:, 0]
    T = res["bat"][:, 2]
    assert T[-1] > 25.0
    R0 = bat.param[5]
    # equilibrio: T_inf = T_amb + I^2*R(T_inf)*R_th (ecuacion lineal en T)
    T_inf = (25.0 + 40.0 ** 2 * R0 * 2.0
             * (1.0 - 0.004 * 25.0)) / (1.0 - 40.0 ** 2 * R0 * 2.0 * 0.004)
    assert abs(T[-1] - T_inf) / T_inf < 0.03
    # con la corriente constante, la caida R*i crece con R
    assert abs(vbat[0] - vbat[-1]) > 1e-3


def test_termica_false_temperatura_fija():
    m, bat = _arma(termica=False, i=40.0)
    res = m.run(5000.0, registrar=[bat.salida])
    T = res["bat"][:, 2]
    assert np.all(np.abs(T - 25.0) < 1e-9)


def test_termica_c_igual_numpy():
    m, bat = _arma(i=30.0, C_th=5_000.0, R_th=2.0, alpha=0.004, T_amb=22.0)
    res = m.run(3000.0, registrar=[bat.salida])
    m._resolver()
    rec_idx, _ = resolver_registro(m, [bat.salida])
    datos = simular(m.bloques, m.dt, 3000.0, rec_idx,
                    max_iter=m.max_iter, tol=m.tol, w_opt=m.w_opt,
                    orden_estatico=m._orden_estatico())
    for k in range(3):
        np.testing.assert_allclose(datos[k], res["bat"][:, k],
                                   rtol=1e-9, atol=1e-8)


def test_sensor_temperatura():
    m, bat = _arma(i=40.0, T_amb=28.0)
    st = bat.sensorTemperatura()
    res = m.run(1.0, registrar=[st])
    assert res["T"][0] == 28.0


def test_validaciones():
    with pytest.raises(ValueError):
        Bateria("b", termica=True, C_th=-1.0)
    with pytest.raises(ValueError):
        Bateria("b", termica=True, R_th=0.0)