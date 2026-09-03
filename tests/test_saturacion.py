"""Saturacion magnetica opcional del PMAC (LUT de flujo de eje d vs Id)."""

import numpy as np
import pytest

from bloques_crysi import (
    FuenteConstante, FuenteTrifasica, MaquinaImanesPermanentes, Modelo,
    Sensor, Puerto,
)
from bloques_crysi.backend_numpy import simular
from .helpers import resolver_registro

# LUT: flujo total de eje d [Wb] vs Id [A]; flujo(0) = lam_m = 0.5
LUT = [(-30, 0.28), (-20, 0.32), (-10, 0.40), (0.0, 0.50), (10, 0.60)]

_LUT_A = np.asarray(LUT)
_PEND = np.diff(_LUT_A[:, 1]) / np.diff(_LUT_A[:, 0])


def _lds_ref(ids):
    """Inductancia incremental por tramos (con extrapolacion de extremos)."""
    ids = float(ids)
    if ids <= _LUT_A[0, 0]:
        return _PEND[0]
    for k in range(len(LUT) - 1):
        if ids <= _LUT_A[k + 1, 0]:
            return _PEND[k]
    return _PEND[-1]


def _flujo_ref(ids):
    """Flujo total de eje d por tramos (con extrapolacion de extremos)."""
    ids = np.asarray(ids, dtype=float)
    f = np.interp(ids, _LUT_A[:, 0], _LUT_A[:, 1])
    lo = ids < _LUT_A[0, 0]
    hi = ids > _LUT_A[-1, 0]
    f = np.where(lo, _LUT_A[0, 1] + _PEND[0] * (ids - _LUT_A[0, 0]), f)
    return np.where(hi, _LUT_A[-1, 1] + _PEND[-1] * (ids - _LUT_A[-1, 0]), f)


def _ref_ids(vd, t_fin, dt, rs=0.1, lut=True):
    """Euler de referencia: ids' = (vd - rs*ids)/Lds(ids) (rotor bloqueado)."""
    n = int(round(t_fin / dt))
    ids = np.empty(n + 1)
    ids[0] = 0.0
    for i in range(n):
        lds = _lds_ref(ids[i]) if lut else 1e-3
        ids[i + 1] = ids[i] + dt * (vd - rs * ids[i]) / lds
    return ids


def par(bloque, tipo, offset, n):
    return Puerto(bloque, tipo, offset, n)


def _sensor_ids(maq):
    return Sensor("ids", maq, "sal", 4, 1, canales=["ids"])


def _sensor_iqs(maq):
    return Sensor("iqs", maq, "sal", 3, 1, canales=["iqs"])


def _arma_maq(saturacion=None, v_cc=-10.0):
    """PMAC con rotor bloqueado (wm=0, th=0) y vd constante = (2/3)v_cc."""
    m = Modelo(dt=1e-5)
    v = m.add(FuenteConstante("v", v_cc))
    g0 = m.add(FuenteConstante("g0", 0.0))
    maq = m.add(MaquinaImanesPermanentes(
        "pmsm", rs=0.1, Ld=1e-3, Lq=1e-3, lam_m=0.5, P=6, J=0.01, Bm=0.001,
        mecanica_interna=False, saturacion=saturacion))
    for k in range(5):
        m.conectar(v.salida if k == 0 else g0.salida,
                   par(maq, "ent", k, 1))
    return m, maq


def _corre_y_compara(v_cc, t_fin, saturacion, lut=True):
    m, maq = _arma_maq(saturacion=saturacion, v_cc=v_cc)
    res = m.run(t_fin, registrar=[_sensor_ids(maq)])
    vd = (2.0 / 3.0) * v_cc
    ref = _ref_ids(vd, t_fin, m.dt, lut=lut)
    np.testing.assert_allclose(res["ids"], ref, rtol=1e-4, atol=1e-9)


def test_saturacion_por_defecto_es_lineal():
    _corre_y_compara(-10.0, 2e-3, None, lut=False)


def test_saturacion_usa_pendiente_de_la_lut():
    _corre_y_compara(-10.0, 2e-3, LUT)


def test_saturacion_extrapola_fuera_de_rango():
    # vd = -200 V: Id sale del rango inferior de la LUT (extrapola
    # con la pendiente del tramo extremo: 0.004 H)
    _corre_y_compara(-300.0, 2e-3, LUT)


def test_saturacion_c_igual_numpy():
    m, maq = _arma_maq(saturacion=LUT)
    reg = [_sensor_ids(maq), maq.sensorPar()]
    res = m.run(2e-3, registrar=reg)
    m._resolver()
    rec_idx, _ = resolver_registro(m, reg)
    datos = simular(m.bloques, m.dt, 2e-3, rec_idx,
                    max_iter=m.max_iter, tol=m.tol, w_opt=m.w_opt,
                    orden_estatico=m._orden_estatico())
    for k, nombre in enumerate(("ids", "Te")):
        np.testing.assert_allclose(datos[k], res[nombre],
                                   rtol=1e-12, atol=1e-9)


def test_saturacion_par_consistente_con_flujo():
    # PMAC girando con carga: Te medido debe cumplir
    # Te = 1.5*(P/2)*iqs*(flujo(ids) - Lq*ids) con el flujo de la LUT
    m = Modelo(dt=1e-5)
    red = m.add(FuenteTrifasica("red", 100.0, 10.0))
    tl = m.add(FuenteConstante("tl", 0.3))
    maq = m.add(MaquinaImanesPermanentes(
        "pmsm", rs=0.1, Ld=1e-3, Lq=1e-3, lam_m=0.5, P=6, J=0.01,
        Bm=0.001, saturacion=LUT))
    m.conectar(red.salida, maq.terminales)
    m.conectar(tl.salida, maq.T_L)
    res = m.run(0.5, registrar=[_sensor_ids(maq), _sensor_iqs(maq),
                                maq.sensorPar()])
    esperado = 1.5 * (6 / 2) * res["iqs"] * (
        _flujo_ref(res["ids"]) - 1e-3 * res["ids"])
    np.testing.assert_allclose(res["Te"], esperado, rtol=1e-9, atol=1e-9)


def test_saturacion_validacion_puntos():
    with pytest.raises(ValueError, match="creciente"):
        _arma_maq(saturacion=[(0.0, 0.5), (-1.0, 0.4)])
    with pytest.raises(ValueError, match="2 puntos"):
        _arma_maq(saturacion=[(0.0, 0.5)])