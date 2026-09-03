"""Validación de sensores: C == numpy y relaciones físicas de los sensores."""

import numpy as np

from bloques_crysi import (
    Modelo, FuenteTrifasica, FuenteConstante, Scope, MedidorPotencia,
    MaquinaImanesPermanentes,
)

from .helpers import correr_c_y_numpy, flat_c


def _arma(dt):
    m = Modelo(dt=dt)
    src = m.add(FuenteTrifasica("red", 310.0, 50.0))
    tl = m.add(FuenteConstante("tl", 5.0))
    maq = m.add(MaquinaImanesPermanentes(
        "pmsm", rs=0.1, Ld=1e-3, Lq=1e-3, lam_m=0.1, P=6, J=0.01, Bm=0.001))
    m.conectar(src.salida, maq.terminales)
    m.conectar(tl.salida, maq.T_L)
    return m, maq


def _compara(modelo, t_fin, registrar, rtol=1e-9, atol=1e-10):
    res_c, arr_np, n_steps = correr_c_y_numpy(modelo, t_fin, registrar)
    arr_c = flat_c(res_c, registrar, modelo, n_steps)
    assert arr_c.shape == arr_np.shape
    np.testing.assert_allclose(arr_c, arr_np, rtol=rtol, atol=atol)


def _sensores_todos(maq):
    return [maq.sensor3V(), maq.sensor3I(), maq.sensorVelocidad(),
            maq.sensorPosicion(), maq.sensorPosicionElectrica(),
            maq.sensorPar(), maq.sensorCorrienteD(), maq.sensorCorrienteQ()]


def test_todos_los_sensores_c_numpy_euler():
    m, maq = _arma(1e-4)
    _compara(m, 0.3, _sensores_todos(maq), atol=1e-8)


def test_todos_los_sensores_c_numpy_rk4():
    m, maq = _arma(2e-4)
    m.metodo = 1
    _compara(m, 0.3, _sensores_todos(maq), atol=1e-8)


def test_nombres_canales():
    m, maq = _arma(1e-4)
    res = m.run(0.5, registrar=_sensores_todos(maq))
    # sensor de voltaje: clave "V", array (n_steps, 3) con va,vb,vc
    assert "V" in res and res["V"].shape[1] == 3
    assert "I" in res and res["I"].shape[1] == 3
    for k in ("wm", "th_rm", "th_e", "Te", "ids", "iqs"):
        assert k in res
    assert res["wm"].ndim == 1


def test_sensores_3I_consistencia():
    # corrientes de fase equilibradas: ia+ib+ic = 0 (sin neutro) en todo instante
    m, maq = _arma(1e-4)
    res = m.run(0.5, registrar=[maq.sensor3I()])
    I = np.asarray(res["I"])
    suma = I[:, 0] + I[:, 1] + I[:, 2]
    np.testing.assert_allclose(suma, 0.0, atol=1e-3)
    # desfase de 120 grados entre fases
    ia, ib = I[:, 0], I[:, 1]
    c = np.dot(ia, ib) / (np.linalg.norm(ia) * np.linalg.norm(ib) + 1e-12)
    assert abs(np.rad2deg(np.arccos(np.clip(c, -1, 1))) - 120.0) < 5.0


def test_sensores_velocidad_posicion():
    # d(th)/dt = wm => th[k] - th[k-1] ~ wm*dt (Euler, valor pre-actualizacion)
    m, maq = _arma(1e-4)
    res = m.run(0.2, registrar=[maq.sensorVelocidad(), maq.sensorPosicion()])
    wm = np.asarray(res["wm"]); th = np.asarray(res["th_rm"])
    dt = m.dt
    np.testing.assert_allclose(np.diff(th), wm[:-1] * dt, atol=1e-3)


def test_sensores_3V_amplitud():
    # fuente 310 V pico-fase; sensor lee en terminales
    m, maq = _arma(1e-4)
    res = m.run(0.1, registrar=[maq.sensor3V()])
    a = np.asarray(res["V"])[:, 0]
    assert abs(np.max(np.abs(a)) - 310.0) < 1e-6


def test_scope_hereda_metadatos_de_sensores():
    # el Scope toma ["ia","ib","ic"] del sensor3I en vez de "I[0]"...
    m, maq = _arma(1e-4)
    sc = m.add(Scope("sc", anchos=[1, 3], mostrar=False))
    m.conectar(maq.sensorVelocidad(), sc.canales[0])
    m.conectar(maq.sensor3I(), sc.canales[1])
    assert sc._etiquetas(4) == ["wm[0]", "ia", "ib", "ic"]
    # los guiones explicitos tienen prioridad sobre la metadata
    sc2 = m.add(Scope("sc2", anchos=[3], mostrar=False,
                      guiones=["a", "b", "c"]))
    m.conectar(maq.sensor3I(), sc2.canales[0])
    assert sc2._etiquetas(3) == ["a", "b", "c"]


def test_medidor_salida_lleva_canales():
    med = MedidorPotencia("med")
    assert med.salida.canales == ["P_e", "Q_e", "P_m"]
    med_dc = MedidorPotencia("med_dc", fases=1)
    assert med_dc.salida.canales == ["P_e", "P_m"]
    sc = Scope("sc", anchos=[3], mostrar=False, bloqueo=False)
    assert sc.bloqueo is False
    m = Modelo(dt=1e-4)
    m.add(med)
    m.add(sc)
    m.conectar(med.salida, sc.canales[0])
    assert sc._etiquetas(3) == ["P_e", "Q_e", "P_m"]