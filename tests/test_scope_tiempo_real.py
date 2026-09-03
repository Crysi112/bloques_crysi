"""ScopeTiempoReal: la ruta en vivo (iterar por bloques) devuelve el mismo
resultado que Modelo.run() clasico, con timestamps unicos."""

import numpy as np

from bloques_crysi import (
    Modelo, FuenteTrifasica, FuenteConstante, MaquinaInduccion,
    ScopeTiempoReal,
)

WB = 2 * np.pi * 60.0


def _arma():
    m = Modelo(dt=1e-4)
    red = m.add(FuenteTrifasica("red", amplitud=220.0 * np.sqrt(2.0 / 3.0),
                                frecuencia=60.0))
    tl = m.add(FuenteConstante("tl", 0.0))
    maq = m.add(MaquinaInduccion(
        "mi", rs=0.435, rr=0.816, Lm=26.13 / WB, Lls=0.754 / WB,
        Llr=0.754 / WB, P=4, J=0.089))
    m.conectar(red.salida, maq.terminales)
    m.conectar(tl.salida, maq.T_L)
    return m, maq


def test_scope_tiempo_real_igual_a_run():
    m, maq = _arma()
    res = m.run(0.5, registrar=[maq.sensorVelocidad(), maq.sensorPar(),
                                maq.sensor3I()])
    m2, maq2 = _arma()
    sc = m2.add(ScopeTiempoReal("sc", anchos=[1, 1, 3], mostrar=False))
    m2.conectar(maq2.sensorVelocidad(), sc.canales[0])
    m2.conectar(maq2.sensorPar(), sc.canales[1])
    m2.conectar(maq2.sensor3I(), sc.canales[2])
    res2 = m2.run(0.5, registrar=[sc, maq2.sensorVelocidad(),
                                  maq2.sensorPar(), maq2.sensor3I()])
    np.testing.assert_allclose(res2["wm"], res["wm"], rtol=1e-6, atol=1e-9)
    np.testing.assert_allclose(res2["Te"], res["Te"], rtol=1e-6, atol=1e-9)
    np.testing.assert_allclose(res2["I"], res["I"], rtol=1e-6, atol=1e-9)
    assert len(res2.t) == len(res.t) == int(round(0.5 / 1e-4)) + 1
    assert np.unique(res2.t).size == len(res2.t)


def test_scope_tiempo_real_sin_ventana_no_bloquea():
    m, maq = _arma()
    sc = m.add(ScopeTiempoReal("sc", anchos=[1], mostrar=False,
                               esperar=True))
    m.conectar(maq.sensorVelocidad(), sc.canales[0])
    res = m.run(0.05, registrar=[sc, maq.sensorVelocidad()])
    assert res["wm"][-1] > 10.0
    assert np.isfinite(res["wm"]).all()