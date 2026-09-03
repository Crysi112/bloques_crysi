"""InterruptorIdeal y DiodoIdeal (conductancia variable, Peaceman-Snyder)."""

import numpy as np
import pytest

from bloques_crysi import (
    DiodoIdeal, FuenteConstante, FuenteSeno, InterruptorIdeal, Modelo,
    PulsoRectangular,
)
from bloques_crysi.backend_numpy import simular
from .helpers import resolver_registro


def par(bloque, tipo, offset, n):
    """Sub-puerto: porcion continua [offset, offset+n) de un puerto del bloque."""
    from bloques_crysi import Puerto
    return Puerto(bloque, tipo, offset, n)


def _arma_interruptor(g=1.0, dt=1e-4):
    """10 V -> interruptor; I_sw a tierra (V_c = 0)."""
    m = Modelo(dt=dt)
    v = m.add(FuenteConstante("v", 10.0))
    gnd = m.add(FuenteConstante("gnd", 0.0))
    gb = m.add(FuenteConstante("g", g))
    sw = m.add(InterruptorIdeal("sw", R_on=1e-3, R_off=1e6))
    m.conectar(gb.salida, sw.control)
    m.conectar(v.salida, par(sw, "ent", 1, 1))
    m.conectar(gnd.salida, par(sw, "ent", 2, 1))
    return m, sw


def test_interruptor_abre_y_cierra():
    m, sw = _arma_interruptor(g=1.0)
    res = m.run(0.01, registrar=[sw])
    np.testing.assert_allclose(res["sw"][:, 0], 10.0 / 1e-3, rtol=1e-9)
    m2, sw2 = _arma_interruptor(g=0.0)
    res2 = m2.run(0.01, registrar=[sw2])
    np.testing.assert_allclose(res2["sw"][:, 0], 10.0 / 1e6, rtol=1e-9)
    np.testing.assert_allclose(res2["sw"][:, 1], 10.0, rtol=1e-9)


def test_interruptor_pwm_promedio():
    m = Modelo(dt=1e-4)
    v = m.add(FuenteConstante("v", 12.0))
    gnd = m.add(FuenteConstante("gnd", 0.0))
    pwm = m.add(PulsoRectangular("pwm", amplitud=1.0, periodo=0.01,
                                 duty=0.5))
    sw = m.add(InterruptorIdeal("sw", R_on=1e-3, R_off=1e6))
    m.conectar(pwm.salida, sw.control)
    m.conectar(v.salida, par(sw, "ent", 1, 1))
    m.conectar(gnd.salida, par(sw, "ent", 2, 1))
    res = m.run(0.5, registrar=[sw])
    # 50 % duty sobre R_on => promedio 0.5*12/R_on (cambio de estado
    # en el instante del paso; el error relativo queda < 1 %)
    prom = res["sw"][:, 0].mean()
    esperado = 0.5 * 12.0 / 1e-3
    assert abs(prom - esperado) / esperado < 1e-2


def test_interruptor_numpy_igual_a_c():
    m, sw = _arma_interruptor(g=1.0)
    res = m.run(0.02, registrar=[sw])
    m._resolver()
    rec_idx, _ = resolver_registro(m, [sw])
    datos = simular(m.bloques, m.dt, 0.02, rec_idx,
                    max_iter=m.max_iter, tol=m.tol, w_opt=m.w_opt,
                    orden_estatico=m._orden_estatico())
    np.testing.assert_allclose(datos[0], res["sw"][:, 0], rtol=1e-12, atol=1e-9)


def _arma_diodo(dt=1e-4):
    """Seno 10 V@50 Hz -> diodo (V_c = 0)."""
    m = Modelo(dt=dt)
    v = m.add(FuenteSeno("v", amplitud=10.0, frecuencia=50.0))
    gnd = m.add(FuenteConstante("gnd", 0.0))
    d = m.add(DiodoIdeal("d", R_on=1e-3, R_off=1e6))
    m.conectar(v.salida, par(d, "ent", 0, 1))
    m.conectar(gnd.salida, par(d, "ent", 1, 1))
    return m, v, d


def test_diodo_conduce_solo_positivo():
    m, v, d = _arma_diodo()
    res = m.run(0.06, registrar=[v, d])
    i, vs = res["d"][:, 0], res["d"][:, 1]
    # no invierte de forma descontrolada: a lo sumo la banda de histeresis
    assert (i > -1.0 - 1e-9).all()
    pos = vs > 1e-3
    np.testing.assert_allclose(i[pos], vs[pos] / 1e-3, rtol=1e-6)
    neg = vs < -1e-3                       # fuera de la banda: bloquea
    np.testing.assert_allclose(i[neg], vs[neg] / 1e6, rtol=1e-6, atol=1e-9)


def test_diodo_histeresis_retencion():
    # V dentro de la banda de histeresis: conserva el estado previo
    m = Modelo(dt=1e-4)
    v = m.add(FuenteConstante("v", 5e-4))       # < V_f + h = 1e-3
    gnd = m.add(FuenteConstante("gnd", 0.0))
    d = m.add(DiodoIdeal("d", R_on=1e-3, R_off=1e6, V_f=0.0,
                         histeresis=1e-3))
    m.conectar(v.salida, par(d, "ent", 0, 1))
    m.conectar(gnd.salida, par(d, "ent", 1, 1))
    res = m.run(0.01, registrar=[d])
    np.testing.assert_allclose(res["d"][:, 0], 0.0, atol=1e-6)  # queda abierto
    # y con V por encima de la banda conduce
    m2 = Modelo(dt=1e-4)
    v2 = m2.add(FuenteConstante("v2", 2e-3))
    gnd2 = m2.add(FuenteConstante("gnd2", 0.0))
    d2 = m2.add(DiodoIdeal("d2", R_on=1e-3, R_off=1e6))
    m2.conectar(v2.salida, par(d2, "ent", 0, 1))
    m2.conectar(gnd2.salida, par(d2, "ent", 1, 1))
    res2 = m2.run(0.01, registrar=[d2])
    # la muestra t=0 emite con el estado inicial (abierto): desde t=dt conduce
    np.testing.assert_allclose(res2["d2"][1:, 0], 2.0, rtol=1e-9)


def test_diodo_numpy_igual_a_c():
    m, v, d = _arma_diodo()
    res = m.run(0.04, registrar=[d])
    m._resolver()
    rec_idx, _ = resolver_registro(m, [d])
    datos = simular(m.bloques, m.dt, 0.04, rec_idx,
                    max_iter=m.max_iter, tol=m.tol, w_opt=m.w_opt,
                    orden_estatico=m._orden_estatico())
    np.testing.assert_allclose(datos[0], res["d"][:, 0], rtol=1e-12, atol=1e-9)


def test_interruptor_en_lazo_estatico():
    """ref -> Suma -> interruptor (g=1) -> Ganancia -> Suma: el interruptor
    participa del lazo algebraico; I_sw = ref/(R_on + G)."""
    from bloques_crysi import Ganancia, Suma
    m = Modelo(dt=1e-4)
    ref = m.add(FuenteConstante("ref", 10.0))
    gnd = m.add(FuenteConstante("gnd", 0.0))
    suma = m.add(Suma("err", (1.0, -1.0)))
    sw = m.add(InterruptorIdeal("sw", R_on=1.0, R_off=1e6))
    g = m.add(Ganancia("g", 0.5))
    gb = m.add(FuenteConstante("gb", 1.0))
    m.conectar(ref.salida, par(suma, "ent", 0, 1))
    m.conectar(g.salida, par(suma, "ent", 1, 1))
    m.conectar(suma.salida, par(sw, "ent", 1, 1))
    m.conectar(gnd.salida, par(sw, "ent", 2, 1))
    m.conectar(gb.salida, sw.control)
    m.conectar(par(sw, "sal", 0, 1), par(g, "ent", 0, 1))
    res = m.run(0.01, registrar=[sw])
    # la muestra t=0 emite antes del primer paso del lazo algebraico
    np.testing.assert_allclose(res["sw"][1:, 0], 10.0 / 1.5, rtol=1e-6)
    np.testing.assert_allclose(res["sw"][1:, 1], 10.0 / 1.5, rtol=1e-6)