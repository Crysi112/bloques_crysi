"""Validación del PLL trifásico: C == numpy y enganche de fase/frecuencia."""

import numpy as np
import pytest

from bloques_crysi import Modelo, FuenteTrifasica, PLLTrifasico

from .helpers import correr_c_y_numpy, flat_c

W_FF = 2 * np.pi * 50.0


def _modelo_pll(dt=1e-4, Kp=5.0, Ki=50.0, f_ff=50.0, theta0=0.0, t_fin=1.0):
    m = Modelo(dt=dt)
    ft = m.add(FuenteTrifasica("ft", amplitud=311.0, frecuencia=50.0))
    pll = m.add(PLLTrifasico("pll", Kp=Kp, Ki=Ki, f_ff=f_ff, theta0=theta0))
    m.conectar(ft.salida, pll.entrada)
    return m, pll


def test_pll_c_numpy():
    m, pll = _modelo_pll()
    res_c, arr_np, n_steps = correr_c_y_numpy(m, 1.0, [pll.salida])
    arr_c = flat_c(res_c, [pll.salida], m, n_steps)
    np.testing.assert_array_equal(arr_c, arr_np)


def test_pll_engancha_frecuencia():
    m, pll = _modelo_pll()
    res = m.run(1.0, registrar=[pll.salida])
    t, w, th = res.t, res["pll"][:, 0], res["pll"][:, 1]
    # w -> 2*pi*50 rad/s
    assert abs(w[-1] - W_FF) < 0.1
    # th avanza a w_ff (paso estable, sin deriva)
    dth = np.diff(th) / m.dt
    assert np.max(np.abs(dth[-1000:] - W_FF)) < 1e-3
    # fase enganchada: th = wt - pi/2 (con el offset de un paso del
    # registro post-update: th_s = wt(t_{s+1}) - pi/2)
    ideal = W_FF * (t + m.dt) - np.pi / 2
    assert np.max(np.abs(th[-1000:] - ideal[-1000:])) < 1e-3


def test_pll_theta0_inicial():
    m, pll = _modelo_pll(theta0=1.0)
    res = m.run(1.0, registrar=[pll.salida])
    w, th = res["pll"][:, 0], res["pll"][:, 1]
    assert abs(w[-1] - W_FF) < 0.1
    ideal = W_FF * (res.t + m.dt) - np.pi / 2
    assert np.max(np.abs(th[-1000:] - ideal[-1000:])) < 1e-3


def test_pll_gains_lentos():
    m, pll = _modelo_pll(Kp=1.0, Ki=10.0, t_fin=2.0)
    res = m.run(2.0, registrar=[pll.salida])
    w = res["pll"][:, 0]
    assert abs(w[-1] - W_FF) < 0.5


def test_pll_f_ff_negativa():
    with pytest.raises(ValueError):
        PLLTrifasico("pll", f_ff=0.0)