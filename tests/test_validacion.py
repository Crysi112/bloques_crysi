"""Validación C == backend numpy (mismos coefs, misma topología)."""

import numpy as np
import pytest

from bloques_crysi import (
    Modelo, FuenteConstante, FuenteEscalon, FuenteRampa, FuenteSeno,
    FuenteTrifasica, Ganancia, Suma, Integrador, PID, Clarke, InvClarke,
    Park, InvPark, FuncionTransferencia,
    Puerto,
)

from .helpers import correr_c_y_numpy, flat_c


def par(bloque, tipo, offset, n):
    return Puerto(bloque, tipo, offset, n)


def _compara(modelo, t_fin, registrar, rtol=1e-9, atol=1e-10):
    res_c, arr_np, n_steps = correr_c_y_numpy(modelo, t_fin, registrar)
    arr_c = flat_c(res_c, registrar, modelo, n_steps)
    assert arr_c.shape == arr_np.shape
    np.testing.assert_allclose(arr_c, arr_np, rtol=rtol, atol=atol)
    return arr_c, arr_np


def test_fuentes_estatico():
    m = Modelo(dt=1e-3)
    c1 = m.add(FuenteConstante("c1", 3.5))
    g = m.add(Ganancia("g", 2.0))
    m.conectar(c1.salida, g.entrada)
    res_c, arr_np, ns = correr_c_y_numpy(m, 0.1, [g])
    arr_c = flat_c(res_c, [g], m, ns)
    np.testing.assert_allclose(arr_c, arr_np, atol=1e-12)
    assert np.allclose(arr_c, 7.0)


def test_sumador_con_signos():
    m = Modelo(dt=1e-3)
    a = m.add(FuenteConstante("a", 10.0))
    b = m.add(FuenteConstante("b", 4.0))
    s = m.add(Suma("s", (1.0, -1.0)))
    m.conectar(a.salida, par(s, "ent", 0, 1))
    m.conectar(b.salida, par(s, "ent", 1, 1))
    _compara(m, 0.05, [s], atol=1e-12)
    res = m.run(0.05, registrar=[s])
    assert np.allclose(res["s"], 6.0)


def test_integrador_rampa():
    m = Modelo(dt=1e-4)
    r = m.add(FuenteRampa("r", 2.0))
    integ = m.add(Integrador("y", 1.0))
    m.conectar(r.salida, integ.entrada)
    t_fin = 0.5
    res_c, arr_np, ns = correr_c_y_numpy(m, t_fin, [integ])
    arr_c = flat_c(res_c, [integ], m, ns)
    np.testing.assert_allclose(arr_c, arr_np, atol=1e-10)
    t = np.arange(ns) * m.dt
    np.testing.assert_allclose(arr_c.ravel(), 1.0 + t**2, atol=5e-4)


def test_tf_tustin_rampa():
    # G(s) = 1/(s+1): respuesta a rampa = t - 1 + exp(-t)
    m = Modelo(dt=1e-4)
    r = m.add(FuenteRampa("r", 1.0))
    tf = m.add(FuncionTransferencia("fr", num=[1.0], den=[1.0, 1.0]))
    m.conectar(r.salida, tf.entrada)
    res_c, arr_np, ns = correr_c_y_numpy(m, 1.0, [tf])
    arr_c = flat_c(res_c, [tf], m, ns)
    np.testing.assert_allclose(arr_c, arr_np, atol=1e-10)
    t = np.arange(ns) * m.dt
    ref = t - 1.0 + np.exp(-t)
    np.testing.assert_allclose(arr_c.ravel(), ref, rtol=2e-2, atol=1e-3)


def test_pid_pasos():
    m = Modelo(dt=1e-4)
    esc = m.add(FuenteEscalon("ref", 1.0, 0.05))
    pid = m.add(PID("pid", Kp=5.0, Ki=20.0, Kd=0.1, Tf=0.001, u_min=-10, u_max=10))
    m.conectar(esc.salida, pid.entrada)
    _compara(m, 1.0, [pid])


def test_cadena_clarke_park_inversas():
    m = Modelo(dt=1e-4)
    f3 = m.add(FuenteTrifasica("red", 310.0, 50.0))
    cl = m.add(Clarke("cl"))
    pk = m.add(Park("pk"))
    ipk = m.add(InvPark("ipk"))
    icl = m.add(InvClarke("icl"))
    th = m.add(FuenteRampa("th", 314.0))
    m.conectar(f3.salida, cl.entrada)
    m.conectar(cl.salida, par(pk, "ent", 0, 2))
    m.conectar(th.salida, par(pk, "ent", 2, 1))
    m.conectar(pk.salida, par(ipk, "ent", 0, 2))
    m.conectar(th.salida, par(ipk, "ent", 2, 1))
    m.conectar(ipk.salida, par(icl, "ent", 0, 2))
    _compara(m, 0.1, [icl])


def test_park_referencia_conocida():
    """Park de una terna balanceada con theta = wt deja un eje constante (=100)."""
    m = Modelo(dt=1e-5)
    f3 = m.add(FuenteTrifasica("red", 100.0, 50.0))
    th = m.add(FuenteRampa("th", 2 * np.pi * 50))
    cl = m.add(Clarke("cl"))
    pk = m.add(Park("pk"))
    m.conectar(f3.salida, cl.entrada)
    m.conectar(cl.salida, par(pk, "ent", 0, 2))
    m.conectar(th.salida, par(pk, "ent", 2, 1))
    res = m.run(0.04, registrar=[pk])
    d, q = res["pk"][:, 0], res["pk"][:, 1]
    # un eje queda ≈ 0 (constante) y el otro en |coordenada|=amplitud
    ejes = np.abs([d, q])
    mag = np.sqrt(d**2 + q**2)
    assert np.allclose(np.max(ejes, axis=0), 100.0, atol=1e-6)
    assert np.allclose(mag, 100.0, atol=1e-6)
    assert np.min(np.min(ejes, axis=0)) < 1e-6


def test_seno_fase():
    m = Modelo(dt=1e-4)
    s = m.add(FuenteSeno("s", 2.0, 10.0, 0.3, 1.0))
    res_c, arr_np, ns = correr_c_y_numpy(m, 0.1, [s])
    arr_c = flat_c(res_c, [s], m, ns)
    np.testing.assert_allclose(arr_c, arr_np, atol=1e-12)
    t = np.arange(ns) * m.dt
    ref = 1.0 + 2.0 * np.sin(2 * np.pi * 10 * t + 0.3)
    np.testing.assert_allclose(arr_c.ravel(), ref, atol=1e-8)