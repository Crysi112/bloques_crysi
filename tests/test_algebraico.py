"""Lazos de realimentacion (algebraicos) con PID y validacion C == numpy."""

import numpy as np
import pytest

from bloques_crysi import (
    Modelo, FuenteConstante, FuenteTrifasica, Suma, Ganancia, PID,
    MaquinaImanesPermanentes,
)
from bloques_crysi.backend_numpy import simular
from .helpers import correr_c_y_numpy, flat_c, resolver_registro


def par(bloque, tipo, offset, n):
    """Sub-puerto: porcion continua [offset, offset+n) de un puerto del bloque."""
    from bloques_crysi import Puerto
    return Puerto(bloque, tipo, offset, n)


def _lazo_estatico(G=0.5):
    """y = G*(ref - y)  =>  y = G*ref/(1+G). Lazo 100% estatico."""
    m = Modelo(dt=1e-3)
    ref = m.add(FuenteConstante("ref", 2.0))
    suma = m.add(Suma("err", (1.0, -1.0)))
    k = m.add(Ganancia("k", G))
    m.conectar(ref.salida, par(suma, "ent", 0, 1))
    m.conectar(k.salida, par(suma, "ent", 1, 1))
    m.conectar(suma.salida, k.entrada)
    return m, ref, suma, k


def test_lazo_estatico_simple():
    m, ref, suma, k = _lazo_estatico(0.5)
    res = m.run(0.05, registrar=[k])
    y = float(res["k"][-1])
    assert abs(y - 0.5 * 2.0 / 1.5) < 1e-6


def test_lazo_estatico_c_numpy():
    m, ref, suma, k = _lazo_estatico(0.5)
    reg = [ref, suma, k]
    res_c, arr_np, ns = correr_c_y_numpy(m, 0.05, reg)
    arr_c = flat_c(res_c, reg, m, ns)
    assert arr_c.shape == arr_np.shape
    np.testing.assert_allclose(arr_c, arr_np, rtol=1e-8, atol=1e-9)


def _arma_pid_velocidad(dt, wref=60.0):
    """Lazo de velocidad con PMAC (par de control via T_L) y terminales a neutro.

    Estructura: ref - wm -> PID -> Te* (a traves de T_L). El control del par
    se inyecta como par de carga electrico en la PMAC para cerrar el lazo
    cinetico (valida el lazo algebraico PID->realimentacion de velocidad).
    """
    m = Modelo(dt=dt)
    ref = m.add(FuenteConstante("ref", wref))
    suma = m.add(Suma("err", (1.0, -1.0)))
    pid = m.add(PID("pi", Kp=20.0, Ki=60.0, u_min=-500.0, u_max=500.0))
    inv = m.add(Ganancia("inv", -1.0))
    maq = m.add(MaquinaImanesPermanentes(
        "pmsm", rs=10.0, Ld=5e-3, Lq=5e-3, lam_m=0.5, P=4, J=0.05, Bm=0.01))
    neutro = m.add(FuenteTrifasica("neut", 0.0, 50.0))
    m.conectar(neutro.salida, maq.terminales)   # terminales a 0 V
    m.conectar(ref.salida, par(suma, "ent", 0, 1))
    m.conectar(maq.sensorVelocidad(), par(suma, "ent", 1, 1))
    m.conectar(suma.salida, pid.entrada)
    m.conectar(pid.salida, inv.entrada)
    m.conectar(inv.salida, maq.T_L)
    return m, maq, pid


def test_pid_lazo_c_numpy():
    m, maq, pid = _arma_pid_velocidad(1e-4, 60.0)
    reg = [maq.sensorVelocidad(), maq.sensorPar(), pid]
    res_c, arr_np, ns = correr_c_y_numpy(m, 0.3, reg)
    arr_c = flat_c(res_c, reg, m, ns)
    assert arr_c.shape == arr_np.shape
    np.testing.assert_allclose(arr_c, arr_np, rtol=1e-8, atol=1e-9)


def test_pid_alcanza_regimen():
    m, maq, pid = _arma_pid_velocidad(1e-4, 60.0)
    res = m.run(1.5, registrar=[maq.sensorVelocidad(), maq.sensorPar()])
    assert abs(res["wm"][-1] - 60.0) < 2.0


def test_pid_anti_windup_saturado():
    # referencia muy por encima de lo alcanzable: el PID debe saturarse
    m, maq, pid = _arma_pid_velocidad(1e-4, 500.0)
    res = m.run(2.0, registrar=[pid])
    u = np.asarray(res["pi"])
    assert np.max(u) <= 500.0 + 1e-6
    assert np.min(u) >= -500.0 - 1e-6


def _lazo_divergente():
    """y = 1 + 2*y (Ganancia 2.0 realimentada a la suma): GS diverge."""
    m = Modelo(dt=1e-3)
    ref = m.add(FuenteConstante("ref", 1.0))
    suma = m.add(Suma("err", (1.0, 1.0)))
    k = m.add(Ganancia("k", 2.0))
    m.conectar(ref.salida, par(suma, "ent", 0, 1))
    m.conectar(k.salida, par(suma, "ent", 1, 1))
    m.conectar(suma.salida, k.entrada)
    return m, ref, suma, k


def test_lazo_divergente_lanza_error_c():
    m, ref, suma, k = _lazo_divergente()
    with pytest.raises(RuntimeError, match="no convergio"):
        m.run(0.05, registrar=[k])


def test_lazo_divergente_lanza_error_numpy():
    m, ref, suma, k = _lazo_divergente()
    m._resolver()
    rec_idx, _ = resolver_registro(m, [k])
    with pytest.raises(RuntimeError, match="no convergio"):
        simular(m.bloques, m.dt, 0.05, rec_idx,
                max_iter=m.max_iter, tol=m.tol, w_opt=m.w_opt,
                orden_estatico=m._orden_estatico())


# ---------------------------------------------------- feedthrough (1.3)

def _lazo_tf(num, den):
    """ref -> suma -> TF -> Ganancia -> suma (lazo cerrado)."""
    from bloques_crysi import FuncionTransferencia
    m = Modelo(dt=1e-3)
    ref = m.add(FuenteConstante("ref", 1.0))
    suma = m.add(Suma("err", (1.0, -1.0)))
    tf = m.add(FuncionTransferencia("tf", num, den))
    g = m.add(Ganancia("g", 1.0))
    m.conectar(ref.salida, par(suma, "ent", 0, 1))
    m.conectar(g.salida, par(suma, "ent", 1, 1))
    m.conectar(suma.salida, tf.entrada)
    m.conectar(tf.salida, g.entrada)
    return m, tf


def test_lazo_tf_feedthrough_detectado():
    # num/den de igual grado => b0 != 0 => lazo algebraico oculto
    m, tf = _lazo_tf([1.0, 2.0], [1.0, 2.0])
    with pytest.raises(ValueError, match="feedthrough"):
        m.run(0.1, registrar=[tf])


def test_lazo_tf_sin_feedthrough_corre():
    # den de grado mayor => b0 == 0 => sin feedthrough => sin error
    m, tf = _lazo_tf([1.0], [1.0, 2.0])
    res = m.run(0.1, registrar=[tf])
    assert np.isfinite(np.asarray(res["tf"])).all()


def test_lazo_pid_sin_planta_detectado():
    # PID (Kp) realimentado directo: no hay estados que rompan el lazo
    m = Modelo(dt=1e-3)
    ref = m.add(FuenteConstante("ref", 1.0))
    suma = m.add(Suma("err", (1.0, -1.0)))
    pid = m.add(PID("pid", Kp=5.0, Ki=0.0))
    m.conectar(ref.salida, par(suma, "ent", 0, 1))
    m.conectar(pid.salida, par(suma, "ent", 1, 1))
    m.conectar(suma.salida, pid.entrada)
    with pytest.raises(ValueError, match="feedthrough"):
        m.run(0.1, registrar=[pid])


def test_lazo_pid_integral_puro_corre():
    # Ki solamente: sin feedthrough (u depende solo del estado)
    m = Modelo(dt=1e-3)
    ref = m.add(FuenteConstante("ref", 1.0))
    suma = m.add(Suma("err", (1.0, -1.0)))
    pid = m.add(PID("pid", Kp=0.0, Ki=5.0, Kd=0.0))
    m.conectar(ref.salida, par(suma, "ent", 0, 1))
    m.conectar(pid.salida, par(suma, "ent", 1, 1))
    m.conectar(suma.salida, pid.entrada)
    res = m.run(0.1, registrar=[pid])
    assert np.isfinite(np.asarray(res["pid"])).all()