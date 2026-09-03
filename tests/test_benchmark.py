"""Benchmark C vs backend numpy: el nucleo C debe ser ~10-100x mas rapido."""

import time

import numpy as np

from bloques_crysi import (
    Modelo, FuenteConstante, FuenteTrifasica,
    MaquinaImanesPermanentes,
)
from .helpers import correr_c_y_numpy, resolver_registro
from bloques_crysi.backend_numpy import simular


def _arma_pmac_1ms():
    m = Modelo(dt=1e-4)
    src = m.add(FuenteTrifasica("red", 310.0, 50.0))
    tl = m.add(FuenteConstante("tl", 5.0))
    n = 4
    maqs = [m.add(MaquinaImanesPermanentes(
        f"pmsm{i}", rs=0.1, Ld=1e-3, Lq=1e-3, lam_m=0.1, P=6, J=0.01, Bm=0.001))
        for i in range(n)]
    for maq in maqs:
        m.conectar(src.salida, maq.terminales)
        m.conectar(tl.salida, maq.T_L)
    reg = [maq.sensorVelocidad() for maq in maqs]
    return m, reg


def test_c_mas_rapido_que_numpy():
    m, reg = _arma_pmac_1ms()
    m._resolver()
    rec_idx, _ = resolver_registro(m, reg)

    t_fin = 1.0
    n_steps = int(round(t_fin / m.dt)) + 1

    def corrida_c():
        t0 = time.perf_counter()
        m.run(t_fin, registrar=reg)
        return time.perf_counter() - t0

    def corrida_np():
        t0 = time.perf_counter()
        simular(m.bloques, m.dt, t_fin, rec_idx, metodo=0,
                max_iter=m.max_iter, tol=m.tol, w_opt=m.w_opt,
                orden_estatico=m._orden_estatico())
        return time.perf_counter() - t0

    T = 3
    for _ in range(2):  # warm-up
        corrida_c(); corrida_np()
    t_c = min(corrida_c() for _ in range(T))
    t_np = min(corrida_np() for _ in range(T))

    assert t_c < t_np, "C no debe ser mas lento que numpy"
    factor = t_np / max(t_c, 1e-9)
    assert factor > 3.0, f"aceleracion C vs numpy solo {factor:.1f}x (deseado >3x)"
    print(f"  C: {t_c*1e3:.2f} ms | numpy: {t_np*1e3:.2f} ms | {factor:.0f}x")