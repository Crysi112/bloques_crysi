"""Fisica multidominio: bloques termicos (MasaTermica, ResistenciaTermica)
y de transmision mecanica (Engranaje, EjeFlexible, Embrague).
Validacion analitica + C vs numpy.
"""

import numpy as np
import pytest

from bloques_crysi import (
    EjeFlexible, Embrague, Engranaje, FuenteConstante, FuenteEscalon,
    FuenteRampa, MasaTermica, Modelo, ResistenciaTermica,
)
from .helpers import correr_c_y_numpy, flat_c


def _compara_c_numpy(m, t_fin, registrar, atol=1e-12):
    res_c, arr_np, n_steps = correr_c_y_numpy(m, t_fin, registrar)
    flat = flat_c(res_c, registrar, m, n_steps)
    np.testing.assert_allclose(flat, arr_np, rtol=0.0, atol=atol)
    return res_c


# ---------------------------------------------------------- termicos
def test_masa_termica_integra_flujo_constante():
    # Q = 100 W, C = 50 J/C -> dT/dt = 2 C/s (Euler exacto para constante)
    m = Modelo(dt=0.01)
    q = m.add(FuenteConstante("q", 100.0))
    masa = m.add(MasaTermica("masa", C_th=50.0, T_inicial=10.0))
    m.conectar(q.salida, masa.entradas[0])
    res = _compara_c_numpy(m, 2.0, [masa.salida])
    y = np.asarray(res["masa"]).ravel()
    t = np.arange(len(y)) * m.dt
    np.testing.assert_allclose(y, 10.0 + 2.0 * t, rtol=1e-12, atol=1e-12)


def test_masa_termica_con_ambiente():
    # Q=0, T_amb=25, R_amb=1, C=10: T' = -(T-25)/10
    # Euler: T_k = 25 + (T0-25)*(1 - dt/(R*C))^k (cerrado exacto)
    m = Modelo(dt=0.01)
    q = m.add(FuenteConstante("q", 0.0))
    masa = m.add(MasaTermica("masa", C_th=10.0, T_inicial=30.0,
                             T_amb=25.0, R_amb=1.0))
    m.conectar(q.salida, masa.entradas[0])
    res = _compara_c_numpy(m, 20.0, [masa.salida])
    y = np.asarray(res["masa"]).ravel()
    k = np.arange(len(y))
    esperado = 25.0 + 5.0 * (1.0 - m.dt / 10.0) ** k
    np.testing.assert_allclose(y, esperado, rtol=1e-12, atol=1e-12)
    assert y[-1] < 27.0  # descendiendo hacia el ambiente


def test_resistencia_termica_flujo():
    m = Modelo(dt=0.01)
    t1 = m.add(FuenteConstante("t1", 50.0))
    t2 = m.add(FuenteConstante("t2", 25.0))
    r = m.add(ResistenciaTermica("r", R=10.0))
    m.conectar(t1.salida, r.T1)
    m.conectar(t2.salida, r.T2)
    res = _compara_c_numpy(m, 0.1, [r.salida])
    y = np.asarray(res["r"]).ravel()
    np.testing.assert_allclose(y, 2.5, rtol=1e-12, atol=1e-12)


def test_circuito_termico_masa_resistencia():
    # T1=100 (escalon) -> R=10 -> masa C=100 (T0=25):
    # Q = (100 - T)/10 ; dT/dt = (100 - T)/1000
    # Euler: T_k = 100 - 75*(1 - dt/1000)^k (cerrado exacto)
    m = Modelo(dt=0.5)
    src = m.add(FuenteEscalon("src", valor_final=100.0, t_paso=0.0))
    r = m.add(ResistenciaTermica("r", R=10.0))
    masa = m.add(MasaTermica("masa", C_th=100.0, T_inicial=25.0))
    m.conectar(src.salida, r.T1)
    m.conectar(masa.salida, r.T2)
    m.conectar(r.salida, masa.entradas[0])
    res = _compara_c_numpy(m, 3000.0, [masa.salida])
    y = np.asarray(res["masa"]).ravel()
    k = np.arange(len(y))
    esperado = 100.0 - 75.0 * (1.0 - m.dt / 1000.0) ** k
    np.testing.assert_allclose(y, esperado, rtol=1e-12, atol=1e-12)
    assert y[-1] > 95.0  # subiendo hacia el equilibrio


# ---------------------------------------------------------- transmision
def test_engranaje_conserva_potencia():
    # a=2: w2 = 2*w1 = 200 ; T2 = T1/2 = 2.5
    m = Modelo(dt=0.01)
    w1 = m.add(FuenteConstante("w1", 100.0))
    t1 = m.add(FuenteConstante("t1", 5.0))
    eng = m.add(Engranaje("eng", relacion=2.0))
    m.conectar(w1.salida, eng.w1)
    m.conectar(t1.salida, eng.T1)
    res = _compara_c_numpy(m, 0.1, [eng.salida])
    y = np.asarray(res["eng"])
    np.testing.assert_allclose(y[:, 0], 200.0, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(y[:, 1], 2.5, rtol=1e-12, atol=1e-12)
    # potencia conservada: T1*w1 == T2*w2
    np.testing.assert_allclose(y[:, 1] * y[:, 0], 500.0, rtol=1e-12, atol=1e-12)


def test_eje_flexible_torsion_lineal():
    # w1=10, w2=0, K=5, B=0: th1 = 10t -> T = 50t (Euler exacto)
    m = Modelo(dt=0.001)
    w1 = m.add(FuenteConstante("w1", 10.0))
    w2 = m.add(FuenteConstante("w2", 0.0))
    eje = m.add(EjeFlexible("eje", K=5.0))
    m.conectar(w1.salida, eje.w1)
    m.conectar(w2.salida, eje.w2)
    res = _compara_c_numpy(m, 1.0, [eje.salida])
    y = np.asarray(res["eje"]).ravel()
    t = np.arange(len(y)) * m.dt
    np.testing.assert_allclose(y, 50.0 * t, rtol=1e-12, atol=1e-12)


def test_eje_flexible_con_amortiguacion():
    # K=5, B=2: T = 50t + 2*10 = 50t + 20
    m = Modelo(dt=0.001)
    w1 = m.add(FuenteConstante("w1", 10.0))
    w2 = m.add(FuenteConstante("w2", 0.0))
    eje = m.add(EjeFlexible("eje", K=5.0, B=2.0))
    m.conectar(w1.salida, eje.w1)
    m.conectar(w2.salida, eje.w2)
    res = _compara_c_numpy(m, 1.0, [eje.salida])
    y = np.asarray(res["eje"]).ravel()
    t = np.arange(len(y)) * m.dt
    np.testing.assert_allclose(y, 50.0 * t + 20.0, rtol=1e-12, atol=1e-12)


def test_embrague_saturacion_de_par():
    m = Modelo(dt=0.01)
    t = m.add(FuenteConstante("t", 50.0))
    c = m.add(FuenteEscalon("c", valor_final=1.0, t_paso=0.0))
    emb = m.add(Embrague("emb", T_max=20.0))
    m.conectar(t.salida, emb.entrada)
    m.conectar(c.salida, emb.control)
    res = _compara_c_numpy(m, 0.1, [emb.salida])
    y = np.asarray(res["emb"]).ravel()
    np.testing.assert_allclose(y, 20.0, rtol=1e-12, atol=1e-12)


def test_embrague_desembragado_no_transmite():
    m = Modelo(dt=0.01)
    t = m.add(FuenteConstante("t", 50.0))
    c = m.add(FuenteConstante("c", 0.0))
    emb = m.add(Embrague("emb", T_max=20.0))
    m.conectar(t.salida, emb.entrada)
    m.conectar(c.salida, emb.control)
    res = _compara_c_numpy(m, 0.1, [emb.salida])
    y = np.asarray(res["emb"]).ravel()
    np.testing.assert_allclose(y, 0.0, rtol=0.0, atol=1e-12)


def test_embrague_deja_pasar_par_menor():
    m = Modelo(dt=0.01)
    t = m.add(FuenteConstante("t", 10.0))
    c = m.add(FuenteEscalon("c", valor_final=1.0, t_paso=0.0))
    emb = m.add(Embrague("emb", T_max=20.0))
    m.conectar(t.salida, emb.entrada)
    m.conectar(c.salida, emb.control)
    res = _compara_c_numpy(m, 0.1, [emb.salida])
    y = np.asarray(res["emb"]).ravel()
    np.testing.assert_allclose(y, 10.0, rtol=1e-12, atol=1e-12)


def test_embrague_rampa_entrada_y_control():
    # par de entrada sube 0..100, control=1 desde t=0: y = min(t*100, 30)
    m = Modelo(dt=0.001)
    t = m.add(FuenteRampa("t", pendiente=100.0))
    c = m.add(FuenteEscalon("c", valor_final=1.0, t_paso=0.0))
    emb = m.add(Embrague("emb", T_max=30.0))
    m.conectar(t.salida, emb.entrada)
    m.conectar(c.salida, emb.control)
    res = _compara_c_numpy(m, 1.0, [emb.salida])
    y = np.asarray(res["emb"]).ravel()
    tvec = np.arange(len(y)) * m.dt
    np.testing.assert_allclose(y, np.minimum(tvec * 100.0, 30.0),
                               rtol=1e-12, atol=1e-12)


def test_validaciones_fisica():
    with pytest.raises(ValueError):
        MasaTermica("m", C_th=0.0)
    with pytest.raises(ValueError):
        ResistenciaTermica("r", R=0.0)
    with pytest.raises(ValueError):
        Engranaje("e", relacion=0.0)
    with pytest.raises(ValueError):
        EjeFlexible("e", K=-1.0)
    with pytest.raises(ValueError):
        Embrague("emb", T_max=-1.0)