"""Tejido conectivo: mux/demux, tablas 1D/2D/3D, logicos y relacionales.

Cada bloque se valida contra valores analiticos y contra el backend
numpy de referencia (los dos nucleos deben coincidir bit a bit).
"""

import numpy as np
import pytest

from bloques_crysi import (
    Demultiplexor, FuenteConstante, FuenteRampa, Logico, Modelo,
    Multiplexor, Relacional, Tabla1D, Tabla2D, Tabla3D,
)
from .helpers import correr_c_y_numpy, flat_c


def _compara_c_numpy(m, t_fin, registrar, atol=1e-12):
    res_c, arr_np, n_steps = correr_c_y_numpy(m, t_fin, registrar)
    flat = flat_c(res_c, registrar, m, n_steps)
    np.testing.assert_allclose(flat, arr_np, rtol=0.0, atol=atol)
    return res_c


# ---------------------------------------------------------------- mux/demux
def test_mux_demux_roundtrip():
    m = Modelo(dt=1e-3)
    v = [m.add(FuenteConstante(f"v{k}", 1.5 * k - 2.0)) for k in range(3)]
    mux = m.add(Multiplexor("mux", n_canales=3))
    demux = m.add(Demultiplexor("demux", n_canales=3))
    for k in range(3):
        m.conectar(v[k].salida, mux.entradas[k])
    m.conectar(mux.salida, demux.entrada)
    res = _compara_c_numpy(m, 0.1, [demux])
    esperado = np.array([1.5 * k - 2.0 for k in range(3)])
    np.testing.assert_allclose(res["demux"][0], esperado,
                               rtol=0.0, atol=1e-12)
    # todas las muestras constantes
    np.testing.assert_allclose(np.diff(res["demux"], axis=0), 0.0,
                               rtol=0.0, atol=1e-12)


def test_demux_reparte_bus_por_canal():
    m = Modelo(dt=1e-3)
    v = [m.add(FuenteConstante(f"v{k}", 10.0 * k + 1.0)) for k in range(3)]
    mux = m.add(Multiplexor("mux", n_canales=3))
    demux = m.add(Demultiplexor("demux", n_canales=3))
    for k in range(3):
        m.conectar(v[k].salida, mux.entradas[k])
    m.conectar(mux.salida, demux.entrada)
    res = _compara_c_numpy(m, 0.1, [demux.salidas[1]])
    # el canal 1 del bus (11.0) llega intacto a la salida 1 del demux
    np.testing.assert_allclose(np.asarray(res["demux"]).ravel(), 11.0,
                               rtol=0.0, atol=1e-12)


def test_mux_n_canales_invalidos():
    with pytest.raises(ValueError, match="n_canales"):
        Multiplexor("mux", 0)
    with pytest.raises(ValueError, match="n_canales"):
        Demultiplexor("demux", 0)


# ---------------------------------------------------------------- Tabla1D
def test_tabla1d_interpola_lineal():
    m = Modelo(dt=1e-3)
    rampa = m.add(FuenteRampa("u", pendiente=1.0))
    tabla = m.add(Tabla1D("t1", [0.0, 1.0, 2.0, 3.0], [10.0, 20.0, 40.0, 80.0]))
    m.conectar(rampa.salida, tabla.entrada)
    res = _compara_c_numpy(m, 3.0, [tabla.salida])
    y = np.asarray(res["t1"]).ravel()
    t = np.arange(len(y)) * m.dt
    esperado = np.where(t <= 1.0, 10.0 + 10.0 * t,
                        np.where(t <= 2.0, 20.0 + 20.0 * (t - 1.0),
                                 40.0 + 40.0 * (t - 2.0)))
    np.testing.assert_allclose(y, esperado, rtol=1e-12, atol=1e-12)


def test_tabla1d_satura_fuera_de_rango():
    m = Modelo(dt=1e-3)
    c = m.add(FuenteConstante("u", 10.0))
    tabla = m.add(Tabla1D("t1", [0.0, 1.0, 2.0], [5.0, 6.0, 7.0]))
    m.conectar(c.salida, tabla.entrada)
    res = _compara_c_numpy(m, 0.05, [tabla.salida])
    np.testing.assert_allclose(res["t1"], 7.0, rtol=0.0, atol=1e-12)


def test_tabla1d_validacion():
    with pytest.raises(ValueError, match="igual longitud"):
        Tabla1D("t", [0.0, 1.0], [1.0])
    with pytest.raises(ValueError, match="2 puntos"):
        Tabla1D("t", [0.0], [1.0])
    with pytest.raises(ValueError, match="creciente"):
        Tabla1D("t", [0.0, 1.0, 0.5], [1.0, 2.0, 3.0])


# ---------------------------------------------------------------- Tabla2D
def test_tabla2d_bilineal():
    z = [[0.0, 10.0, 20.0],     # y = 0
         [100.0, 110.0, 120.0]]  # y = 1
    m = Modelo(dt=1e-3)
    u1 = m.add(FuenteConstante("u1", 0.5))
    u2 = m.add(FuenteConstante("u2", 0.5))
    tabla = m.add(Tabla2D("t2", [0.0, 1.0, 2.0], [0.0, 1.0], z))
    m.conectar(u1.salida, tabla.entrada1)
    m.conectar(u2.salida, tabla.entrada2)
    res = _compara_c_numpy(m, 0.05, [tabla.salida])
    # u1=0.5 -> 5 en y=0 y 105 en y=1; u2=0.5 -> 55
    np.testing.assert_allclose(res["t2"], 55.0, rtol=0.0, atol=1e-12)


def test_tabla2d_en_punto_de_grilla():
    z = [[0.0, 10.0, 20.0], [100.0, 110.0, 120.0]]
    m = Modelo(dt=1e-3)
    u1 = m.add(FuenteConstante("u1", 1.0))
    u2 = m.add(FuenteConstante("u2", 1.0))
    tabla = m.add(Tabla2D("t2", [0.0, 1.0, 2.0], [0.0, 1.0], z))
    m.conectar(u1.salida, tabla.entrada1)
    m.conectar(u2.salida, tabla.entrada2)
    res = _compara_c_numpy(m, 0.05, [tabla.salida])
    np.testing.assert_allclose(res["t2"], 110.0, rtol=0.0, atol=1e-12)


def test_tabla2d_satura_fuera_de_rango():
    z = [[0.0, 10.0, 20.0], [100.0, 110.0, 120.0]]
    m = Modelo(dt=1e-3)
    u1 = m.add(FuenteConstante("u1", -3.0))
    u2 = m.add(FuenteConstante("u2", 5.0))
    tabla = m.add(Tabla2D("t2", [0.0, 1.0, 2.0], [0.0, 1.0], z))
    m.conectar(u1.salida, tabla.entrada1)
    m.conectar(u2.salida, tabla.entrada2)
    res = _compara_c_numpy(m, 0.05, [tabla.salida])
    # u1=-3 -> col x=0; u2=5 -> fila y=1 -> z[1][0] = 100
    np.testing.assert_allclose(res["t2"], 100.0, rtol=0.0, atol=1e-12)


def test_tabla2d_validacion():
    with pytest.raises(ValueError, match="2 puntos"):
        Tabla2D("t", [0.0], [0.0, 1.0], [[1.0], [2.0]])
    with pytest.raises(ValueError, match="creciente"):
        Tabla2D("t", [0.0, 0.5, 0.4], [0.0, 1.0], [[0.0, 1.0, 2.0],
                                                   [3.0, 4.0, 5.0]])
    with pytest.raises(ValueError, match="tabla debe tener"):
        Tabla2D("t", [0.0, 1.0], [0.0, 1.0], [[0.0, 1.0, 2.0]])


# ---------------------------------------------------------------- Tabla3D
def test_tabla3d_trilineal():
    z = [[[0.0, 10.0], [20.0, 30.0]],    # k=0: filas y=0,1; cols x=0,1
         [[30.0, 40.0], [50.0, 60.0]]]   # k=1
    m = Modelo(dt=1e-3)
    u1 = m.add(FuenteConstante("u1", 0.5))
    u2 = m.add(FuenteConstante("u2", 0.5))
    u3 = m.add(FuenteConstante("u3", 0.5))
    tabla = m.add(Tabla3D("t3", [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], z))
    m.conectar(u1.salida, tabla.entrada1)
    m.conectar(u2.salida, tabla.entrada2)
    m.conectar(u3.salida, tabla.entrada3)
    res = _compara_c_numpy(m, 0.05, [tabla.salida])
    # promedio de las 8 esquinas: (0+10+20+30+30+40+50+60)/8 = 30
    np.testing.assert_allclose(res["t3"], 30.0, rtol=0.0, atol=1e-12)


def test_tabla3d_punto_interior():
    z = [[[0.0, 10.0], [20.0, 30.0]],
         [[30.0, 40.0], [50.0, 60.0]]]
    m = Modelo(dt=1e-3)
    u1 = m.add(FuenteConstante("u1", 0.25))
    u2 = m.add(FuenteConstante("u2", 1.0))   # tope de y: fila y=1
    u3 = m.add(FuenteConstante("u3", 0.5))
    tabla = m.add(Tabla3D("t3", [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], z))
    m.conectar(u1.salida, tabla.entrada1)
    m.conectar(u2.salida, tabla.entrada2)
    m.conectar(u3.salida, tabla.entrada3)
    res = _compara_c_numpy(m, 0.05, [tabla.salida])
    # en y=1: u1=0.25 -> 22.5 (k=0) y 52.5 (k=1); promedio en z -> 37.5
    np.testing.assert_allclose(res["t3"], 37.5, rtol=0.0, atol=1e-12)


def test_tabla3d_validacion():
    with pytest.raises(ValueError, match="creciente"):
        Tabla3D("t", [0.0, 0.0], [0.0, 1.0], [0.0, 1.0],
                [[[0.0, 0.0], [0.0, 0.0]], [[0.0, 0.0], [0.0, 0.0]]])
    with pytest.raises(ValueError, match="tabla debe tener"):
        Tabla3D("t", [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [[[0.0]]])


# ---------------------------------------------------------------- Logico
def _arma_logico(opcion, entradas, n_entradas=None, umbral=0.5):
    m = Modelo(dt=1e-3)
    bl = m.add(Logico("l", opcion=opcion,
                      n_entradas=n_entradas if n_entradas else len(entradas),
                      umbral=umbral))
    for k, val in enumerate(entradas):
        c = m.add(FuenteConstante(f"c{k}", val))
        m.conectar(c.salida, bl.entradas[k])
    return m, bl


def test_logico_tabla_de_verdad():
    casos = [("AND", (1.0, 1.0), 1.0), ("AND", (1.0, 0.0), 0.0),
             ("OR", (0.0, 0.0), 0.0), ("OR", (0.0, 1.0), 1.0),
             ("NAND", (1.0, 1.0), 0.0), ("NOR", (0.0, 0.0), 1.0),
             ("XOR", (1.0, 1.0), 0.0), ("XOR", (1.0, 0.0), 1.0),
             ("XNOR", (1.0, 0.0), 0.0), ("NOT", (0.0,), 1.0)]
    for opcion, entradas, esperado in casos:
        m, bl = _arma_logico(opcion, entradas)
        res = _compara_c_numpy(m, 0.05, [bl.salida])
        np.testing.assert_allclose(res["l"], esperado, rtol=0.0, atol=1e-12,
                                   err_msg=opcion)


def test_logico_multi_entrada_y_umbral():
    m, bl = _arma_logico("AND", (0.6, 0.7, 0.9), n_entradas=3)
    res = _compara_c_numpy(m, 0.05, [bl.salida])
    np.testing.assert_allclose(res["l"], 1.0, rtol=0.0, atol=1e-12)
    m2, bl2 = _arma_logico("AND", (0.4, 0.6, 0.9), n_entradas=3, umbral=0.7)
    res2 = _compara_c_numpy(m2, 0.05, [bl2.salida])
    np.testing.assert_allclose(res2["l"], 0.0, rtol=0.0, atol=1e-12)


def test_logico_validacion():
    with pytest.raises(ValueError, match="opcion"):
        Logico("l", opcion="SUMA")
    with pytest.raises(ValueError, match="n_entradas"):
        Logico("l", opcion="AND", n_entradas=0)


# ---------------------------------------------------------------- Relacional
def _arma_relacional(opcion, a, b, tol=0.0):
    m = Modelo(dt=1e-3)
    bl = m.add(Relacional("r", opcion=opcion, tol=tol))
    ca = m.add(FuenteConstante("a", a))
    cb = m.add(FuenteConstante("b", b))
    m.conectar(ca.salida, bl.a)
    m.conectar(cb.salida, bl.b)
    return m, bl


def test_relacional_tabla_de_verdad():
    casos = [("<", (2.0, 3.0), 1.0), ("<", (3.0, 2.0), 0.0),
             ("<=", (3.0, 3.0), 1.0), (">", (3.0, 2.0), 1.0),
             (">=", (2.0, 3.0), 0.0), ("==", (2.0, 2.0), 1.0),
             ("==", (2.0, 3.0), 0.0), ("!=", (2.0, 3.0), 1.0)]
    for opcion, (a, b), esperado in casos:
        m, bl = _arma_relacional(opcion, a, b)
        res = _compara_c_numpy(m, 0.05, [bl.salida])
        np.testing.assert_allclose(res["r"], esperado, rtol=0.0, atol=1e-12,
                                   err_msg=opcion)


def test_relacional_tolerancia():
    m, bl = _arma_relacional("==", 2.0, 2.0 + 1e-6, tol=1e-5)
    res = _compara_c_numpy(m, 0.05, [bl.salida])
    np.testing.assert_allclose(res["r"], 1.0, rtol=0.0, atol=1e-12)
    m2, bl2 = _arma_relacional("!=", 2.0, 2.0 + 1e-6, tol=1e-9)
    res2 = _compara_c_numpy(m2, 0.05, [bl2.salida])
    np.testing.assert_allclose(res2["r"], 1.0, rtol=0.0, atol=1e-12)


def test_relacional_validacion():
    with pytest.raises(ValueError, match="opcion"):
        Relacional("r", opcion="~")
    with pytest.raises(ValueError, match="tol"):
        Relacional("r", tol=-1.0)


def test_relacional_con_rampa_dinamica():
    # u1 = t, u2 = 1.5 fijo: '>=' emite 1 una vez que t alcanza 1.5
    m = Modelo(dt=0.01)
    rampa = m.add(FuenteRampa("u", pendiente=1.0))
    c = m.add(FuenteConstante("ref", 1.5))
    rel = m.add(Relacional("r", opcion=">="))
    m.conectar(rampa.salida, rel.a)
    m.conectar(c.salida, rel.b)
    res = _compara_c_numpy(m, 3.0, [rel.salida])
    y = np.asarray(res["r"]).ravel()
    t = np.arange(len(y)) * m.dt
    np.testing.assert_allclose(y, (t >= 1.5).astype(float),
                               rtol=0.0, atol=1e-12)


# ----------------------------------------------------------------- multiplicador
def test_multiplicador_producto():
    from bloques_crysi import Multiplicador

    m = Modelo(dt=1e-3)
    a = m.add(FuenteConstante("a", 2.5))
    b = m.add(FuenteRampa("b", pendiente=3.0))      # b(t) = 3*t
    mux = m.add(Multiplexor("mux", n_canales=2))
    mult = m.add(Multiplicador("mult"))
    m.conectar(a.salida, mux.entradas[0])
    m.conectar(b.salida, mux.entradas[1])
    m.conectar(mux.salida, mult.entrada)
    res = _compara_c_numpy(m, 0.4, [mult.salida])
    y = np.asarray(res["mult"]).ravel()
    t = np.arange(len(y)) * m.dt
    np.testing.assert_allclose(y, 2.5 * 3.0 * t, rtol=0.0, atol=1e-12)