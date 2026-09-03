"""Entorno: FuenteCSV (interpolacion de puntos desde archivo) e inyeccion
de fallos (programado por tiempo y por evento).
"""

import numpy as np
import pytest

from bloques_crysi import (
    FalloEvento, FalloProgramado, FuenteConstante, FuenteCSV,
    FuenteEscalon, FuenteRampa, Modelo, PulsoRectangular,
)
from .helpers import correr_c_y_numpy, flat_c


def _compara_c_numpy(m, t_fin, registrar, atol=1e-12):
    res_c, arr_np, n_steps = correr_c_y_numpy(m, t_fin, registrar)
    flat = flat_c(res_c, registrar, m, n_steps)
    np.testing.assert_allclose(flat, arr_np, rtol=0.0, atol=atol)
    return res_c


def _escribir_csv(tmp_path, texto, nombre="rampa.csv"):
    archivo = tmp_path / nombre
    archivo.write_text(texto, encoding="utf-8")
    return str(archivo)


# ---------------------------------------------------------- FuenteCSV
def test_fuente_csv_interpola_linealmente(tmp_path):
    # rampa y = 10t
    archivo = _escribir_csv(
        tmp_path, "t,y\n0,0\n1,10\n2,20\n", nombre="rampa.csv")
    m = Modelo(dt=0.01)
    csv = m.add(FuenteCSV("csv", archivo))
    res = _compara_c_numpy(m, 2.0, [csv.salida])
    y = np.asarray(res["csv"]).ravel()
    t = np.arange(len(y)) * m.dt
    np.testing.assert_allclose(y, 10.0 * t, rtol=1e-12, atol=1e-12)


def test_fuente_csv_satura_fuera_de_rango(tmp_path):
    archivo = _escribir_csv(
        tmp_path, "0;0\n1;10\n2;20\n", nombre="sat.csv")
    m = Modelo(dt=0.01)
    csv = m.add(FuenteCSV("csv", archivo))
    res = _compara_c_numpy(m, 3.0, [csv.salida])
    y = np.asarray(res["csv"]).ravel()
    t = np.arange(len(y)) * m.dt
    esperado = np.clip(10.0 * t, 0.0, 20.0)
    np.testing.assert_allclose(y, esperado, rtol=1e-12, atol=1e-12)


def test_fuente_csv_sin_interpolar_mantiene_valor(tmp_path):
    # escalon: mantiene el ultimo valor hasta el proximo punto
    archivo = _escribir_csv(
        tmp_path, "0,0\n1,10\n2,20\n", nombre="hold.csv")
    m = Modelo(dt=0.01)
    csv = m.add(FuenteCSV("csv", archivo, interpolar=False))
    res = _compara_c_numpy(m, 2.0, [csv.salida])
    y = np.asarray(res["csv"]).ravel()
    t = np.arange(len(y)) * m.dt
    esperado = np.where(t < 1.0, 0.0, np.where(t < 2.0, 10.0, 20.0))
    np.testing.assert_allclose(y, esperado, rtol=1e-12, atol=1e-12)


def test_fuente_csv_ignora_encabezados_y_filas_raras(tmp_path):
    archivo = _escribir_csv(
        tmp_path, "tiempo,senal\nhola,mundo\n0,0\n1,5\n2,10\n",
        nombre="raro.csv")
    m = Modelo(dt=0.01)
    csv = m.add(FuenteCSV("csv", archivo))
    res = _compara_c_numpy(m, 2.0, [csv.salida])
    y = np.asarray(res["csv"]).ravel()
    t = np.arange(len(y)) * m.dt
    np.testing.assert_allclose(y, 5.0 * t, rtol=1e-12, atol=1e-12)


def test_fuente_csv_validaciones(tmp_path):
    with pytest.raises(FileNotFoundError):
        FuenteCSV("csv", str(tmp_path / "no_existe.csv"))
    archivo = _escribir_csv(tmp_path, "0,0\n", nombre="unpunto.csv")
    with pytest.raises(ValueError):
        FuenteCSV("csv", archivo)
    archivo = _escribir_csv(tmp_path, "0,0\n0,5\n1,10\n", nombre="repetido.csv")
    with pytest.raises(ValueError):
        FuenteCSV("csv", archivo)


# ---------------------------------------------------------- fallos
def test_fallo_programado_reemplaza_desde_instante():
    m = Modelo(dt=0.01)
    u = m.add(FuenteEscalon("u", valor_final=5.0, t_paso=0.0))
    fallo = m.add(FalloProgramado("fallo", t_fallo=1.0, valor=0.0))
    m.conectar(u.salida, fallo.entrada)
    res = _compara_c_numpy(m, 2.0, [fallo.salida])
    y = np.asarray(res["fallo"]).ravel()
    t = np.arange(len(y)) * m.dt
    esperado = np.where(t < 1.0, 5.0, 0.0)
    np.testing.assert_allclose(y, esperado, rtol=1e-12, atol=1e-12)


def test_fallo_programado_modo_suma():
    m = Modelo(dt=0.01)
    u = m.add(FuenteRampa("u", pendiente=2.0))
    fallo = m.add(FalloProgramado("fallo", t_fallo=1.0, valor=-10.0, modo=1))
    m.conectar(u.salida, fallo.entrada)
    res = _compara_c_numpy(m, 2.0, [fallo.salida])
    y = np.asarray(res["fallo"]).ravel()
    t = np.arange(len(y)) * m.dt
    esperado = np.where(t < 1.0, 2.0 * t, 2.0 * t - 10.0)
    np.testing.assert_allclose(y, esperado, rtol=1e-12, atol=1e-12)


def test_fallo_evento_dispara_por_umbral():
    # disparo alto en [0,0.5) (pulso 0..1 periodo 1 duty 0.5), bajo en [0.5,1)
    m = Modelo(dt=0.01)
    u = m.add(FuenteEscalon("u", valor_final=7.0, t_paso=0.0))
    disp = m.add(PulsoRectangular("disp", amplitud=1.0, periodo=1.0,
                                  duty=0.5, offset=0.0))
    fallo = m.add(FalloEvento("fallo", umbral=0.5, valor=-3.0))
    m.conectar(u.salida, fallo.senal)
    m.conectar(disp.salida, fallo.disparo)
    res = _compara_c_numpy(m, 1.0, [fallo.salida])
    y = np.asarray(res["fallo"]).ravel()
    t = np.arange(len(y)) * m.dt
    esperado = np.where((t % 1.0) < 0.5, -3.0, 7.0)
    np.testing.assert_allclose(y, esperado, rtol=1e-12, atol=1e-12)


def test_fallo_evento_no_dispara_sin_trigger():
    m = Modelo(dt=0.01)
    u = m.add(FuenteEscalon("u", valor_final=7.0, t_paso=0.0))
    disp = m.add(FuenteConstante("disp", 0.0))
    fallo = m.add(FalloEvento("fallo", umbral=0.5, valor=-3.0))
    m.conectar(u.salida, fallo.senal)
    m.conectar(disp.salida, fallo.disparo)
    res = _compara_c_numpy(m, 1.0, [fallo.salida])
    y = np.asarray(res["fallo"]).ravel()
    np.testing.assert_allclose(y, 7.0, rtol=1e-12, atol=1e-12)


def test_fallo_validaciones():
    with pytest.raises(ValueError):
        FalloProgramado("f", t_fallo=-1.0, valor=0.0)
    with pytest.raises(ValueError):
        FalloProgramado("f", t_fallo=0.0, valor=0.0, modo=2)
    with pytest.raises(ValueError):
        FalloEvento("f", umbral=0.5, valor=0.0, modo=3)