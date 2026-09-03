"""Validación de bloques de control: Saturar, Relay y PulsoRectangular."""

import numpy as np

from bloques_crysi import (
    Modelo, FuenteRampa, FuenteSeno, Ganancia, Scope,
    Saturar, Relay, PulsoRectangular,
)
from .helpers import correr_c_y_numpy, flat_c


def _compara(modelo, t_fin, registrar, rtol=1e-9, atol=1e-12, metodo=0):
    modelo.metodo = metodo
    res_c, arr_np, n_steps = correr_c_y_numpy(modelo, t_fin, registrar)
    arr_c = flat_c(res_c, registrar, modelo, n_steps)
    assert arr_c.shape == arr_np.shape, (arr_c.shape, arr_np.shape)
    np.testing.assert_allclose(arr_c, arr_np, rtol=rtol, atol=atol)
    return arr_c, arr_np


# ---------------------------------------------------------- Saturar

def test_saturar_c_numpy():
    m = Modelo(dt=0.01)
    r = m.add(FuenteRampa("r", pendiente=1.0, t_inicio=0.0, offset=-0.5))
    s = m.add(Saturar("s", u_min=0.0, u_max=1.0))
    m.conectar(r.salida, s.entrada)
    _compara(m, 2.0, [s.salida])


def test_saturar_clampa():
    m = Modelo(dt=0.01)
    r = m.add(FuenteRampa("r", pendiente=1.0, t_inicio=0.0, offset=-0.5))
    s = m.add(Saturar("s", u_min=0.0, u_max=1.0))
    m.conectar(r.salida, s.entrada)
    arr, _ = _compara(m, 2.0, [s.salida])
    y = arr[0]
    assert y.min() >= -1e-12 and y.max() <= 1.0 + 1e-12
    assert y[0] == 0.0                 # entrada negativa -> u_min
    assert y[-1] == 1.0                # entrada > 1 -> u_max


def test_saturar_passthrough():
    m = Modelo(dt=0.01)
    sn = m.add(FuenteSeno("sn", amplitud=0.5, frecuencia=1.0))
    s = m.add(Saturar("s", u_min=-1.0, u_max=1.0))
    m.conectar(sn.salida, s.entrada)
    arr, _ = _compara(m, 1.0, [s.salida])
    np.testing.assert_allclose(arr[0],
                               0.5 * np.sin(2 * np.pi * np.arange(101) * 0.01),
                               atol=1e-6)


# ---------------------------------------------------------- Relay

def test_relay_c_numpy():
    m = Modelo(dt=1e-3)
    sn = m.add(FuenteSeno("sn", amplitud=1.0, frecuencia=1.0))
    ry = m.add(Relay("ry", umbral_on=0.2, umbral_off=-0.2))
    m.conectar(sn.salida, ry.entrada)
    _compara(m, 1.5, [ry.salida])


def test_relay_histeresis():
    # con histeresis el encendido ocurre en arcsin(0.2)/2pi y el apagado
    # en la bajada cruzando -0.2 (pi + arcsin(0.2))/2pi
    m = Modelo(dt=1e-4)
    sn = m.add(FuenteSeno("sn", amplitud=1.0, frecuencia=1.0))
    ry = m.add(Relay("ry", umbral_on=0.2, umbral_off=-0.2))
    m.conectar(sn.salida, ry.entrada)
    arr, _ = _compara(m, 1.2, [ry.salida])
    y = arr[0]
    assert set(np.unique(y)) <= {0.0, 1.0}
    t = np.arange(len(y)) * 1e-4
    on = np.where(np.diff(y.astype(int)) == 1)[0] * 1e-4
    off = np.where(np.diff(y.astype(int)) == -1)[0] * 1e-4
    t_on = np.arcsin(0.2) / (2 * np.pi)        # 0.0322...
    t_off = (np.pi + np.arcsin(0.2)) / (2 * np.pi)  # 0.5322...
    assert abs(on[0] - t_on) < 5e-3
    assert abs(off[0] - t_off) < 5e-3
    # entre umbrales el estado se mantiene (no vuelve a 0 en [on,off])
    assert y[int(0.1 / 1e-4)] == 1.0


def test_relay_mantiene_estado_sin_cruce():
    # entrada constante en la zona muerta: nunca conmuta
    m = Modelo(dt=0.01)
    m2 = m
    sn = m2.add(FuenteSeno("sn", amplitud=0.1, frecuencia=0.1))
    ry = m2.add(Relay("ry", umbral_on=0.2, umbral_off=-0.2))
    m2.conectar(sn.salida, ry.entrada)
    arr, _ = _compara(m2, 2.0, [ry.salida])
    assert (arr[0] == 0.0).all()


# ---------------------------------------------------------- Pulso

def test_pulso_c_numpy():
    m = Modelo(dt=0.01)
    p = m.add(PulsoRectangular("p", amplitud=2.0, periodo=1.0, duty=0.25))
    _compara(m, 2.0, [p.salida])


def test_pulso_duty_y_media():
    m = Modelo(dt=1e-3)
    p = m.add(PulsoRectangular("p", amplitud=2.0, periodo=1.0, duty=0.25,
                               offset=0.5))
    arr, _ = _compara(m, 4.0, [p.salida])
    y = arr[0]
    assert set(np.unique(y)) <= {0.5, 2.5}
    # media sobre periodos enteros: offset + amp*duty
    assert abs(y.mean() - (0.5 + 2.0 * 0.25)) < 1e-3
    # los bordes estan en t = duty (primer periodo empieza en t=0)
    t = np.arange(len(y)) * 1e-3
    flancos = np.where(np.diff(y.astype(float)) != 0)[0] * 1e-3
    np.testing.assert_allclose(flancos, [0.25, 1.0, 1.25, 2.0, 2.25, 3.0,
                                         3.25], atol=2e-3)


def test_pulso_fase():
    # fase 0.5 con duty 0.5: el pulso se corre medio periodo -> arranca off
    m = Modelo(dt=1e-3)
    p = m.add(PulsoRectangular("p", amplitud=1.0, periodo=1.0, duty=0.5,
                               fase=0.5))
    arr, _ = _compara(m, 1.0, [p.salida])
    assert arr[0][0] == 0.0
    assert abs(arr[0].mean() - 0.5) < 1e-3


def test_pulso_validaciones():
    import pytest
    with pytest.raises(ValueError):
        PulsoRectangular("p", periodo=0.0)
    with pytest.raises(ValueError):
        PulsoRectangular("p", duty=0.0)
    with pytest.raises(ValueError):
        PulsoRectangular("p", duty=1.5)
    with pytest.raises(ValueError):
        Saturar("s", u_min=2.0, u_max=1.0)
    with pytest.raises(ValueError):
        Relay("r", umbral_on=-1.0, umbral_off=1.0)


# ---------------------------------------------------------- Display

def test_display_imprime_valor_final(capsys):
    from bloques_crysi import Display
    m = Modelo(dt=0.01)
    sn = m.add(FuenteSeno("sn", amplitud=1.0, frecuencia=0.25))
    d = m.add(Display("d", formato="%.6f"))
    m.conectar(sn.salida, d.entrada)
    m.run(1.0, registrar=[])
    out = capsys.readouterr().out
    assert "d:" in out
    esperado = float(np.sin(2 * np.pi * 0.25 * 1.0))
    assert abs(float(out.split(":")[-1]) - esperado) < 1e-4


def test_display_c_numpy():
    from bloques_crysi import Display
    m = Modelo(dt=0.01)
    sn = m.add(FuenteSeno("sn", amplitud=1.0, frecuencia=1.0))
    d = m.add(Display("d"))
    m.conectar(sn.salida, d.entrada)
    arr, _ = _compara(m, 1.0, [d.salida])
    np.testing.assert_allclose(arr[0],
                               np.sin(2 * np.pi * np.arange(101) * 0.01),
                               atol=1e-6)


# ---------------------------------------------------------- guardar_csv

def test_guardar_csv_1d(tmp_path):
    from bloques_crysi import FuenteConstante, Scope
    m = Modelo(dt=0.1)
    f = m.add(FuenteConstante("f", 2.5))
    sc = m.add(Scope("sc", mostrar=False))
    m.conectar(f.salida, sc.canales[0])
    res = m.run(0.3, registrar=[sc])
    ruta = tmp_path / "res.csv"
    res.guardar_csv(str(ruta))
    data = np.genfromtxt(str(ruta), delimiter=",", names=True)
    assert data.dtype.names == ("t", "sc")
    np.testing.assert_allclose(data["t"], np.arange(4) * 0.1)
    np.testing.assert_allclose(data["sc"], np.full(4, 2.5))


def test_guardar_csv_2d(tmp_path):
    m = Modelo(dt=0.01)
    sn = m.add(FuenteSeno("sn", amplitud=1.0, frecuencia=1.0))
    sc = m.add(Scope("sc", mostrar=False, anchos=[1, 1],
                     guiones=["seno", "seno_x2"]))
    g = m.add(Ganancia("g", 2.0))
    m.conectar(sn.salida, sc.canales[0])
    m.conectar(sn.salida, g.entrada)
    m.conectar(g.salida, sc.canales[1])
    res = m.run(1.0, registrar=[sc])
    ruta = tmp_path / "res2.csv"
    res.guardar_csv(str(ruta))
    data = np.genfromtxt(str(ruta), delimiter=",", names=True)
    # genfromtxt sanea "sc[0]" -> "sc0"
    assert data.dtype.names == ("t", "sc0", "sc1")
    np.testing.assert_allclose(data["sc1"], 2 * data["sc0"], atol=1e-9)