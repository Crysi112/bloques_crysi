"""Cerebro: limitador de rapidez, filtros configurables, retenedor
disparado y maquina de estados. Validacion analitica + C vs numpy.
"""

import numpy as np
import pytest

from bloques_crysi import (
    FiltroNotch, FiltroPasoAlto, FiltroPasoBajo, FuenteEscalon,
    FuenteSeno, LimitadorRapidez, MaquinaEstados, Modelo,
    PulsoRectangular, RetenedorDisparado,
)
from .helpers import correr_c_y_numpy, flat_c


def _compara_c_numpy(m, t_fin, registrar, atol=1e-12):
    res_c, arr_np, n_steps = correr_c_y_numpy(m, t_fin, registrar)
    flat = flat_c(res_c, registrar, m, n_steps)
    np.testing.assert_allclose(flat, arr_np, rtol=0.0, atol=atol)
    return res_c


# ---------------------------------------------------------- limitador
def test_limitador_rapidez_pendiente_saturada():
    # escalon 0 -> 10 con subida 2/s: y = min(2t, 10)
    m = Modelo(dt=1e-3)
    esc = m.add(FuenteEscalon("u", valor_final=10.0, t_paso=0.0))
    lim = m.add(LimitadorRapidez("lim", subida=2.0, bajada=5.0))
    m.conectar(esc.salida, lim.entrada)
    res = _compara_c_numpy(m, 6.0, [lim.salida])
    y = np.asarray(res["lim"]).ravel()
    t = np.arange(len(y)) * m.dt
    np.testing.assert_allclose(y, np.minimum(2.0 * t, 10.0),
                               rtol=1e-9, atol=1e-9)


def test_limitador_rapidez_bajada():
    # escalon 10 -> 0 con bajada 1/s: y = max(10 - t, 0)
    m = Modelo(dt=1e-3)
    esc = m.add(FuenteEscalon("u", valor_final=0.0, t_paso=0.0,
                              valor_inicial=10.0))
    lim = m.add(LimitadorRapidez("lim", subida=10.0, bajada=1.0,
                                 valor_inicial=10.0))
    m.conectar(esc.salida, lim.entrada)
    res = _compara_c_numpy(m, 12.0, [lim.salida])
    y = np.asarray(res["lim"]).ravel()
    t = np.arange(len(y)) * m.dt
    np.testing.assert_allclose(y, np.maximum(10.0 - t, 0.0),
                               rtol=1e-9, atol=1e-9)


def test_limitador_rapidez_segunda_etapa_sigue_a_la_primera():
    # rampa con pendiente 0.5/s: por debajo del limite 2/s, sigue sin saturar
    m = Modelo(dt=1e-3)
    esc = m.add(FuenteEscalon("u", valor_final=3.0, t_paso=0.0))
    lim = m.add(LimitadorRapidez("lim", subida=2.0, bajada=2.0))
    m.conectar(esc.salida, lim.entrada)
    res = _compara_c_numpy(m, 2.0, [lim.salida])
    y = np.asarray(res["lim"]).ravel()
    t = np.arange(len(y)) * m.dt
    np.testing.assert_allclose(y, np.minimum(2.0 * t, 3.0),
                               rtol=1e-9, atol=1e-9)


def test_limitador_rapidez_validacion():
    with pytest.raises(ValueError, match="subida"):
        LimitadorRapidez("lim", subida=-1.0)
    with pytest.raises(ValueError, match="bajada"):
        LimitadorRapidez("lim", bajada=-1.0)


# ---------------------------------------------------------- filtros
def test_filtro_pasobajo1_respuesta_escalon():
    # LPF1 fc=10 Hz: y(t) -> 1 - exp(-wc t). La TF del nucleo arranca
    # con historia de entrada en cero (un escalon de "medio paso"): la
    # comparacion analitica se hace una vez asentado el transitorio.
    m = Modelo(dt=1e-4)
    esc = m.add(FuenteEscalon("u", valor_final=1.0))
    f = m.add(FiltroPasoBajo("f", fc=10.0, orden=1))
    m.conectar(esc.salida, f.entrada)
    res = _compara_c_numpy(m, 0.3, [f.salida], atol=1e-8)
    y = np.asarray(res["f"]).ravel()
    t = np.arange(len(y)) * m.dt
    wc = 2.0 * np.pi * 10.0
    np.testing.assert_allclose(y[0], 0.0, rtol=0.0, atol=1e-12)
    mask = t >= 0.02
    np.testing.assert_allclose(y[mask], 1.0 - np.exp(-wc * t[mask]),
                               rtol=1e-3, atol=1e-3)
    np.testing.assert_allclose(y[-1], 1.0, rtol=0.0, atol=1e-4)


def test_filtro_pasobajo2_ganancia_dc():
    m = Modelo(dt=1e-4)
    esc = m.add(FuenteEscalon("u", valor_final=1.0))
    f = m.add(FiltroPasoBajo("f", fc=20.0, zeta=0.707, orden=2))
    m.conectar(esc.salida, f.entrada)
    res = _compara_c_numpy(m, 0.5, [f.salida], atol=1e-9)
    y = np.asarray(res["f"]).ravel()
    np.testing.assert_allclose(y[-1], 1.0, rtol=1e-6, atol=1e-6)


def test_filtro_pasobajo2_subamortiguado_pico():
    # zeta=0.2: el escalon debe superar 1 (sobreimpulso ~52 %)
    m = Modelo(dt=1e-5)
    esc = m.add(FuenteEscalon("u", valor_final=1.0))
    f = m.add(FiltroPasoBajo("f", fc=5.0, zeta=0.2, orden=2))
    m.conectar(esc.salida, f.entrada)
    res = _compara_c_numpy(m, 1.0, [f.salida], atol=1e-9)
    y = np.asarray(res["f"]).ravel()
    pico = np.max(y)
    # sobreimpulso teorico: exp(-pi zeta / sqrt(1-zeta^2)) = 0.5267
    np.testing.assert_allclose(pico, 1.0 + 0.5267, rtol=1e-2, atol=1e-2)


def test_filtro_notch_rechaza_fn():
    m = Modelo(dt=1e-5)
    seno = m.add(FuenteSeno("u", amplitud=1.0, frecuencia=50.0))
    f = m.add(FiltroNotch("f", fn=50.0, zeta=0.05))
    m.conectar(seno.salida, f.entrada)
    res = _compara_c_numpy(m, 1.0, [f.salida], atol=1e-12)
    y = np.asarray(res["f"]).ravel()
    ultimo = y[-int(0.05 / m.dt):]  # ultimos 50 ms
    np.testing.assert_allclose(np.max(np.abs(ultimo)), 0.0,
                               rtol=0.0, atol=1e-3)


def test_filtro_notch_deja_pasar_alejado():
    # a 1 Hz (lejos de fn=50 Hz) la ganancia debe ser ~1
    m = Modelo(dt=1e-4)
    seno = m.add(FuenteSeno("u", amplitud=1.0, frecuencia=1.0))
    f = m.add(FiltroNotch("f", fn=50.0, zeta=0.05))
    m.conectar(seno.salida, f.entrada)
    res = _compara_c_numpy(m, 2.0, [f.salida], atol=1e-12)
    y = np.asarray(res["f"]).ravel()
    ultimo = y[-int(0.5 / m.dt):]
    np.testing.assert_allclose(np.max(np.abs(ultimo)), 1.0,
                               rtol=1e-2, atol=1e-2)


def test_filtro_pasoalto_rechaza_dc():
    m = Modelo(dt=1e-4)
    esc = m.add(FuenteEscalon("u", valor_final=1.0))
    f = m.add(FiltroPasoAlto("f", fc=10.0, orden=1))
    m.conectar(esc.salida, f.entrada)
    res = _compara_c_numpy(m, 0.3, [f.salida], atol=1e-8)
    y = np.asarray(res["f"]).ravel()
    t = np.arange(len(y)) * m.dt
    wc = 2.0 * np.pi * 10.0
    np.testing.assert_allclose(y[0], 0.0, rtol=0.0, atol=1e-12)
    mask = t >= 0.02
    np.testing.assert_allclose(y[mask], np.exp(-wc * t[mask]),
                               rtol=1e-2, atol=1e-3)
    np.testing.assert_allclose(y[-1], 0.0, rtol=0.0, atol=1e-6)


def test_filtro_validacion():
    with pytest.raises(ValueError, match="fc"):
        FiltroPasoBajo("f", fc=0.0)
    with pytest.raises(ValueError, match="orden"):
        FiltroPasoBajo("f", fc=10.0, orden=3)
    with pytest.raises(ValueError, match="fn"):
        FiltroNotch("f", fn=-1.0)


# ---------------------------------------------------------- retenedor
def test_retenedor_captura_en_flancos():
    # senal constante 5, trigger con pulsos: captura 5 en el primer flanco
    m = Modelo(dt=0.01)
    sen = m.add(FuenteEscalon("u", valor_final=5.0, t_paso=0.0))
    trig = m.add(PulsoRectangular("tr", amplitud=1.0, periodo=1.0,
                                  duty=0.5, offset=0.0))
    ret = m.add(RetenedorDisparado("ret", umbral=0.5))
    m.conectar(sen.salida, ret.senal)
    m.conectar(trig.salida, ret.trigger)
    res = _compara_c_numpy(m, 3.0, [ret.salida])
    y = np.asarray(res["ret"]).ravel()
    esperado = np.full_like(y, 5.0)
    esperado[0] = 0.0  # emision inicial: aun no hubo flanco
    np.testing.assert_allclose(y, esperado, rtol=0.0, atol=1e-12)


def test_retenedor_mantiene_valor_entre_flancos():
    # la senal pasa a 3 en t=0.5; los flancos del trigger caen en t=1 y t=2:
    # la salida salta de 0 a 3 en el primer flanco posterior a t=0.5
    m = Modelo(dt=0.01)
    sen = m.add(FuenteEscalon("u", valor_final=3.0, t_paso=0.5))
    trig = m.add(PulsoRectangular("tr", amplitud=1.0, periodo=1.0,
                                  duty=0.5, offset=0.0))
    ret = m.add(RetenedorDisparado("ret", umbral=0.5))
    m.conectar(sen.salida, ret.senal)
    m.conectar(trig.salida, ret.trigger)
    res = _compara_c_numpy(m, 3.0, [ret.salida])
    y = np.asarray(res["ret"]).ravel()
    t = np.arange(len(y)) * m.dt
    esperado = np.where(t < 1.0, 0.0, 3.0)
    np.testing.assert_allclose(y, esperado, rtol=0.0, atol=1e-12)


def test_retenedor_sin_flancos_mantiene_el_valor_inicial():
    m = Modelo(dt=0.01)
    sen = m.add(FuenteEscalon("u", valor_final=3.0, t_paso=0.0))
    trig = m.add(FuenteEscalon("t", valor_final=0.0, t_paso=0.0))
    ret = m.add(RetenedorDisparado("ret", umbral=0.5, valor_inicial=-2.0))
    m.conectar(sen.salida, ret.senal)
    m.conectar(trig.salida, ret.trigger)
    res = _compara_c_numpy(m, 1.0, [ret.salida])
    y = np.asarray(res["ret"]).ravel()
    np.testing.assert_allclose(y, -2.0, rtol=0.0, atol=1e-12)


def test_retenedor_validacion():
    ret = RetenedorDisparado("ret", umbral=0.5)
    assert ret.estados_iniciales[1] == -1.0  # sin flanco previo


# ---------------------------------------------------------- maquina de estados
def test_maquina_estados_transicion_por_umbral():
    # 0 -> 1 cuando u0 >= 2; 1 -> 2 cuando u0 >= 4. Escalon a 5: en el
    # primer paso salta a 1 y en el segundo a 2 (se quedan en 2)
    m = Modelo(dt=0.01)
    maq = m.add(MaquinaEstados("maq", n_estados=3, n_entradas=1,
                               transiciones=[(0, 1, 0, ">=", 2.0),
                                             (1, 2, 0, ">=", 4.0)]))
    esc = m.add(FuenteEscalon("u", valor_final=5.0, t_paso=0.0))
    m.conectar(esc.salida, maq.entradas[0])
    res = _compara_c_numpy(m, 0.5, [maq.salida])
    y = np.asarray(res["maq"]).ravel()
    assert y[0] == 0.0          # emision inicial con el estado inicial
    assert y[1] == 1.0          # primera transicion en el primer paso
    np.testing.assert_allclose(y[2:], 2.0, rtol=0.0, atol=1e-12)


def test_maquina_estados_secuencia_con_pulsos():
    # u: 6 en [0,1), 0 en [1,4), 6 en [4,5), 0 en [5,6).
    # transiciones: u<1 -> 1; u<=0 -> 2; u>5 -> 0.
    # secuencia de estados: 0* -> 1 -> 2* -> 0* -> 1 -> 2*
    m = Modelo(dt=0.01)
    maq = m.add(MaquinaEstados("maq", n_estados=3, n_entradas=1,
                               transiciones=[(0, 1, 0, "<", 1.0),
                                             (1, 2, 0, "<=", 0.0),
                                             (2, 0, 0, ">", 5.0)]))
    pulsos = m.add(PulsoRectangular("u", amplitud=6.0, periodo=4.0,
                                    duty=0.25, offset=0.0))
    m.conectar(pulsos.salida, maq.entradas[0])
    res = _compara_c_numpy(m, 6.0, [maq.salida])
    y = np.asarray(res["maq"]).ravel()
    assert y[0] == 0.0
    i1 = np.argmax(y == 1.0)      # 0->1 en t=1.0
    i2 = np.argmax(y == 2.0)      # 1->2 en t=1.01
    i0 = np.argmax(y[100:] == 0.0) + 100  # 2->0 en t=4.0/4.01
    i1b = np.argmax(y[i0:] == 1.0) + i0   # 0->1 en t=5.0/5.01
    i2b = np.argmax(y[i1b:] == 2.0) + i1b # 1->2 en t=5.01/5.02
    np.testing.assert_allclose([i1, i2], [100, 101], rtol=0.0, atol=1)
    np.testing.assert_allclose([i0, i1b, i2b], [400, 500, 501],
                               rtol=0.0, atol=1)
    assert np.all(y[1:i1] == 0.0)
    assert np.all(y[i1 + 1:i2] == 0.0) or i1 + 1 == i2
    assert np.all(y[i2:i0] == 2.0)
    assert np.all(y[i0:i1b] == 0.0)
    assert np.all(y[i1b:i2b] == 1.0)
    assert np.all(y[i2b:] == 2.0)


def test_maquina_estados_igualdad_y_desigualdad():
    m = Modelo(dt=0.01)
    maq = m.add(MaquinaEstados("maq", n_estados=2, n_entradas=1,
                               transiciones=[(0, 1, 0, "==", 3.0)]))
    esc = m.add(FuenteEscalon("u", valor_final=3.0, t_paso=0.0))
    m.conectar(esc.salida, maq.entradas[0])
    res = _compara_c_numpy(m, 0.2, [maq.salida])
    y = np.asarray(res["maq"]).ravel()
    assert y[0] == 0.0
    np.testing.assert_allclose(y[1:], 1.0, rtol=0.0, atol=1e-12)

    m2 = Modelo(dt=0.01)
    maq2 = m2.add(MaquinaEstados("maq2", n_estados=2, n_entradas=1,
                                 transiciones=[(0, 1, 0, "!=", 2.0)]))
    esc2 = m2.add(FuenteEscalon("u", valor_final=3.0, t_paso=0.0))
    m2.conectar(esc2.salida, maq2.entradas[0])
    res2 = _compara_c_numpy(m2, 0.2, [maq2.salida])
    y2 = np.asarray(res2["maq2"]).ravel()
    np.testing.assert_allclose(y2[1:], 1.0, rtol=0.0, atol=1e-12)


def test_maquina_estados_prioridad():
    # dos transiciones aplicables desde el estado 0: gana la primera
    m = Modelo(dt=0.01)
    maq = m.add(MaquinaEstados("maq", n_estados=3, n_entradas=1,
                               transiciones=[(0, 1, 0, ">", 0.0),
                                             (0, 2, 0, ">", 0.0)]))
    esc = m.add(FuenteEscalon("u", valor_final=1.0, t_paso=0.0))
    m.conectar(esc.salida, maq.entradas[0])
    res = _compara_c_numpy(m, 0.2, [maq.salida])
    y = np.asarray(res["maq"]).ravel()
    np.testing.assert_allclose(y[1:], 1.0, rtol=0.0, atol=1e-12)


def test_maquina_estados_dos_entradas():
    # u0 >= 1 -> estado 1 (usa la entrada 0); u1 >= 2 -> estado 2 (entrada 1)
    m = Modelo(dt=0.01)
    maq = m.add(MaquinaEstados("maq", n_estados=3, n_entradas=2,
                               transiciones=[(0, 1, 0, ">=", 1.0),
                                             (1, 2, 1, ">=", 2.0)]))
    e0 = m.add(FuenteEscalon("e0", valor_final=1.0, t_paso=0.0))
    e1 = m.add(FuenteEscalon("e1", valor_final=2.0, t_paso=0.0))
    m.conectar(e0.salida, maq.entradas[0])
    m.conectar(e1.salida, maq.entradas[1])
    res = _compara_c_numpy(m, 0.3, [maq.salida])
    y = np.asarray(res["maq"]).ravel()
    assert y[1] == 1.0
    np.testing.assert_allclose(y[2:], 2.0, rtol=0.0, atol=1e-12)


def test_maquina_estados_validacion():
    with pytest.raises(ValueError, match="estado_inicial"):
        MaquinaEstados("maq", n_estados=2, n_entradas=1, transiciones=[],
                       estado_inicial=5)
    with pytest.raises(ValueError, match="fuera de rango"):
        MaquinaEstados("maq", n_estados=2, n_entradas=1,
                       transiciones=[(0, 5, 0, ">", 1.0)])
    with pytest.raises(ValueError, match="indice_senal"):
        MaquinaEstados("maq", n_estados=2, n_entradas=1,
                       transiciones=[(0, 1, 7, ">", 1.0)])
    with pytest.raises(ValueError, match="condicion"):
        MaquinaEstados("maq", n_estados=2, n_entradas=1,
                       transiciones=[(0, 1, 0, "~", 1.0)])