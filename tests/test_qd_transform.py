"""Test de regresion para el bloque TransformadaQD (OP_QD)."""

import numpy as np
from bloques_crysi import (
    Modelo, FuenteTrifasica, FuenteConstante, FuenteRampa,
    MaquinaImanesPermanentes, TransformadaQD,
)
from bloques_crysi.backend_numpy import simular
from tests.helpers import correr_c_y_numpy, flat_c, resolver_registro


def _arma_qd_con_maquina(dt=1e-4):
    """Modelo con PMAC + TransformadaQD alimentada por sensores de la maquina."""
    m = Modelo(dt=dt)
    src = m.add(FuenteTrifasica("red", amplitud=310.0, frecuencia=50.0))
    tl = m.add(FuenteConstante("par", 5.0))
    maq = m.add(MaquinaImanesPermanentes("pmsm", rs=0.1, Ld=1e-3, Lq=1e-3,
                                         lam_m=0.1, P=6, J=0.01))
    qd = m.add(TransformadaQD("qd"))
    m.conectar(src.salida, maq.terminales)
    m.conectar(tl.salida, maq.T_L)
    m.conectar(maq.sensor3V(), qd.vabc)
    m.conectar(maq.sensor3I(), qd.iabc)
    m.conectar(maq.sensorPosicionElectrica(), qd.th)
    return m, qd, maq


def test_qd_maquina_coincide():
    """El bloque QD debe coincidir con los sensores iqs/ids internos de la PMAC."""
    m, qd, maq = _arma_qd_con_maquina()
    t_fin = 0.3
    registrar = [qd.salida, maq.sensorCorrienteQ(), maq.sensorCorrienteD()]
    res = m.run(t_fin, registrar=registrar)
    i = 2000  # t = 0.2 s
    iqs_bloque = res["qd"][i][2]
    ids_bloque = res["qd"][i][3]
    iqs_maq = res["iqs"][i]
    ids_maq = res["ids"][i]
    assert abs(iqs_bloque - iqs_maq) < 1e-6
    assert abs(ids_bloque - ids_maq) < 1e-6


def test_qd_c_numpy_coincide():
    """Backend C y numpy deben producir resultados identicos para TransformadaQD."""
    m, qd, maq = _arma_qd_con_maquina()
    t_fin = 0.3
    registrar = [qd.salida]
    res_c, arr_np, n_steps = correr_c_y_numpy(m, t_fin, registrar)
    arr_c = flat_c(res_c, registrar, m, n_steps)
    assert arr_c.shape == arr_np.shape
    np.testing.assert_allclose(arr_c, arr_np, rtol=1e-8, atol=1e-9)


def _arma_qd_fuentes(dt=1e-4):
    """Modelo con fuente trifasica pura + TransformadaQD (sin maquina)."""
    m = Modelo(dt=dt)
    src = m.add(FuenteTrifasica("red", amplitud=310.0, frecuencia=50.0))
    ramp = m.add(FuenteRampa("th", 2 * np.pi * 50.0, 0.0))
    qd = m.add(TransformadaQD("qd"))
    m.conectar(src.salida, qd.vabc)
    m.conectar(src.salida, qd.iabc)
    m.conectar(ramp.salida, qd.th)
    return m, qd, src, ramp


def test_qd_fuentes_manual():
    """TransformadaQD con fuentes puras debe coincidir con el calculo manual exacto."""
    m, qd, src, ramp = _arma_qd_fuentes()
    t_fin = 0.3
    registrar = [qd.salida, src.salida, ramp.salida]
    res = m.run(t_fin, registrar=registrar)
    i = 2000
    # Transformada manual a partir de las senales registradas
    va = res["red"][i][0]
    vb = res["red"][i][1]
    vc = res["red"][i][2]
    ia = res["red"][i][0]
    ib = res["red"][i][1]
    ic = res["red"][i][2]
    th = float(res["th"][i])
    alv = (2.0 / 3.0) * (va - 0.5 * vb - 0.5 * vc)
    bev = (vb - vc) / np.sqrt(3.0)
    ali = (2.0 / 3.0) * (ia - 0.5 * ib - 0.5 * ic)
    bei = (ib - ic) / np.sqrt(3.0)
    dv = alv * np.cos(th) + bev * np.sin(th)
    qv = -alv * np.sin(th) + bev * np.cos(th)
    di = ali * np.cos(th) + bei * np.sin(th)
    qi = -ali * np.sin(th) + bei * np.cos(th)
    vqs_b, vds_b = res["qd"][i][0], res["qd"][i][1]
    iqs_b, ids_b = res["qd"][i][2], res["qd"][i][3]
    assert abs(vqs_b - qv) < 1e-6
    assert abs(vds_b - dv) < 1e-6
    assert abs(iqs_b - qi) < 1e-6
    assert abs(ids_b - di) < 1e-6


def test_qd_fuentes_c_numpy():
    """Backend C y numpy coinciden con fuentes puras."""
    m, qd, src, ramp = _arma_qd_fuentes()
    t_fin = 0.3
    registrar = [qd.salida]
    res_c, arr_np, n_steps = correr_c_y_numpy(m, t_fin, registrar)
    arr_c = flat_c(res_c, registrar, m, n_steps)
    assert arr_c.shape == arr_np.shape
    np.testing.assert_allclose(arr_c, arr_np, rtol=1e-8, atol=1e-9)


def test_qd_dos_bloques_independientes():
    """Dos bloques QD con entradas distintas no deben interferir."""
    m = Modelo(dt=1e-4)
    src = m.add(FuenteTrifasica("red", amplitud=310.0, frecuencia=50.0))
    tl = m.add(FuenteConstante("par", 5.0))
    maq = m.add(MaquinaImanesPermanentes("pmsm", rs=0.1, Ld=1e-3, Lq=1e-3,
                                         lam_m=0.1, P=6, J=0.01))
    ramp = m.add(FuenteRampa("th2", 2 * np.pi * 50.0, 0.0))
    qd = m.add(TransformadaQD("qd"))
    qd2 = m.add(TransformadaQD("qd2"))
    m.conectar(src.salida, maq.terminales)
    m.conectar(tl.salida, maq.T_L)
    m.conectar(maq.sensor3V(), qd.vabc)
    m.conectar(maq.sensor3I(), qd.iabc)
    m.conectar(maq.sensorPosicionElectrica(), qd.th)
    m.conectar(src.salida, qd2.vabc)
    m.conectar(src.salida, qd2.iabc)
    m.conectar(ramp.salida, qd2.th)
    t_fin = 0.3
    registrar = [qd.salida, qd2.salida]
    res_c, arr_np, n_steps = correr_c_y_numpy(m, t_fin, registrar)
    arr_c = flat_c(res_c, registrar, m, n_steps)
    assert arr_c.shape == arr_np.shape
    np.testing.assert_allclose(arr_c, arr_np, rtol=1e-8, atol=1e-9)


if __name__ == "__main__":
    test_qd_maquina_coincide()
    test_qd_c_numpy_coincide()
    test_qd_fuentes_manual()
    test_qd_fuentes_c_numpy()
    test_qd_dos_bloques_independientes()
    print("Todos los tests de TransformadaQD pasaron.")