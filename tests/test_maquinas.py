"""Validación de las 4 máquinas: C == numpy (Euler y RK4) y chequecos físicos."""

import numpy as np
import pytest

from bloques_crysi import (
    Modelo, FuenteTrifasica, FuenteConstante,
    MaquinaImanesPermanentes, MaquinaInduccion, MaquinaSincrona,
    MaquinaCorrienteContinua, MaquinaDCImanesPermanentes,
)

from .helpers import correr_c_y_numpy, flat_c


def _vref_menor(W):
    # voltaje de linea-media por fase para frecuencia f
    return W / (np.sqrt(2) * np.sqrt(3))


def _arma_pmac(dt):
    m = Modelo(dt=dt)
    src = m.add(FuenteTrifasica("red", 310.0, 50.0))
    tl = m.add(FuenteConstante("tl", 5.0))
    maq = m.add(MaquinaImanesPermanentes(
        "pmsm", rs=0.1, Ld=1e-3, Lq=1e-3, lam_m=0.1, P=6, J=0.01, Bm=0.001))
    m.conectar(src.salida, maq.terminales)
    m.conectar(tl.salida, maq.T_L)
    return m, maq


def test_pmac_euler_c_numpy():
    m, maq = _arma_pmac(1e-4)
    _compara(m, 0.5, [maq.sensorVelocidad(), maq.sensorPar(), maq.sensor3I()])


def test_pmac_rk4_c_numpy():
    m, maq = _arma_pmac(2e-4)
    m.metodo = 1
    _compara(m, 0.5, [maq.sensorVelocidad(), maq.sensor3I(), maq.sensor3V()])


def test_pmac_arranca():
    m, maq = _arma_pmac(1e-4)
    res = m.run(0.5, registrar=[maq.sensorVelocidad(), maq.sensorPar()])
    assert res["wm"][-1] > 0
    # en regimen la velocidad creciente con par > T_L inicialmente
    assert np.max(np.abs(res["Te"])) > 5.0


def _arma_mi(dt):
    m = Modelo(dt=dt)
    src = m.add(FuenteTrifasica("red", 310.0, 50.0))
    tl = m.add(FuenteConstante("tl", 0.0))
    maq = m.add(MaquinaInduccion(
        "mi", rs=0.5, rr=0.4, Lm=0.1, Lls=0.005, Llr=0.005, P=4, J=0.5, Bm=0.01))
    m.conectar(src.salida, maq.terminales)
    m.conectar(tl.salida, maq.T_L)
    return m, maq


def test_mi_euler_c_numpy():
    m, maq = _arma_mi(1e-4)
    _compara(m, 0.5, [maq.sensorVelocidad(), maq.sensorPar(), maq.sensor3I()])


def test_mi_rk4_c_numpy():
    m, maq = _arma_mi(2e-4)
    m.metodo = 1
    _compara(m, 0.5, [maq.sensorVelocidad(), maq.sensor3I()])


def test_mi_acelera():
    m, maq = _arma_mi(1e-4)
    res = m.run(1.0, registrar=[maq.sensorVelocidad(), maq.sensorPar()])
    assert res["wm"][-1] > 20


def _arma_sinc(dt):
    m = Modelo(dt=dt)
    src = m.add(FuenteTrifasica("red", 310.0, 50.0))
    vfd = m.add(FuenteConstante("vfd", 10.0))
    tl = m.add(FuenteConstante("tl", 0.0))
    maq = m.add(MaquinaSincrona(
        "ms", rs=0.3, rfd=0.05, rkq1=0.1, rkq2=0.1, rkd=0.1,
        Lls=0.002, Lmq=0.08, Llkq1=0.01, Llkq2=0.005,
        Lmd=0.1, Llf=0.02, Llkd=0.005, P=4, J=1.0, Bm=0.01))
    m.conectar(src.salida, maq.terminales)
    m.conectar(vfd.salida, maq.vfd)
    m.conectar(tl.salida, maq.T_L)
    return m, maq


def test_sinc_euler_c_numpy():
    m, maq = _arma_sinc(1e-4)
    _compara(m, 0.5, [maq.sensorVelocidad(), maq.sensor3I(), maq.sensorPar()])


def test_sinc_rk4_c_numpy():
    m, maq = _arma_sinc(2e-4)
    m.metodo = 1
    _compara(m, 0.5, [maq.sensorVelocidad(), maq.sensor3I()])


def test_sinc_sincroniza_hueco():
    # sin campo (vfd=0) no debe haber par medio; la maquina oscila
    m = Modelo(dt=1e-4)
    src = m.add(FuenteTrifasica("red", 100.0, 50.0))
    vfd = m.add(FuenteConstante("vfd", 0.0))
    tl = m.add(FuenteConstante("tl", 0.0))
    maq = m.add(MaquinaSincrona(
        "ms", rs=0.3, rfd=0.05, rkq1=0.1, rkq2=0.1, rkd=0.1,
        Lls=0.002, Lmq=0.08, Llkq1=0.01, Llkq2=0.005,
        Lmd=0.1, Llf=0.02, Llkd=0.005, P=4, J=1.0, Bm=0.01))
    m.conectar(src.salida, maq.terminales)
    m.conectar(vfd.salida, maq.vfd)
    m.conectar(tl.salida, maq.T_L)
    res = m.run(1.0, registrar=[maq.sensorVelocidad()])
    # velocidad media baja o negativa (no sincroniza / frena por rozamiento)
    assert res["wm"][-1] < 200


def _arma_cc(dt):
    m = Modelo(dt=1e-5)
    va = m.add(FuenteConstante("va", 120.0))
    vf = m.add(FuenteConstante("vf", 100.0))
    tl = m.add(FuenteConstante("tl", 0.0))
    maq = m.add(MaquinaCorrienteContinua(
        "cc", r_a=1.0, L_a=0.01, r_f=100.0, L_f=10.0, L_AF=1.5, J=0.5))
    m.conectar(va.salida, maq.entrada)
    m.conectar(vf.salida, maq.campo)
    m.conectar(tl.salida, maq.T_L)
    return m, maq


def test_cc_euler_c_numpy():
    m, maq = _arma_cc(1e-5)
    _compara(m, 0.5, [maq.sensorCorriente(), maq.sensorVelocidad(), maq.sensorPar()])


def test_cc_rk4_c_numpy():
    m, maq = _arma_cc(2e-5)
    m.metodo = 1
    _compara(m, 0.5, [maq.sensorCorriente(), maq.sensorVelocidad(), maq.sensorPar()])


def test_cc_te_fisico():
    m, maq = _arma_cc(1e-5)
    res = m.run(0.5, registrar=[maq.sensorCorriente(), maq.sensorCampo(),
                                maq.sensorPar(), maq.sensorVelocidad()])
    # Te = LAF * if * ia exacto
    np.testing.assert_allclose(
        res["Te"], 1.5 * res["if"] * res["ia"], rtol=1e-6, atol=1e-9)
    # en regimen (dia/dt -> 0): residual de la ecuacion de armadura pequeno
    n = len(res["ia"])
    residuo = np.abs(120.0 - 1.0 * res["ia"] - 1.5 * res["if"] * res["wm"])
    assert residuo[-1] < 2.0
    assert res["wm"][-1] > 0


def test_cc_ext_c_numpy():
    m = Modelo(dt=1e-5)
    from bloques_crysi.puertos import Puerto
    va = m.add(FuenteConstante("va", 0.0))
    vf = m.add(FuenteConstante("vf", 10.0))
    wm = m.add(FuenteConstante("wm_ext", 50.0))
    th = m.add(FuenteConstante("th_ext", 0.2))
    maq = m.add(MaquinaCorrienteContinua(
        "cc", r_a=0.5, L_a=0.01, r_f=1.0, L_f=0.1, L_AF=0.5, J=0.02,
        mecanica_interna=False))
    m.conectar(va.salida, maq.entrada)
    m.conectar(vf.salida, maq.campo)
    m.conectar(wm.salida, Puerto(maq, "ent", 2, 1))
    m.conectar(th.salida, Puerto(maq, "ent", 3, 1))
    _compara(m, 0.5, [maq.sensorCorriente(), maq.sensorCampo(),
                      maq.sensorVelocidad(), maq.sensorPar(),
                      maq.sensorEa(), maq.sensorVoltajeTerminal()])


def test_cc_ea_vt_fisico():
    m, maq = _arma_cc(1e-5)
    res = m.run(0.5, registrar=[maq.sensorCorriente(), maq.sensorCampo(),
                                maq.sensorVelocidad(), maq.sensorEa(),
                                maq.sensorVoltajeTerminal()])
    # Ea = LAF * if * wm exacto
    np.testing.assert_allclose(
        res["Ea"], 1.5 * res["if"] * res["wm"], rtol=1e-6, atol=1e-9)
    # V_t = Ea + ra*ia exacto
    np.testing.assert_allclose(
        res["V_t"], res["Ea"] + 1.0 * res["ia"], rtol=1e-6, atol=1e-9)
    # en regimen V_t = va impuesto
    assert abs(res["V_t"][-1] - 120.0) < 2.0
    assert res["Ea"][-1] > 0.0


def _compara(modelo, t_fin, registrar, rtol=1e-9, atol=1e-10):
    from .helpers import correr_c_y_numpy, flat_c
    res_c, arr_np, n_steps = correr_c_y_numpy(modelo, t_fin, registrar)
    arr_c = flat_c(res_c, registrar, modelo, n_steps)
    assert arr_c.shape == arr_np.shape
    np.testing.assert_allclose(arr_c, arr_np, rtol=rtol, atol=atol)


# ------------------------------------------- DC de imanes permanentes (PMDC)

def _arma_dcpm(dt, ext=False):
    m = Modelo(dt=dt)
    va = m.add(FuenteConstante("va", 50.0))
    maq = m.add(MaquinaDCImanesPermanentes(
        "pmdc", r_a=0.5, L_a=0.01, Kt=1.0, J=0.02, Bm=0.001,
        mecanica_interna=not ext))
    m.conectar(va.salida, maq.entrada)
    if ext:
        wm = m.add(FuenteConstante("wm_ext", 30.0))
        th = m.add(FuenteConstante("th_ext", 0.2))
        from bloques_crysi.puertos import Puerto
        m.conectar(wm.salida, Puerto(maq, "ent", 1, 1))
        m.conectar(th.salida, Puerto(maq, "ent", 2, 1))
    else:
        tl = m.add(FuenteConstante("tl", 1.0))
        m.conectar(tl.salida, maq.T_L)
    return m, maq


def test_dcpm_euler_c_numpy():
    m, maq = _arma_dcpm(1e-4)
    _compara(m, 0.5, [maq.sensorCorriente(), maq.sensorVelocidad(),
                      maq.sensorPar(), maq.sensorEa(),
                      maq.sensorVoltajeTerminal()])


def test_dcpm_rk4_c_numpy():
    m, maq = _arma_dcpm(2e-4)
    m.metodo = 1
    _compara(m, 0.5, [maq.sensorCorriente(), maq.sensorVelocidad(),
                      maq.sensorPar(), maq.sensorEa(),
                      maq.sensorVoltajeTerminal()])


def test_dcpm_ext_c_numpy():
    m, maq = _arma_dcpm(1e-4, ext=True)
    _compara(m, 0.5, [maq.sensorCorriente(), maq.sensorPar(),
                      maq.sensorEa(), maq.sensorVelocidad()])


def test_dcpm_equivale_cc_excitada():
    # PMDC con Kt == CC con excitación separada L_AF*if (vf=rf=1 -> if=1):
    # deben ser idénticas (Ea = Kt*wm en ambos casos). Se fija el estado
    # inicial del campo if(0)=1 para que if quede constante en todo el run.
    def _corre(bloque, con_campo):
        m = Modelo(dt=1e-5)
        va = m.add(FuenteConstante("va", 50.0))
        tl = m.add(FuenteConstante("tl", 1.0))
        maq = m.add(bloque)
        m.conectar(va.salida, maq.entrada)
        if con_campo:
            vf = m.add(FuenteConstante("vf", 1.0))
            m.conectar(vf.salida, maq.campo)
        m.conectar(tl.salida, maq.T_L)
        return m.run(0.5, registrar=[maq.sensorCorriente(), maq.sensorVelocidad(),
                                     maq.sensorPar(), maq.sensorEa(),
                                     maq.sensorVoltajeTerminal()])

    cc = MaquinaCorrienteContinua(
        "cc", r_a=0.5, L_a=0.01, r_f=1.0, L_f=0.1, L_AF=1.0, J=0.02, Bm=0.001)
    cc.estados_iniciales = [0.0, 1.0, 0.0, 0.0]
    pm = MaquinaDCImanesPermanentes(
        "pmdc", r_a=0.5, L_a=0.01, Kt=1.0, J=0.02, Bm=0.001)
    res_cc = _corre(cc, con_campo=True)
    res_pm = _corre(pm, con_campo=False)
    for ky in ("ia", "wm", "Te", "Ea", "V_t"):
        np.testing.assert_allclose(res_pm[ky], res_cc[ky], rtol=1e-9, atol=1e-9)


def test_dcpm_fisica():
    m, maq = _arma_dcpm(1e-5)
    res = m.run(0.5, registrar=[maq.sensorCorriente(), maq.sensorVelocidad(),
                                maq.sensorEa(), maq.sensorVoltajeTerminal(),
                                maq.sensorPar()])
    # Ea = Kt*wm y V_t = Ea + ra*ia exactos por construccion
    np.testing.assert_allclose(res["Ea"], 1.0 * res["wm"], rtol=1e-6, atol=1e-9)
    np.testing.assert_allclose(res["V_t"], res["Ea"] + 0.5 * res["ia"],
                               rtol=1e-6, atol=1e-9)
    # regimen: Te = Kt*ia > TL (friccion) y Ea ~ va - ra*ia
    assert np.mean(res["Te"][-50:]) > 1.0
    assert abs(res["V_t"][-1] - 50.0) < 2.0