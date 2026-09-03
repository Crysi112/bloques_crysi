import numpy as np
from bloques_crysi import (Modelo, MaquinaImanesPermanentes, MaquinaCorrienteContinua,
                           FuenteTrifasica, FuenteConstante)
from .helpers import correr_c_y_numpy, flat_c


def _arma_pmac_cc_acopladas(dt=1e-5):
    m = Modelo(dt=dt)
    red = m.add(FuenteTrifasica("red", amplitud=310.0, frecuencia=50.0))
    pmac = m.add(MaquinaImanesPermanentes("pmac", rs=0.1, Ld=1e-3, Lq=1e-3,
                                          lam_m=0.1, P=4, J=0.01,
                                          mecanica_interna=False))
    m.conectar(red.salida, pmac.terminales)
    v_cc = m.add(FuenteConstante("vcc", 0.0))
    vf_cc = m.add(FuenteConstante("vfcc", 10.0))
    cc = m.add(MaquinaCorrienteContinua("cc", r_a=0.5, L_a=0.01, r_f=1.0,
                                        L_f=0.1, L_AF=0.5, J=0.02,
                                        mecanica_interna=False))
    m.conectar(v_cc.salida, cc.entrada)
    m.conectar(vf_cc.salida, cc.campo)
    eje = m.acoplar_maquinas(pmac, cc, J_eq=0.03, Bm_eq=0.001)
    return m, pmac, cc, eje


def test_acople_eje_c_numpy():
    # el eje mecanico (OP_EJE_MECANICO) debe integrarse igual en el backend
    # numpy que en el C (antes numpy no lo registraba: wm quedaba en cero)
    m, pmac, cc, eje = _arma_pmac_cc_acopladas()
    s_w = pmac.sensorVelocidad()
    s_w.nombre = "w"
    s_T = cc.sensorPar()
    s_T.nombre = "T_cc"
    reg = [s_w, s_T, eje.salida]
    res_c, arr_np, n_steps = correr_c_y_numpy(m, 0.1, reg)
    arr_c = flat_c(res_c, reg, m, n_steps)
    np.testing.assert_allclose(arr_c, arr_np, rtol=1e-9, atol=1e-10)
    # y con el fix el eje realmente gira (no se queda en 0)
    assert np.max(np.abs(arr_c[0])) > 10.0

def test_acoplamiento_mecanico_pmac_cc():
    """Conecta una PMAC a la red (actúa de motor) acoplada a una CC en vacío (actúa de generador o carga inercial)."""
    m = Modelo(dt=1e-5)
    
    red = m.add(FuenteTrifasica("red", amplitud=310.0, frecuencia=50.0))
    # PMAC como motor
    pmac = m.add(MaquinaImanesPermanentes("pmac", rs=0.1, Ld=1e-3, Lq=1e-3, lam_m=0.1, P=4, J=0.01, mecanica_interna=False))
    m.conectar(red.salida, pmac.terminales)
    
    # CC como generador (va=0, excitada)
    v_cc = m.add(FuenteConstante("vcc", 0.0))
    vf_cc = m.add(FuenteConstante("vfcc", 10.0))
    cc = m.add(MaquinaCorrienteContinua("cc", r_a=0.5, L_a=0.01, r_f=1.0, L_f=0.1, L_AF=0.5, J=0.02, mecanica_interna=False))
    m.conectar(v_cc.salida, cc.entrada)
    m.conectar(vf_cc.salida, cc.campo)
    
    # Acople mecanico: J_eq = J_pmac + J_cc = 0.03
    eje = m.acoplar_maquinas(pmac, cc, J_eq=0.03, Bm_eq=0.001)
    
    s_w_pmac = pmac.sensorVelocidad()
    s_w_pmac.nombre = "w_pmac"
    s_w_cc = cc.sensorVelocidad()
    s_w_cc.nombre = "w_cc"
    s_T_pmac = pmac.sensorPar()
    s_T_pmac.nombre = "T_pmac"
    s_T_cc = cc.sensorPar()
    s_T_cc.nombre = "T_cc"
    
    res = m.run(0.2, registrar=[s_w_pmac, s_w_cc, s_T_pmac, s_T_cc])
    
    w_pmac = res["w_pmac"]
    w_cc = res["w_cc"]
    Te_pmac = res["T_pmac"]
    Te_cc = res["T_cc"]
    
    # Verificamos que compartan la misma velocidad exactamente
    np.testing.assert_allclose(w_pmac, w_cc, atol=1e-9)
    
    # Verificamos que la maquina CC gira y por tanto tiene un Te_cc negativo o cero, pero como gira con la PMAC generará un par opuesto si va=0
    assert w_cc[-1] > 10.0 # PMAC lo acelera
    assert Te_cc[-1] < -0.1 # La CC frena a la PMAC
    

def test_cc_generador_ea_vt():
    """La CC acoplada opera como generador: Ea = LAF*if*wm y V_t = Ea + ra*ia < Ea."""
    m = Modelo(dt=1e-5)
    red = m.add(FuenteTrifasica("red", amplitud=310.0, frecuencia=50.0))
    pmac = m.add(MaquinaImanesPermanentes("pmac", rs=0.1, Ld=1e-3, Lq=1e-3,
                                          lam_m=0.1, P=4, J=0.01,
                                          mecanica_interna=False))
    m.conectar(red.salida, pmac.terminales)
    v_cc = m.add(FuenteConstante("vcc", 0.0))
    vf_cc = m.add(FuenteConstante("vfcc", 10.0))
    cc = m.add(MaquinaCorrienteContinua("cc", r_a=0.5, L_a=0.01, r_f=1.0,
                                        L_f=0.1, L_AF=0.5, J=0.02,
                                        mecanica_interna=False))
    m.conectar(v_cc.salida, cc.entrada)
    m.conectar(vf_cc.salida, cc.campo)
    m.acoplar_maquinas(pmac, cc, J_eq=0.03, Bm_eq=0.001)

    s_ia = cc.sensorCorriente()
    s_ia.nombre = "ia_cc"
    s_w = cc.sensorVelocidad()
    s_w.nombre = "w_cc"
    s_ea = cc.sensorEa()
    s_ea.nombre = "Ea_cc"
    s_vt = cc.sensorVoltajeTerminal()
    s_vt.nombre = "Vt_cc"

    res = m.run(0.2, registrar=[s_ia, s_w, s_ea, s_vt,
                                cc.sensorCampo()])

    np.testing.assert_allclose(
        res["Ea_cc"], 0.5 * res["if"] * res["w_cc"], rtol=1e-6, atol=1e-9)
    np.testing.assert_allclose(
        res["Vt_cc"], res["Ea_cc"] + 0.5 * res["ia_cc"], rtol=1e-6, atol=1e-9)
    # generadora: corriente negativa, Ea > 0 y V_t = Ea - ra*|ia| < Ea
    assert res["ia_cc"][-1] < 0.0
    assert res["Ea_cc"][-1] > 0.0
    assert res["Vt_cc"][-1] < res["Ea_cc"][-1]
    
