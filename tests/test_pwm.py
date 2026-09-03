import numpy as np
from bloques_crysi import Modelo, GeneradorPWM, GeneradorSPWM, GeneradorSVPWM, FuenteConstante
from tests.helpers import correr_c_y_numpy

def test_pwm_1f():
    # Modelo monofásico de duty cycle constante
    m = Modelo(dt=1e-6)
    ref = m.add(FuenteConstante("ref", 0.75))
    pwm = m.add(GeneradorPWM("pwm", fsw=1000.0, dead_time=0.0))
    m.conectar(ref.salida, pwm.d)
    
    res_c, arr_np, _ = correr_c_y_numpy(m, 0.005, [pwm.sensorDisparo()])
    
    # El duty es 0.75, así que el promedio de la señal debe ser ~0.75
    s_c = res_c["S_pwm"]
    s_np = arr_np[0]
    
    assert np.allclose(s_c, s_np)
    
    duty_medido = np.mean(s_c)
    assert np.isclose(duty_medido, 0.75, atol=0.05)


def test_pwm_spwm():
    m = Modelo(dt=1e-6)
    ref = m.add(FuenteConstante("ref", 0.8)) # índice de modulación 0.8
    spwm = m.add(GeneradorSPWM("spwm", f_out=50.0, fsw=2000.0))
    m.conectar(ref.salida, spwm.m)
    
    # Registramos sólo los disparos
    res_c, arr_np, _ = correr_c_y_numpy(m, 0.02, [spwm.sensorDisparos()])
    
    sa_c = res_c["S_spwm"][:, 0]
    sa_np = arr_np[0]
    assert np.allclose(sa_c, sa_np)
    
    # Para m=0.8, en un semiciclo positivo de la moduladora, el duty promedio de Sa 
    # es alto. Como es onda senoidal centrada, el valor medio en 1 ciclo es 0.5
    assert np.isclose(np.mean(sa_c), 0.5, atol=0.01)


def test_pwm_svpwm():
    m = Modelo(dt=1e-6)
    va = m.add(FuenteConstante("va", 100.0))
    vb = m.add(FuenteConstante("vb", 0.0))
    svpwm = m.add(GeneradorSVPWM("svpwm", Vdc=600.0, fsw=2000.0))
    
    m.conectar(va.salida, svpwm.v_alpha)
    m.conectar(vb.salida, svpwm.v_beta)
    
    res_c, arr_np, _ = correr_c_y_numpy(m, 0.005, [svpwm.sensorDisparos()])
    
    sa_c = res_c["S_svpwm"][:, 0]
    sa_np = arr_np[0]
    
    assert np.allclose(sa_c, sa_np)
