"""Validación de módulos de potencia: C == numpy y física (todo promediado)."""

import numpy as np
import pytest

from bloques_crysi import (
    Modelo, Puerto, FuenteConstante, FuenteSeno, FuenteTrifasica,
    ConvertidorBuck, ConvertidorBoost, ConvertidorBuckBoost,
    RectificadorTrifasico, InversorTrifasico, CargaRLTrifasica,
    Transformador, MedidorPotencia, MaquinaImanesPermanentes,
    MaquinaCorrienteContinua,
)
from .helpers import correr_c_y_numpy, flat_c

def _compara(modelo, t_fin, registrar, rtol=1e-9, atol=1e-10, metodo=0):
    modelo.metodo = metodo
    res_c, arr_np, n_steps = correr_c_y_numpy(modelo, t_fin, registrar)
    arr_c = flat_c(res_c, registrar, modelo, n_steps)
    assert arr_c.shape == arr_np.shape, (arr_c.shape, arr_np.shape)
    np.testing.assert_allclose(arr_c, arr_np, rtol=rtol, atol=atol)
    return arr_c, arr_np


# ---------------------------------------------------------- DC-DC

def _arma_dcdc(Cls, vin=100.0, d=0.6, L=5e-3, C=1e-3, R=10.0):
    m = Modelo(dt=1e-5)
    v = m.add(FuenteConstante("vin", vin))
    dc = m.add(FuenteConstante("d", d))
    conv = m.add(Cls("conv", L=L, C=C, R=R))
    m.conectar(v.salida, conv.entrada)
    m.conectar(dc.salida, conv.d)
    return m, conv


def test_buck_c_numpy_euler():
    m, c = _arma_dcdc(ConvertidorBuck)
    _compara(m, 0.2, [c.salida])


def test_buck_c_numpy_rk4():
    m, c = _arma_dcdc(ConvertidorBuck)
    _compara(m, 0.2, [c.salida], metodo=1)


def test_buck_fisica():
    # regimen: L diL/dt=0 -> vout = d*vin ; iL = vout/R
    m, c = _arma_dcdc(ConvertidorBuck, vin=100.0, d=0.6)
    arr, _ = _compara(m, 0.2, [c.salida], atol=1e-8)
    vout = arr[0, -1]
    iL = arr[1, -1]
    assert abs(vout - 60.0) < 0.5
    assert abs(iL - 6.0) < 0.1


def test_boost_fisica():
    # vout = vin/(1-d) = 100/0.4 = 250 ; iL = vout/(R*(1-d))
    m, c = _arma_dcdc(ConvertidorBoost, vin=100.0, d=0.6)
    arr, _ = _compara(m, 0.5, [c.salida], atol=1e-8)
    assert abs(arr[0, -1] - 250.0) < 2.0
    assert abs(arr[1, -1] - 250.0 / (10.0 * (1 - 0.6))) < 0.5


def test_buckboost_c_numpy():
    m, c = _arma_dcdc(ConvertidorBuckBoost)
    _compara(m, 0.3, [c.salida])


# ---------------------------------------------------------- AC-DC

def _arma_rect(v_amp=310.0, C=1e-3, R=100.0, Rint=0.05, t=0.3):
    m = Modelo(dt=1e-4)
    src = m.add(FuenteTrifasica("red", v_amp, 50.0))
    rect = m.add(RectificadorTrifasico("rect", C=C, R=R, Rint=Rint))
    m.conectar(src.salida, rect.entrada)
    return m, rect


def test_rect_c_numpy():
    m, rect = _arma_rect()
    _compara(m, 0.2, [rect.salida])


def test_rect_fisica():
    # regimen: vdc ~ sqrt(3)*Vf (max linea-linea) con poca carga
    m, rect = _arma_rect(v_amp=100.0, C=1e-2, R=1000.0, Rint=0.05)
    arr, _ = _compara(m, 0.6, [rect.salida], atol=1e-8)
    vdc = arr[0, -1]
    assert abs(vdc - np.sqrt(3) * 100.0) < 10.0, vdc


def test_rect_envolvente_instantanea():
    """Reconstruye vrec(t) a partir de idc y vdc (vdc = ich*Rint + vC en
    conduccion continua) y lo compara punto a punto contra la envolvente
    fisica real max(va,vb,vc) - min(va,vb,vc), NO contra max(vab,vbc,vca).

    Carga pesada (R chico) para forzar conduccion continua (idc > 0 en
    todo el ciclo) y poder despejar vrec sin ambiguedad en cada muestra.
    Con la formula max(vab,vbc,vca) el error punta a 15-20% de Vpk en la
    mitad del ciclo; con la envolvente correcta el error es solo ruido
    numerico de integracion.
    """
    Vpk, f, Rint = 100.0, 50.0, 0.02
    m, rect = _arma_rect(v_amp=Vpk, C=2e-3, R=1.0, Rint=Rint, t=0.3)
    m.dt = 1e-5
    res = m.run(0.3, registrar=[rect.salida])
    t = res.t
    vdc = res["rect"][:, 0]
    idc = res["rect"][:, 1]

    k = t > 0.28  # ultimo ciclo, en regimen permanente
    assert idc[k].min() > 1e-3, "la carga elegida no fuerza conduccion continua"

    va = Vpk * np.sin(2 * np.pi * f * t[k])
    vb = Vpk * np.sin(2 * np.pi * f * t[k] - 2 * np.pi / 3)
    vc = Vpk * np.sin(2 * np.pi * f * t[k] + 2 * np.pi / 3)
    vrec_teorico = np.maximum.reduce([va, vb, vc]) - np.minimum.reduce([va, vb, vc])
    vrec_reconstruido = idc[k] * Rint + vdc[k]

    err_max = np.max(np.abs(vrec_reconstruido - vrec_teorico))
    assert err_max < 5.0, (
        f"envolvente del rectificador no coincide con max-min "
        f"(error max = {err_max:.2f} V); revisar formula de vrec "
        f"(debe ser max(va,vb,vc)-min(va,vb,vc), no max(vab,vbc,vca))"
    )


def test_rect_rizado_carga_real():
    """Con parametros de bus tipicos (311 Vpk, C=1000uF, R=50 ohm) el
    rizado pico-pico en regimen debe quedar acotado a lo que predice la
    envolvente fisica real. Con la formula max(vab,vbc,vca) el rizado
    sale ~2x mas grande de lo fisicamente correcto.
    """
    m, rect = _arma_rect(v_amp=311.0, C=1000e-6, R=50.0, Rint=1e-3, t=0.08)
    m.dt = 1e-6
    res = m.run(0.08, registrar=[rect.salida])
    t, vdc = res.t, res["rect"][:, 0]
    k = t > 0.06  # ultimo ciclo (16.6 ms a 60 Hz, dejamos margen)
    ripple = vdc[k].max() - vdc[k].min()
    assert ripple < 30.0, f"rizado de {ripple:.1f} Vpp, esperado < 30 V"


# ---------------------------------------------------------- DC-AC

def _arma_inv(conmutada=False, f_out=60.0, fsw=5000.0, mod=1.0, Lf=5e-3,
              Cf=20e-6, R=50.0, dt=5e-6, t_fin=0.05):
    m = Modelo(dt=dt)
    bus = m.add(FuenteConstante("vdc", 400.0))
    inv = m.add(InversorTrifasico("inv", f_out=f_out, fsw=fsw,
                                  m_start=mod, m_end=mod, Lf=Lf, Cf=Cf, R=R,
                                  conmutada=conmutada))
    m.conectar(bus.salida, inv.entrada)
    return m, inv, t_fin


def test_inv_promediado_c_numpy():
    m, inv, tf = _arma_inv(conmutada=False, t_fin=0.05)
    _compara(m, tf, [inv.salida])


def test_inv_promediado_fisica():
    # vC ~ m*vdc/2 = 200 V pico (fundamental tras LC), cargado en R
    m, inv, tf = _arma_inv(conmutada=False, t_fin=0.2)
    arr, _ = _compara(m, tf, [inv.salida])
    vCa = arr[0, -4000:]  # 1.2 periodos a 60 Hz (dt=5us)
    amp = np.max(np.abs(vCa - vCa.mean()))
    assert 150.0 < amp < 240.0, amp


def test_inv_conmutado_c_numpy_euler():
    # a fsw=5000 la comparacion con dt=5e-6 resuelve la portadora
    m, inv, tf = _arma_inv(conmutada=True, fsw=5000.0, t_fin=0.04)
    _compara(m, tf, [inv.salida])


def test_inv_conmutado_tiene_armonicos():
    # conmutado: la salida tras LC tiene ripple de fsw (mas que el promedio)
    m, inv, tf = _arma_inv(conmutada=True, t_fin=0.1)
    res = m.run(tf, registrar=[inv.sensorVoltajesSalida()])
    v = np.asarray(res["V_f"])[:, 0]
    n = len(v)
    # valor pico-pico del ripple tras la portadora (filtrado, acotado)
    w = v[-200:]
    pp = np.max(w) - np.min(w)
    assert pp < 100.0  # filtrado: ripple acotado

def test_encadenado_rect_inv():
    """Rectificador -> bus -> inversor (promediado): el bus alimenta al inversor."""
    from bloques_crysi import Puerto
    m = Modelo(dt=1e-4)
    src = m.add(FuenteTrifasica("red", 230.0, 50.0))
    rect = m.add(RectificadorTrifasico("rect", C=1e-2, R=100.0, Rint=0.05))
    inv = m.add(InversorTrifasico("inv", f_out=50.0, fsw=5000.0,
                                  m_start=0.9, m_end=0.9, Lf=5e-3,
                                  Cf=20e-6, R=30.0, conmutada=False))
    m.conectar(src.salida, rect.entrada)
    m.conectar(Puerto(rect, "sal", 0, 1), inv.entrada)  # vdc del bus al inversor
    reg = [rect.salida, inv.salida]
    res_c, arr_np, ns = correr_c_y_numpy(m, 0.2, reg)
    arr_c = flat_c(res_c, reg, m, ns)
    np.testing.assert_allclose(arr_c, arr_np, rtol=1e-9, atol=1e-10)
    vdc = np.asarray(res_c["rect"])[:, 0]
    vCa = np.asarray(res_c["inv"])[:, 0]
    assert vdc[-1] > 250.0
    assert np.max(np.abs(vCa[-1000:])) > 50.0


# ---------------------------------------------------------- Inversor monofásico

def _arma_inv1f(conmutada=False, fsw=10000.0, t_fin=0.1, m=0.8, dt=2e-6):
    from bloques_crysi import InversorMonofasico
    m2 = Modelo(dt=dt)
    bus = m2.add(FuenteConstante("vdc", 400.0))
    inv = m2.add(InversorMonofasico("inv", f_out=50.0, fsw=fsw,
                                    m_start=m, m_end=m,
                                    Lf=5e-3, Cf=20e-6, R=50.0,
                                    conmutada=conmutada))
    m2.conectar(bus.salida, inv.entrada)
    return m2, inv, t_fin


def test_inv1f_promediado_c_numpy():
    m, inv, tf = _arma_inv1f()
    _compara(m, tf, [inv.salida])


def test_inv1f_promediado_rk4():
    m, inv, tf = _arma_inv1f()
    _compara(m, tf, [inv.salida], metodo=1)


def test_inv1f_fisica():
    # vC ~ m*vdc/2 = 160 V pico, senoidal a 50 Hz
    m, inv, tf = _arma_inv1f()
    res = m.run(tf, registrar=[inv.salida])
    t, vC = res.t, res["inv"][:, 0]
    k = t > 0.05
    amp = (np.max(vC[k]) - np.min(vC[k])) / 2
    assert 150.0 < amp < 175.0, amp
    # frecuencia: ~5 ciclos en 0.1 s
    cruces = np.where(np.diff(np.sign(vC[k])) != 0)[0]
    f_est = len(cruces) / 2 / (t[k][-1] - t[k][0])
    assert abs(f_est - 50.0) < 2.0


def test_inv1f_conmutado_c_numpy():
    m, inv, tf = _arma_inv1f(conmutada=True)
    _compara(m, tf, [inv.salida], rtol=1e-6, atol=1e-8)


def test_inv1f_conmutado_fisica():
    # SPWM bipolar: fundamental ~ m*vdc/2, sin componente DC apreciable
    m, inv, tf = _arma_inv1f(conmutada=True)
    res = m.run(tf, registrar=[inv.salida])
    t, vC = res.t, res["inv"][:, 0]
    k = (t >= 0.06) & (t < 0.10)   # 2 ciclos exactos de 50 Hz: media = DC pura
    dc = np.mean(vC[k])
    amp = (np.max(vC[k]) - np.min(vC[k])) / 2
    assert abs(dc) < 5.0, dc          # sin DC (bipolar)
    assert 150.0 < amp < 180.0, amp   # fundamental ~ 161 + ripple


def test_inv1f_rampa():
    # rampa de modulacion: la amplitud crece desde m_start hacia m_end
    from bloques_crysi import InversorMonofasico
    m = Modelo(dt=1e-5)
    bus = m.add(FuenteConstante("vdc", 400.0))
    inv = m.add(InversorMonofasico("inv", f_out=50.0, fsw=10000.0,
                                   m_start=0.2, m_end=1.0, t_ramp=0.05,
                                   Lf=5e-3, Cf=20e-6, R=50.0))
    m.conectar(bus.salida, inv.entrada)
    res = m.run(0.2, registrar=[inv.salida])
    t, vC = res.t, res["inv"][:, 0]
    k1 = t < 0.02          # m(0.02) = 0.52 -> amp ~ 105
    k2 = (t >= 0.18) & (t < 0.2)  # m = 1 -> amp ~ 202
    amp1 = (np.max(vC[k1]) - np.min(vC[k1])) / 2
    amp2 = (np.max(vC[k2]) - np.min(vC[k2])) / 2
    assert amp1 < 130.0, amp1
    assert 185.0 < amp2 < 220.0, amp2
    assert amp2 > amp1 * 1.8


def test_inv3f_conmutado_bipolar_sin_dc():
    # SPWM conmutado del inversor trifasico debe ser bipolar (±vdc/2) como el
    # monofasico: sin componente DC en vCa (el modo 0..vdc inyectaba vdc/2 al
    # filtro LC referido a 0 V). Con m=1, vdc=400: fundamental ~ 200 V.
    m, inv, tf = _arma_inv(conmutada=True, mod=1.0, t_fin=0.2)
    res = m.run(tf, registrar=[inv.salida])
    t, vC = res.t, res["inv"][:, 0]
    k = t >= 0.1667  # 2 ciclos exactos de 60 Hz: la media mide solo DC
    dc = np.mean(vC[k])
    amp = (np.max(vC[k]) - np.min(vC[k])) / 2
    assert abs(dc) < 5.0, f"componente DC de {dc:.1f} V en el modo conmutado"
    assert 150.0 < amp < 260.0, amp


def test_inv3f_rampa_descendente():
    # rampa de modulacion descendente (m 1.0 -> 0.4): la amplitud debe bajar
    # gradualmente; el clamp antiguo (mt > m1 -> m1) la cortaba en un escalon
    m = Modelo(dt=5e-6)
    bus = m.add(FuenteConstante("vdc", 400.0))
    inv = m.add(InversorTrifasico("inv", f_out=60.0, fsw=5000.0,
                                  m_start=1.0, m_end=0.4, t_ramp=0.5,
                                  Lf=5e-3, Cf=20e-6, R=50.0))
    m.conectar(bus.salida, inv.entrada)
    res = m.run(1.0, registrar=[inv.salida])
    t, vC = res.t, res["inv"][:, 0]
    k_rampa = (t >= 0.2) & (t < 0.3)          # m(t) ~ 0.7, LC ya asentado
    k_fin = (t >= 0.7) & (t < 1.0)            # m = 0.4 (rampa termina en 0.5)
    amp_rampa = (np.max(vC[k_rampa]) - np.min(vC[k_rampa])) / 2
    amp_fin = (np.max(vC[k_fin]) - np.min(vC[k_fin])) / 2
    assert amp_rampa > amp_fin * 1.3, (amp_rampa, amp_fin)
    assert 60.0 < amp_fin < 110.0, amp_fin


def test_inv1f_rampa_descendente():
    # lo mismo para el inversor monofasico (m 1.0 -> 0.4)
    from bloques_crysi import InversorMonofasico
    m = Modelo(dt=1e-5)
    bus = m.add(FuenteConstante("vdc", 400.0))
    inv = m.add(InversorMonofasico("inv", f_out=50.0, fsw=10000.0,
                                   m_start=1.0, m_end=0.4, t_ramp=0.5,
                                   Lf=5e-3, Cf=20e-6, R=50.0))
    m.conectar(bus.salida, inv.entrada)
    res = m.run(1.0, registrar=[inv.salida])
    t, vC = res.t, res["inv"][:, 0]
    k_rampa = (t >= 0.2) & (t < 0.3)          # m(t) ~ 0.7, LC ya asentado
    k_fin = (t >= 0.7) & (t < 1.0)            # m = 0.4 (rampa termina en 0.5)
    amp_rampa = (np.max(vC[k_rampa]) - np.min(vC[k_rampa])) / 2
    amp_fin = (np.max(vC[k_fin]) - np.min(vC[k_fin])) / 2
    assert amp_rampa > amp_fin * 1.3, (amp_rampa, amp_fin)
    assert 60.0 < amp_fin < 110.0, amp_fin


# ---------------------------------------------------------- Carga RL trifásica

def _arma_carga_rl(R=10.0, L=10e-3, V=311.0, f=50.0, t_fin=0.5):
    m = Modelo(dt=1e-4)
    ft = m.add(FuenteTrifasica("ft", amplitud=V, frecuencia=f))
    c = m.add(CargaRLTrifasica("c", R=R, L=L))
    m.conectar(ft.salida, c.entrada)
    return m, c, t_fin, V, f


def test_carga_rl_c_numpy_euler():
    m, c, tf, _, _ = _arma_carga_rl()
    _compara(m, tf, [c.salida])


def test_carga_rl_c_numpy_rk4():
    m, c, tf, _, _ = _arma_carga_rl()
    _compara(m, tf, [c.salida], metodo=1)


def test_carga_rl_fisica():
    # regimen: I = V/|Z| con Z = R + j w L ; ia+ib+ic = 0
    m, c, tf, V, f = _arma_carga_rl()
    res = m.run(tf, registrar=[c.salida])
    t, i = res.t, res["c"]
    k = t > 0.3
    Z = np.sqrt(10.0 ** 2 + (2 * np.pi * f * 10e-3) ** 2)
    ip = np.max(np.abs(i[k, 0]))
    assert abs(ip - V / Z) / (V / Z) < 0.02, ip
    # suma de fases = 0 en todo momento (estrella sin neutro)
    suma = i[:, 0] + i[:, 1] + i[:, 2]
    assert np.max(np.abs(suma)) < 1e-9
    # desfase: i atrasada ~ atan(wL/R); el integrador añade retraso ZOH
    # (~h*w con Euler, ~h*w/2 con RK4)
    va = V * np.sin(2 * np.pi * f * t)
    phi = np.arccos(np.clip(
        np.dot(i[k, 0], va[k]) / np.linalg.norm(i[k, 0]) / np.linalg.norm(va[k]),
        -1, 1))
    phi_teo = np.arctan2(2 * np.pi * f * 10e-3, 10.0)
    assert abs(phi - phi_teo) < 0.06, (phi, phi_teo)


# ------------------------------------------------------------- Transformador

def _arma_transf_trifasico():
    m = Modelo(dt=1e-4)
    ft = m.add(FuenteTrifasica("ft", amplitud=311.0, frecuencia=50.0))
    tr = m.add(Transformador("tr", a=2.0, fases=3))
    c = m.add(CargaRLTrifasica("c", R=10.0, L=10e-3))
    m.conectar(ft.salida, Puerto(tr, "ent", 0, 3))
    m.conectar(Puerto(tr, "sal", 0, 3), c.entrada)   # v2 -> carga
    m.conectar(c.salida, Puerto(tr, "ent", 3, 3))     # i2 -> primario
    return m, tr, c


def test_transformador_c_numpy_mono():
    m = Modelo(dt=1e-4)
    f1 = m.add(FuenteSeno("f1", amplitud=220.0, frecuencia=50.0))
    cero = m.add(FuenteConstante("cero", 0.0))
    tr = m.add(Transformador("tr", a=0.5))
    m.conectar(f1.salida, Puerto(tr, "ent", 0, 1))
    m.conectar(cero.salida, Puerto(tr, "ent", 1, 1))
    _compara(m, 0.1, [tr.salida])


def test_transformador_c_numpy_trifasico():
    m, tr, c = _arma_transf_trifasico()
    _compara(m, 0.5, [tr.salida, c.salida])


def test_transformador_fisica():
    # v2 = a*v1 ; i1 = -a*i2 (potencia conservada) ; I2 = V2/|Z|
    m, tr, c = _arma_transf_trifasico()
    res = m.run(0.5, registrar=[tr.salida, c.salida])
    t, out = res.t, res["tr"]
    k = t > 0.3
    va2 = np.max(np.abs(out[k, 0]))
    ia1 = np.max(np.abs(out[k, 3]))
    ia2 = np.max(np.abs(res["c"][k, 0]))
    assert abs(va2 - 2.0 * 311.0) / (2.0 * 311.0) < 0.01, va2
    Z = np.sqrt(10.0 ** 2 + (2 * np.pi * 50.0 * 10e-3) ** 2)
    assert abs(ia2 - va2 / Z) / (va2 / Z) < 0.02, ia2
    assert abs(ia1 - 2.0 * ia2) / (2.0 * ia2) < 0.02, (ia1, ia2)


def test_transformador_a_invalido():
    from bloques_crysi import Transformador
    with pytest.raises(ValueError):
        Transformador("tr", a=0.0)
    with pytest.raises(ValueError):
        Transformador("tr", a=-1.0)


# --------------------------------------------------------- Medidor de potencia

def _arma_medidor_maquina(tl=5.0):
    m = Modelo(dt=1e-5)
    red = m.add(FuenteTrifasica("red", amplitud=310.0, frecuencia=50.0))
    carga = m.add(FuenteConstante("tl", tl))
    maq = m.add(MaquinaImanesPermanentes(
        "pmac", rs=0.1, Ld=1e-3, Lq=1e-3, lam_m=0.1, P=6, J=0.01, Bm=0.001))
    med = m.add(MedidorPotencia("med"))
    m.conectar(red.salida, maq.terminales)
    m.conectar(carga.salida, maq.T_L)
    m.conectar(red.salida, med.entrada)          # tensiones de terminales
    m.conectar(maq.sensor3I(), med.corrientes)
    m.conectar(maq.sensorPar(), Puerto(med, "ent", 6, 1))
    m.conectar(maq.sensorVelocidad(), Puerto(med, "ent", 7, 1))
    return m, med, maq, red


def test_medidor_potencia_c_numpy():
    # fuentes constantes: P_e = va*ia + vb*ib + vc*ic ; P_m = Te*wm
    m = Modelo(dt=1e-3)
    vals = [100.0, -50.0, -50.0, 2.0, -1.0, -1.0, 3.0, 10.0]
    med = m.add(MedidorPotencia("med"))
    for k, v in enumerate(vals):
        f = m.add(FuenteConstante(f"f{k}", v))
        m.conectar(f.salida, Puerto(med, "ent", k, 1))
    arr, _ = _compara(m, 0.05, [med.salida])
    # P_e = 100*2 + (-50)(-1) + (-50)(-1) = 300 ;
    # Q_e = (vab*ic + vbc*ia + vca*ib)/sqrt(3) = (150*(-1) + 0 + (-150)*(-1))/sqrt(3) = 0
    # P_m = 3*10 = 30
    np.testing.assert_allclose(arr[:, -1], [300.0, 0.0, 30.0], rtol=1e-9, atol=1e-9)


def test_medidor_potencia_q_c_numpy():
    # caso con Q != 0: Q_e = (vab*ic + vbc*ia + vca*ib)/sqrt(3)
    m = Modelo(dt=1e-3)
    vals = [100.0, -50.0, -50.0, 0.0, 2.0, -1.0, 3.0, 10.0]
    med = m.add(MedidorPotencia("med"))
    for k, v in enumerate(vals):
        f = m.add(FuenteConstante(f"f{k}", v))
        m.conectar(f.salida, Puerto(med, "ent", k, 1))
    arr, _ = _compara(m, 0.05, [med.salida])
    # P_e = 100*0 + (-50)*2 + (-50)(-1) = -50
    # Q_e = (150*(-1) + 0*0 + (-150)*2)/sqrt(3) = -450/sqrt(3)
    # P_m = 3*10 = 30
    np.testing.assert_allclose(arr[:, -1], [-50.0, -450.0 / np.sqrt(3), 30.0],
                               rtol=1e-9, atol=1e-9)


def test_medidor_potencia_q_fisica():
    # carga RL trifasica: P = 3*Vrms^2*R/|Z|^2 ; Q = 3*Vrms^2*X/|Z|^2
    A, f, R, L = 311.0, 50.0, 10.0, 10e-3
    m = Modelo(dt=1e-5)
    red = m.add(FuenteTrifasica("red", amplitud=A, frecuencia=f))
    c = m.add(CargaRLTrifasica("c", R=R, L=L))
    med = m.add(MedidorPotencia("med"))
    m.conectar(red.salida, c.entrada)
    m.conectar(red.salida, med.entrada)
    m.conectar(c.salida, med.corrientes)
    res = m.run(0.2, registrar=[med.salida])
    k = res.t > 0.1
    Vrms = A / np.sqrt(2)
    Z = np.sqrt(R ** 2 + (2 * np.pi * f * L) ** 2)
    P = 3 * Vrms ** 2 * R / Z ** 2
    Q = 3 * Vrms ** 2 * (2 * np.pi * f * L) / Z ** 2
    assert abs(np.mean(res["med"][k, 0]) - P) / P < 0.01
    assert abs(np.mean(res["med"][k, 1]) - Q) / Q < 0.01
    # Q ~ 0 con carga casi resistiva (X = w*L << R)
    m2 = Modelo(dt=1e-5)
    red2 = m2.add(FuenteTrifasica("red2", amplitud=A, frecuencia=f))
    c2 = m2.add(CargaRLTrifasica("c2", R=R, L=1e-4))
    med2 = m2.add(MedidorPotencia("med2"))
    m2.conectar(red2.salida, c2.entrada)
    m2.conectar(red2.salida, med2.entrada)
    m2.conectar(c2.salida, med2.corrientes)
    res2 = m2.run(0.2, registrar=[med2.salida])
    k2 = res2.t > 0.1
    assert abs(np.mean(res2["med2"][k2, 1])) / np.mean(res2["med2"][k2, 0]) < 0.01


def test_medidor_potencia_maquina_c_numpy():
    # PMAC + red con par de carga: C == numpy con el medidor conectado
    m, med, _, _ = _arma_medidor_maquina()
    _compara(m, 0.2, [med.salida], atol=1e-8)


def test_medidor_potencia_fisica():
    # P_e = va*ia + vb*ib + vc*ic y P_m = Te*wm exactamente por construccion
    m, med, maq, red = _arma_medidor_maquina()
    res = m.run(0.2, registrar=[med.salida, red.salida, maq.sensor3I(),
                                maq.sensorPar(), maq.sensorVelocidad()])
    k = res.t > 0.1
    Pe = np.sum(res["red"][k] * res["I"][k], axis=1)
    np.testing.assert_allclose(res["med"][k, 0], Pe, rtol=1e-9, atol=1e-9)
    Pm = res["Te"][k] * res["wm"][k]
    np.testing.assert_allclose(res["med"][k, 2], Pm, rtol=1e-9, atol=1e-9)
    # en regimen el motor absorbe mas de lo que entrega al eje (perdidas)
    assert np.mean(Pm) > 0.0
    assert np.mean(res["med"][k, 0]) > np.mean(Pm)


def test_medidor_potencia_dc_c_numpy():
    # fases=1: P_e = V_t*ia ; P_m = Te*wm (generador: ambas negativas)
    m = Modelo(dt=1e-3)
    v = m.add(FuenteConstante("v", 12.0))
    i = m.add(FuenteConstante("i", -3.0))
    te = m.add(FuenteConstante("te", -2.0))
    wm = m.add(FuenteConstante("wm", 100.0))
    med = m.add(MedidorPotencia("med", fases=1))
    m.conectar(v.salida, med.entrada)
    m.conectar(i.salida, med.corrientes)
    m.conectar(te.salida, Puerto(med, "ent", 2, 1))
    m.conectar(wm.salida, Puerto(med, "ent", 3, 1))
    arr, _ = _compara(m, 0.05, [med.salida])
    # P_e = 12*(-3) = -36 ; P_m = -2*100 = -200
    np.testing.assert_allclose(arr[:, -1], [-36.0, -200.0], rtol=1e-9, atol=1e-9)


def test_medidor_potencia_cc_c_numpy():
    # CC como motor: C == numpy y en regimen P_e > 0, P_m > 0
    m = Modelo(dt=1e-4)
    va = m.add(FuenteConstante("va", 50.0))
    vf = m.add(FuenteConstante("vf", 5.0))
    tl = m.add(FuenteConstante("tl", 2.0))
    cc = m.add(MaquinaCorrienteContinua(
        "cc", r_a=0.5, L_a=0.01, r_f=1.0, L_f=0.1, L_AF=0.5, J=0.02, Bm=0.001))
    med = m.add(MedidorPotencia("med", fases=1))
    m.conectar(va.salida, cc.entrada)
    m.conectar(vf.salida, cc.campo)
    m.conectar(tl.salida, cc.T_L)
    m.conectar(cc.sensorVoltajeTerminal(), med.entrada)
    m.conectar(cc.sensorCorriente(), med.corrientes)
    m.conectar(cc.sensorPar(), Puerto(med, "ent", 2, 1))
    m.conectar(cc.sensorVelocidad(), Puerto(med, "ent", 3, 1))
    arr, _ = _compara(m, 0.3, [med.salida], atol=1e-8)
    Pe = np.mean(arr[0, -50:])
    Pm = np.mean(arr[1, -50:])
    assert Pe > 0.0 and Pm > 0.0, (Pe, Pm)


def test_medidor_potencia_fases_invalido():
    with pytest.raises(ValueError):
        MedidorPotencia("med", fases=2)