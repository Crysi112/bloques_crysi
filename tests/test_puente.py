"""PuenteInversorTrifasico/Monofasico: etapa de potencia para FOC."""

import numpy as np

from bloques_crysi import (
    FuenteConstante, FuenteTrifasica, Ganancia, GeneradorSVPWM, InvClarke,
    InvPark, MaquinaImanesPermanentes, Modelo, PID, Puerto,
    PuenteInversorMonofasico, PuenteInversorTrifasico, Sensor, Suma,
)
from bloques_crysi.backend_numpy import simular
from .helpers import resolver_registro


def par(bloque, tipo, offset, n):
    return Puerto(bloque, tipo, offset, n)


def _arma_puente_3f(promediado=False, Vdc=600.0, v_s=(1.0, 0.0, 0.0)):
    m = Modelo(dt=1e-5)
    vdc = m.add(FuenteConstante("vdc", Vdc))
    g0 = m.add(FuenteConstante("g0", 0.0))
    s = [m.add(FuenteConstante(f"s{k}", v)) for k, v in enumerate(v_s)]
    pu = m.add(PuenteInversorTrifasico("pu", promediado=promediado))
    m.conectar(vdc.salida, par(pu, "ent", 0, 1))
    for k in range(3):
        m.conectar(s[k].salida, par(pu, "ent", 1 + k, 1))
    return m, pu


def test_puente_trifasico_conmutado_estados():
    casos = {
        (1, 0, 0): (200.0, -100.0, -100.0),
        (1, 1, 0): (100.0, 100.0, -200.0),
        (0, 1, 1): (-200.0, 100.0, 100.0),
        (1, 1, 1): (0.0, 0.0, 0.0),
        (0, 0, 0): (0.0, 0.0, 0.0),
    }
    for (sa, sb, sc), (va, vb, vc) in casos.items():
        m, pu = _arma_puente_3f(Vdc=300.0, v_s=(sa, sb, sc))
        res = m.run(1e-3, registrar=[pu])
        np.testing.assert_allclose(res["pu"][-1], (va, vb, vc),
                                   atol=1e-9)


def test_puente_trifasico_promediado():
    casos = {
        (1.0, -1.0, -1.0): (200.0, -100.0, -100.0),
        (0.5, 0.5, 0.5): (0.0, 0.0, 0.0),
        (0.5, -0.5, 0.0): (75.0, -75.0, 0.0),
        (0.0, 0.0, 0.0): (0.0, 0.0, 0.0),
    }
    for (ma, mb, mc), (va, vb, vc) in casos.items():
        m, pu = _arma_puente_3f(promediado=True, Vdc=300.0,
                                v_s=(ma, mb, mc))
        res = m.run(1e-3, registrar=[pu])
        np.testing.assert_allclose(res["pu"][-1], (va, vb, vc),
                                   atol=1e-9)


def test_puente_monofasico():
    m = Modelo(dt=1e-5)
    vdc = m.add(FuenteConstante("vdc", 48.0))
    sa = m.add(FuenteConstante("sa", 1.0))
    sb = m.add(FuenteConstante("sb", 0.0))
    pu = m.add(PuenteInversorMonofasico("pu"))
    m.conectar(vdc.salida, par(pu, "ent", 0, 1))
    m.conectar(sa.salida, par(pu, "ent", 1, 1))
    m.conectar(sb.salida, par(pu, "ent", 2, 1))
    res = m.run(1e-3, registrar=[pu])
    np.testing.assert_allclose(res["pu"][-1], 48.0, atol=1e-9)
    # promediado: (ma, mb) = (0.5, -0.5) -> Vout = 48*1/2 = 24
    m2 = Modelo(dt=1e-5)
    vdc2 = m2.add(FuenteConstante("vdc2", 48.0))
    ma = m2.add(FuenteConstante("ma", 0.5))
    mb = m2.add(FuenteConstante("mb", -0.5))
    pu2 = m2.add(PuenteInversorMonofasico("pu2", promediado=True))
    m2.conectar(vdc2.salida, par(pu2, "ent", 0, 1))
    m2.conectar(ma.salida, par(pu2, "ent", 1, 1))
    m2.conectar(mb.salida, par(pu2, "ent", 2, 1))
    res2 = m2.run(1e-3, registrar=[pu2])
    np.testing.assert_allclose(res2["pu2"][-1], 24.0, atol=1e-9)


def test_puente_trifasico_c_igual_numpy():
    m, pu = _arma_puente_3f(promediado=True, Vdc=600.0, v_s=(0.4, -0.2, 0.1))
    res = m.run(1e-3, registrar=[pu])
    m._resolver()
    rec_idx, _ = resolver_registro(m, [pu])
    datos = simular(m.bloques, m.dt, 1e-3, rec_idx,
                    max_iter=m.max_iter, tol=m.tol, w_opt=m.w_opt,
                    orden_estatico=m._orden_estatico())
    np.testing.assert_allclose(datos[0], res["pu"][:, 0], rtol=1e-12, atol=1e-9)
    np.testing.assert_allclose(datos[1], res["pu"][:, 1], rtol=1e-12, atol=1e-9)
    np.testing.assert_allclose(datos[2], res["pu"][:, 2], rtol=1e-12, atol=1e-9)


def _arma_foc(promediado=True):
    """FOC completo: PID de velocidad -> PID de corriente -> InvPark ->
    InvClarke -> escala 2/Vdc -> puente (promediado o conmutado) -> PMAC.
    Vdc = 600 V, w_ref = 100 rad/s.
    """
    Vdc, w_ref = 600.0, 100.0
    m = Modelo(dt=1e-5)
    ref_w = m.add(FuenteConstante("ref_w", w_ref))
    id_ref = m.add(FuenteConstante("id_ref", 0.0))
    vdc = m.add(FuenteConstante("vdc", Vdc))
    tl = m.add(FuenteConstante("tl", 0.0))
    maq = m.add(MaquinaImanesPermanentes(
        "pmsm", rs=1.0, Ld=1e-3, Lq=1e-3, lam_m=0.1, P=6, J=0.01,
        Bm=0.001))
    m.conectar(tl.salida, maq.T_L)
    wm = maq.sensorVelocidad()
    iq = Sensor("iq", maq, "sal", 3, 1, canales=["iq"])
    id_ = Sensor("id", maq, "sal", 4, 1, canales=["id"])
    th = Sensor("th", maq, "sal", 7, 1, canales=["th"])
    err_w = m.add(Suma("err_w", (1.0, -1.0)))
    pid_w = m.add(PID("pid_w", Kp=0.2, Ki=8.0, u_min=-30.0, u_max=30.0))
    err_q = m.add(Suma("err_q", (1.0, -1.0)))
    pid_q = m.add(PID("pid_q", Kp=10.0, Ki=100.0,
                      u_min=-Vdc, u_max=Vdc))
    err_d = m.add(Suma("err_d", (1.0, -1.0)))
    pid_d = m.add(PID("pid_d", Kp=10.0, Ki=100.0,
                      u_min=-Vdc, u_max=Vdc))
    ipk = m.add(InvPark("ipk"))
    icl = m.add(InvClarke("icl"))
    k2v = [m.add(Ganancia(f"k2v{k}", 2.0 / Vdc)) for k in range(3)]
    pu = m.add(PuenteInversorTrifasico("pu", promediado=promediado))
    m.conectar(ref_w.salida, par(err_w, "ent", 0, 1))
    m.conectar(wm, par(err_w, "ent", 1, 1))
    m.conectar(err_w.salida, pid_w.entrada)
    m.conectar(pid_w.salida, par(err_q, "ent", 0, 1))
    m.conectar(iq, par(err_q, "ent", 1, 1))
    m.conectar(err_q.salida, pid_q.entrada)
    m.conectar(id_ref.salida, par(err_d, "ent", 0, 1))
    m.conectar(id_, par(err_d, "ent", 1, 1))
    m.conectar(err_d.salida, pid_d.entrada)
    m.conectar(pid_d.salida, par(ipk, "ent", 0, 1))
    m.conectar(pid_q.salida, par(ipk, "ent", 1, 1))
    m.conectar(th, par(ipk, "ent", 2, 1))
    m.conectar(par(ipk, "sal", 0, 1), par(icl, "ent", 0, 1))
    m.conectar(par(ipk, "sal", 1, 1), par(icl, "ent", 1, 1))
    m.conectar(vdc.salida, par(pu, "ent", 0, 1))
    for k in range(3):
        m.conectar(par(icl, "sal", k, 1), k2v[k].entrada)
        m.conectar(k2v[k].salida, par(pu, "ent", 1 + k, 1))
    m.conectar(pu.salida, maq.terminales)
    return m, maq, pu, iq, id_


def test_foc_promediado_alcanza_velocidad():
    m, maq, pu, iq, id_ = _arma_foc(promediado=True)
    res = m.run(1.5, registrar=[maq.sensorVelocidad(), iq, id_])
    wm = res["wm"]
    assert abs(wm[-1] - 100.0) / 100.0 < 0.03
    assert np.all(np.isfinite(wm))
    assert abs(np.mean(res["iq"][-200:])) > 0.0   # hay corriente de par
    # el puente promediado reproduce las tensiones de la planta ideal:
    # corrientes acotadas en regimen
    assert np.max(np.abs(res["iq"])) < 40.0


def test_foc_conmutado_funciona_y_es_finito():
    m, maq, pu, iq, id_ = _arma_foc(promediado=False)
    m2, maq2, pu2, iq2, id2 = _arma_foc(promediado=True)
    res = m.run(0.05, registrar=[maq.sensorVelocidad(), iq])
    res2 = m2.run(0.05, registrar=[maq2.sensorVelocidad(), iq2])
    assert np.all(np.isfinite(res["iq"]))
    assert np.all(np.isfinite(res["wm"]))
    assert res["wm"][-1] > 20.0
    np.testing.assert_allclose(np.mean(res["wm"][-2000:]),
                               np.mean(res2["wm"][-2000:]), rtol=0.5)


def test_foc_con_generador_svpwm():
    """Cadena industrial completa: PID -> InvPark -> SVPWM (disparos
    reales) -> Puente conmutado -> PMAC."""
    Vdc, w_ref = 600.0, 100.0
    m = Modelo(dt=1e-5)
    ref_w = m.add(FuenteConstante("ref_w", w_ref))
    vdc = m.add(FuenteConstante("vdc", Vdc))
    tl = m.add(FuenteConstante("tl", 0.0))
    maq = m.add(MaquinaImanesPermanentes(
        "pmsm", rs=1.0, Ld=1e-3, Lq=1e-3, lam_m=0.1, P=6, J=0.01,
        Bm=0.001))
    m.conectar(tl.salida, maq.T_L)
    wm = maq.sensorVelocidad()
    iq = Sensor("iq", maq, "sal", 3, 1, canales=["iq"])
    th = Sensor("th", maq, "sal", 7, 1, canales=["th"])
    err_w = m.add(Suma("err_w", (1.0, -1.0)))
    pid_w = m.add(PID("pid_w", Kp=0.2, Ki=8.0, u_min=-30.0, u_max=30.0))
    err_q = m.add(Suma("err_q", (1.0, -1.0)))
    pid_q = m.add(PID("pid_q", Kp=10.0, Ki=100.0, u_min=-Vdc, u_max=Vdc))
    ipk = m.add(InvPark("ipk"))
    svpwm = m.add(GeneradorSVPWM("svpwm", Vdc=Vdc, fsw=5000.0))
    pu = m.add(PuenteInversorTrifasico("pu"))
    m.conectar(ref_w.salida, par(err_w, "ent", 0, 1))
    m.conectar(wm, par(err_w, "ent", 1, 1))
    m.conectar(err_w.salida, pid_w.entrada)
    m.conectar(pid_w.salida, par(err_q, "ent", 0, 1))
    m.conectar(iq, par(err_q, "ent", 1, 1))
    m.conectar(err_q.salida, pid_q.entrada)
    m.conectar(pid_q.salida, par(ipk, "ent", 0, 1))
    m.conectar(par(maq, "sal", 4, 1), par(ipk, "ent", 1, 1))  # id -> 0
    m.conectar(th, par(ipk, "ent", 2, 1))
    m.conectar(par(ipk, "sal", 0, 1), svpwm.v_alpha)
    m.conectar(par(ipk, "sal", 1, 1), svpwm.v_beta)
    m.conectar(vdc.salida, par(pu, "ent", 0, 1))
    m.conectar(svpwm.salida, par(pu, "ent", 1, 3))
    m.conectar(pu.salida, maq.terminales)
    res = m.run(0.1, registrar=[maq.sensorVelocidad(), iq])
    assert np.all(np.isfinite(res["iq"]))
    assert res["wm"][-1] > 60.0   # avanza hacia los 100 rad/s