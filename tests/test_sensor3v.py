"""sensor3V(): medir las tensiones aplicadas a los terminales de una maquina."""

import numpy as np

from bloques_crysi import (FuenteConstante, FuenteTrifasica, Modelo, Puerto,
                           Scope, Suma)
from bloques_crysi import MaquinaImanesPermanentes


def _arma():
    m = Modelo(dt=1e-4)
    f = m.add(FuenteTrifasica("f", 100.0, 50.0))
    maq = m.add(MaquinaImanesPermanentes(
        "pmsm", rs=1.0, Ld=1e-3, Lq=1e-3, lam_m=0.1, P=6, J=0.01))
    tl = m.add(FuenteConstante("tl", 0.0))
    m.conectar(tl.salida, maq.T_L)
    m.conectar(f.salida, maq.terminales)
    return m, maq, f


def test_sensor3V_conectable_a_suma():
    # la tension medida en los terminales es exactamente la de la fuente:
    # va - V(medida) = 0 en cada fase
    m, maq, f = _arma()
    s3v = maq.sensor3V()
    sumas = [m.add(Suma(f"sum{k}", signos=[1.0, -1.0])) for k in range(3)]
    for k in range(3):
        m.conectar(Puerto(f, "sal", k, 1), Puerto(sumas[k], "ent", 0, 1))
        m.conectar(Puerto(s3v.bloque, s3v.tipo, s3v.offset + k, 1),
                   Puerto(sumas[k], "ent", 1, 1))
    res = m.run(0.05, registrar=[sumas[0].salida, sumas[1].salida,
                                 sumas[2].salida])
    for k in range(3):
        np.testing.assert_allclose(res[f"sum{k}"], 0.0, atol=1e-6)


def test_sensor3V_registro():
    # registrar el sensor directamente: 3 canales con las tensiones
    m, maq, f = _arma()
    res = m.run(0.05, registrar=[maq.sensor3V()])
    assert res["V"].shape == (501, 3)
    assert res["V"][:, 0].max() > 90.0   # va oscila entre -100 y 100


def test_sensor3V_en_scope():
    m, maq, f = _arma()
    sc = m.add(Scope("s", anchos=[3], mostrar=False))
    m.conectar(maq.sensor3V(), sc.canales[0])
    assert sc.canales_meta[:3] == ["va", "vb", "vc"]
    m.run(0.05, registrar=[sc])
    assert len(sc.in_idx) == 3 and all(k >= 0 for k in sc.in_idx)