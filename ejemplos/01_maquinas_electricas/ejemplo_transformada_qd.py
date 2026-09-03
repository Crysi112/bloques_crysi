import numpy as np

from bloques_crysi import (
    Modelo, FuenteTrifasica, FuenteConstante, FuenteRampa,
    TransformadaQD, Scope, Puerto, MaquinaImanesPermanentes,
)

DT = 1e-4

def par(bloque, tipo, offset, n):
    return Puerto(bloque, tipo, offset, n)

def _test_fuentes():
    m = Modelo(dt=DT)
    src = m.add(FuenteTrifasica("red", amplitud=310.0, frecuencia=50.0))
    ramp = m.add(FuenteRampa("th", 2 * np.pi * 50.0, 0.0))
    qd = m.add(TransformadaQD("qd"))
    m.conectar(src.salida, qd.vabc)
    m.conectar(src.salida, qd.iabc)
    m.conectar(ramp.salida, qd.th)

    sc = m.add(Scope("qd_scope", anchos=[1, 1, 1, 1],
                     guiones=["vqs (V)", "vds (V)", "iqs (A)", "ids (A)"]))
    m.conectar(par(qd, "sal", 0, 1), sc.canales[0])
    m.conectar(par(qd, "sal", 1, 1), sc.canales[1])
    m.conectar(par(qd, "sal", 2, 1), sc.canales[2])
    m.conectar(par(qd, "sal", 3, 1), sc.canales[3])

    res = m.run(t_fin=0.1, registrar=[sc, src.salida, ramp.salida, qd])
    i = 500
    va, vb, vc = res["red"][i]
    ia, ib, ic = res["red"][i]
    th = float(res["th"][i])
    qd_vals = res["qd"][i]

    alv = (2.0 / 3.0) * (va - 0.5 * vb - 0.5 * vc)
    bev = (vb - vc) / np.sqrt(3.0)
    ali = (2.0 / 3.0) * (ia - 0.5 * ib - 0.5 * ic)
    bei = (ib - ic) / np.sqrt(3.0)
    dv = alv * np.cos(th) + bev * np.sin(th)
    qv = -alv * np.sin(th) + bev * np.cos(th)
    di = ali * np.cos(th) + bei * np.sin(th)
    qi = -ali * np.sin(th) + bei * np.cos(th)

    vqs_b, vds_b, iqs_b, ids_b = qd_vals[0], qd_vals[1], qd_vals[2], qd_vals[3]

    print("Test TransformadaQD con fuentes puras")
    print("=" * 50)
    print(f"t = {i*DT:.3f} s, th = {th:.4f} rad")
    print(f"vqs: bloque={vqs_b:.6f}  manual={qv:.6f}  diff={abs(vqs_b-qv):.2e}")
    print(f"vds: bloque={vds_b:.6f}  manual={dv:.6f}  diff={abs(vds_b-dv):.2e}")
    print(f"iqs: bloque={iqs_b:.6f}  manual={qi:.6f}  diff={abs(iqs_b-qi):.2e}")
    print(f"ids: bloque={ids_b:.6f}  manual={di:.6f}  diff={abs(ids_b-di):.2e}")

    tol = 1e-9
    assert abs(vqs_b - qv) < tol
    assert abs(vds_b - dv) < tol
    assert abs(iqs_b - qi) < tol
    assert abs(ids_b - di) < tol
    print("\nOK - TransformadaQD coincide exacto con cálculo manual")
    return res

res_fuentes = _test_fuentes()

def _test_maquina():
    m = Modelo(dt=DT)
    src = m.add(FuenteTrifasica("red", amplitud=310.0, frecuencia=50.0))
    tl = m.add(FuenteConstante("tl", 5.0))
    maq = m.add(MaquinaImanesPermanentes(
        "pmsm", rs=0.1, Ld=1e-3, Lq=1e-3, lam_m=0.1, P=6, J=0.01, Bm=0.001))
    qd = m.add(TransformadaQD("qd"))
    m.conectar(src.salida, maq.terminales)
    m.conectar(tl.salida, maq.T_L)
    m.conectar(maq.sensor3V(), qd.vabc)
    m.conectar(maq.sensor3I(), qd.iabc)
    m.conectar(maq.sensorPosicionElectrica(), qd.th)

    res = m.run(t_fin=0.2, registrar=[
        qd.salida, maq.sensorCorrienteQ(), maq.sensorCorrienteD()
    ])
    i = 1500
    iqs_b, ids_b = res["qd"][i][2], res["qd"][i][3]
    iqs_maq, ids_maq = res["iqs"][i], res["ids"][i]

    print("Test TransformadaQD con sensores de PMAC")
    print("=" * 50)
    print(f"iqs: bloque={iqs_b:.6f}  sensor_maq={iqs_maq:.6f}  diff={abs(iqs_b-iqs_maq):.2e}")
    print(f"ids: bloque={ids_b:.6f}  sensor_maq={ids_maq:.6f}  diff={abs(ids_b-ids_maq):.2e}")

    assert abs(iqs_b - iqs_maq) < 1e-6
    assert abs(ids_b - ids_maq) < 1e-6
    print("OK - TransformadaQD coincide con sensores internos de la máquina")
    return res

res_maquina = _test_maquina()
