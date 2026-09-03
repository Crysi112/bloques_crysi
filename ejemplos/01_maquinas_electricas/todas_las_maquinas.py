import numpy as np
import matplotlib.pyplot as plt
from bloques_crysi import (Modelo, FuenteTrifasica, FuenteConstante,
                           MaquinaInduccion, MaquinaSincrona,
                           MaquinaImanesPermanentes, MaquinaCorrienteContinua)

def ver(wm, te, t, titulo, ws=None):
    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(9, 5))
    axs[0].plot(t, wm, lw=1.2)
    if ws is not None:
        axs[0].axhline(ws, color='r', ls='--', lw=1, label=f'ws = {ws:.1f}')
    axs[0].set_ylabel('velocidad (rad/s)')
    axs[0].legend(); axs[0].grid(alpha=0.3)
    axs[1].plot(t, te, lw=1.2)
    axs[1].axhline(0, color='k', lw=0.8)
    axs[1].set_ylabel('Te (N.m)')
    axs[1].set_xlabel('t (s)')
    axs[1].grid(alpha=0.3)
    fig.suptitle(titulo)
    plt.tight_layout()
    plt.show()

def arma_mi(tl):
    m = Modelo(dt=1e-4)
    red = m.add(FuenteTrifasica('red', 310.0, 50.0))
    t = m.add(FuenteConstante('tl', tl))
    maq = m.add(MaquinaInduccion(
        'mi', rs=0.5, rr=0.4, Lm=0.1, Lls=0.005, Llr=0.005,
        P=4, J=0.5, Bm=0.01))
    m.conectar(red.salida, maq.terminales)
    m.conectar(t.salida, maq.T_L)
    return m, maq

m, maq = arma_mi(0.0)
res = m.run(4.0, registrar=[maq.sensorVelocidad(), maq.sensorPar(),
                            maq.sensor3I()])
ws = 2*np.pi*50/2
ver(res['wm'], res['Te'], res.t, 'MI en vacío (T_L = 0)', ws)
print(f'wm final = {res["wm"][-1]:.2f} rad/s  (ws = {ws:.2f})')
print(f'Te medio = {np.mean(res["Te"][-1000:]):.2f} N.m  '
      f'I pico = {np.abs(res["I"]).max():.1f} A')

for tl, nombre in ((30.0, 'motor'), (-30.0, 'generador')):
    m, maq = arma_mi(tl)
    res = m.run(4.0, registrar=[maq.sensorVelocidad(), maq.sensorPar()])
    ver(res['wm'], res['Te'], res.t, f'MI modo {nombre} (T_L = {tl:+.0f})', ws)
    print(f'T_L = {tl:+.0f} | wm final = {res["wm"][-1]:6.2f} rad/s | '
          f'Te medio = {np.mean(res["Te"][-1000:]):+7.2f} N.m')

def arma_ms(tl):
    m = Modelo(dt=1e-4)
    red = m.add(FuenteTrifasica('red', 310.0, 50.0))
    vfd = m.add(FuenteConstante('vfd', 2.0))
    t = m.add(FuenteConstante('tl', tl))
    maq = m.add(MaquinaSincrona(
        'ms', rs=0.3, rfd=0.5, rkq1=0.1, rkq2=0.1, rkd=0.1,
        Lls=0.002, Lmq=0.08, Llkq1=0.01, Llkq2=0.005,
        Lmd=0.1, Llf=0.02, Llkd=0.005, P=4, J=1.0, Bm=0.01))
    m.conectar(red.salida, maq.terminales)
    m.conectar(vfd.salida, maq.vfd)
    m.conectar(t.salida, maq.T_L)
    return m, maq

ws = 2*np.pi*50/2
m, maq = arma_ms(0.0)
res = m.run(6.0, registrar=[maq.sensorVelocidad(), maq.sensorPar()])
ver(res['wm'], res['Te'], res.t, 'MS en vacío (T_L = 0)', ws)
print(f'wm final = {res["wm"][-1]:.2f} rad/s  (ws = {ws:.2f})  '
      f'Te medio = {np.mean(res["Te"][-2000:]):.2f} N.m')

for tl, nombre in ((5.0, 'motor'), (-20.0, 'generador')):
    m, maq = arma_ms(tl)
    res = m.run(6.0, registrar=[maq.sensorVelocidad(), maq.sensorPar()])
    ver(res['wm'], res['Te'], res.t, f'MS modo {nombre} (T_L = {tl:+.0f})', ws)
    print(f'T_L = {tl:+.0f} | wm final = {res["wm"][-1]:6.2f} rad/s | '
          f'Te medio = {np.mean(res["Te"][-2000:]):+7.2f} N.m')

def arma_pmac(tl):
    m = Modelo(dt=1e-4)
    red = m.add(FuenteTrifasica('red', 310.0, 50.0))
    t = m.add(FuenteConstante('tl', tl))
    maq = m.add(MaquinaImanesPermanentes(
        'pmsm', rs=0.1, Ld=1e-3, Lq=1e-3, lam_m=0.5, P=6, J=0.05, Bm=0.5))
    m.conectar(red.salida, maq.terminales)
    m.conectar(t.salida, maq.T_L)
    return m, maq

ws = 2*np.pi*50/3
m, maq = arma_pmac(0.0)
res = m.run(3.0, registrar=[maq.sensorVelocidad(), maq.sensorPar()])
ver(res['wm'], res['Te'], res.t, 'PMAC en vacío (T_L = 0)', ws)
print(f'wm final = {res["wm"][-1]:.2f} rad/s  (ws = {ws:.2f})  '
      f'Te medio = {np.mean(res["Te"][-1000:]):.2f} N.m')

for tl, nombre in ((40.0, 'motor'), (-60.0, 'generador')):
    m, maq = arma_pmac(tl)
    res = m.run(3.0, registrar=[maq.sensorVelocidad(), maq.sensorPar()])
    ver(res['wm'], res['Te'], res.t, f'PMAC modo {nombre} (T_L = {tl:+.0f})', ws)
    print(f'T_L = {tl:+.0f} | wm final = {res["wm"][-1]:6.2f} rad/s | '
          f'Te medio = {np.mean(res["Te"][-1000:]):+7.2f} N.m')

def arma_cc(tl):
    m = Modelo(dt=1e-5)
    va = m.add(FuenteConstante('va', 120.0))
    vf = m.add(FuenteConstante('vf', 100.0))
    t = m.add(FuenteConstante('tl', tl))
    maq = m.add(MaquinaCorrienteContinua(
        'cc', r_a=1.0, L_a=0.01, r_f=100.0, L_f=10.0, L_AF=1.5, J=0.5))
    m.conectar(va.salida, maq.entrada)
    m.conectar(vf.salida, maq.campo)
    m.conectar(t.salida, maq.T_L)
    return m, maq

m, maq = arma_cc(0.0)
res = m.run(2.0, registrar=[maq.sensorVelocidad(), maq.sensorPar(),
                            maq.sensorCorriente(), maq.sensorCampo()])
ver(res['wm'], res['Te'], res.t, 'CC en vacío (T_L = 0)')
print(f'wm final = {res["wm"][-1]:.2f} rad/s  (vacío, va=120 → 80)')
print(f'ia final = {res["ia"][-1]:+.2f} A   if = {res["if"][-1]:.2f} A')

for tl, nombre in ((10.0, 'motor'), (-10.0, 'generador')):
    m, maq = arma_cc(tl)
    res = m.run(2.0, registrar=[maq.sensorVelocidad(), maq.sensorPar(),
                                maq.sensorCorriente()])
    ver(res['wm'], res['Te'], res.t, f'CC modo {nombre} (T_L = {tl:+.0f})')
    print(f'T_L = {tl:+.0f} | wm final = {res["wm"][-1]:6.2f} rad/s | '
          f'ia final = {res["ia"][-1]:+7.2f} A | '
          f'Te final = {res["Te"][-1]:+7.2f} N.m')
