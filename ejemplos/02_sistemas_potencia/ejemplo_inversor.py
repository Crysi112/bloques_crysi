from bloques_crysi import (
    Modelo, FuenteConstante, InversorTrifasico, Scope, FiltroPasoBajo,
)
from bloques_crysi import Puerto

DT = 1e-6
T_END = 600e-3

Vdc, fsw, f_out = 400.0, 10000.0, 60.0
m_start, m_end, t_ramp = 0.4, 1.0, 1e-6
Lf, Cf, Rload = 5e-3, 20e-6, 100.0

def _resumen(res, medicion, conmutada):
    vCa = res["scope"][:, 0]
    iLa = res["scope"][:, 3]
    vCa_f = res["medicion"] if medicion else None
    ns = int(1.2 / f_out / DT)
    amp = (max(vCa[-ns:]) - min(vCa[-ns:])) / 2
    amp_f = (max(vCa_f[-ns:]) - min(vCa_f[-ns:])) / 2 if medicion else 0.0
    print("Inversor trifasico SPWM con filtro LC")
    print("=" * 58)
    print(f"modulacion        : {'conmutada (PWM)' if conmutada else 'promediada'}")
    print(f"amplitud vCa final: {amp:9.3f} V  (teorico m*Vdc/2 = "
          f"{m_end*Vdc/2:.0f} V)")
    if medicion:
        print(f"amplitud vCa filt : {amp_f:9.3f} V  (medicion con "
              f"FiltroPasoBajo fc=500 Hz)")
    print(f"corriente pico iLa: {max(abs(float(x)) for x in iLa[-800:]):9.3f} A")
    print("nota: las tensiones llevan la componente DC de Vdc/2 (referencia")
    print("      al borne negativo del bus; fase-neutro queda centrada en cero)")

def _arma_inversor(conmutada):
    m = Modelo(dt=DT)
    bus = m.add(FuenteConstante("vdc", Vdc))
    inv = m.add(InversorTrifasico(
        "inv", f_out=f_out, fsw=fsw,
        m_start=m_start, m_end=m_end, t_ramp=t_ramp,
        Lf=Lf, Cf=Cf, R=Rload, conmutada=conmutada))
    m.conectar(bus.salida, inv.entrada)
    sc = m.add(Scope("scope", anchos=[3, 1, 2, 1],
                     guiones=["vCa (V)", "vCb (V)", "vCc (V)",
                              "iLa (A)", "iLb (A)", "iLc (A)"]))
    m.conectar(Puerto(inv, "sal", 0, 3), sc.canales[0])
    m.conectar(Puerto(inv, "sal", 3, 1), sc.canales[1])
    m.conectar(Puerto(inv, "sal", 4, 2), sc.canales[2])
    if conmutada:
        medicion = m.add(FiltroPasoBajo("medicion", fc=500.0, orden=1))
        m.conectar(Puerto(inv, "sal", 0, 1), medicion.entrada)
        m.conectar(medicion.salida, sc.canales[3])
    else:
        medicion = None
    return m, sc, medicion

m, sc, medicion = _arma_inversor(conmutada=False)
res = m.run(t_fin=T_END, registrar=[sc])
_resumen(res, medicion, conmutada=False)

m, sc, medicion = _arma_inversor(conmutada=True)
res = m.run(t_fin=T_END, registrar=[sc, medicion])
_resumen(res, medicion, conmutada=True)
