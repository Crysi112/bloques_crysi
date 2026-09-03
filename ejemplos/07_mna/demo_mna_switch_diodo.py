import numpy as np

try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

from bloques_crysi.modelo import Modelo
from bloques_crysi.bloques import FuenteSeno, FuenteEscalon, PulsoRectangular, Multiplexor, FuenteConstante
from bloques_crysi.mna import Nodo, Resistor, Capacitor, Inductor, VSource, Diodo, Switch, SubredMNA

def test_rc():
    print("Ejecutando Test 1: Circuito RC simple...")
    m = Modelo(dt=1e-5)

    n_in = Nodo("in")
    n_out = Nodo("out")
    gnd = Nodo("0")

    v_src = VSource("Vin", n_in, gnd, idx_u=0)
    r1 = Resistor("R1", n_in, n_out, R=1000.0)
    c1 = Capacitor("C1", n_out, gnd, C=1e-6)

    mna = m.add(SubredMNA("mna_rc", 
        nodos=[n_in, n_out, gnd], 
        componentes=[v_src, r1, c1], 
        dt=m.dt, 
        mediciones_v=[(n_out, gnd)], 
        mediciones_i=[v_src],
        precomputar=True
    ))

    step = m.add(FuenteEscalon("step", valor_final=10.0, t_paso=0.001))
    m.conectar(step.salida, mna.entrada)

    res = m.run(t_fin=0.01, registrar=[mna])

    if not HAS_MATPLOTLIB:
        print("Prueba completada. (Matplotlib no instalado para graficar).")
        return

    t = res.t
    vout = res.get("mna_rc")
    if vout is not None and vout.ndim == 2:
        vout = vout[:, 0]

    plt.figure(figsize=(8, 4))
    plt.plot(t, vout, label="V_out (MNA)", linewidth=2)

    t_step = t[t >= 0.001] - 0.001
    v_analitico = np.zeros_like(t)
    v_analitico[t >= 0.001] = 10.0 * (1.0 - np.exp(-t_step / 0.001))

    plt.plot(t, v_analitico, '--', label="V_out (Analítico)", linewidth=2)
    plt.title("Circuito RC - Respuesta al Escalón")
    plt.xlabel("Tiempo [s]")
    plt.ylabel("Tensión [V]")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def test_rectificador():
    print("Ejecutando Test 2: Rectificador de media onda...")
    m = Modelo(dt=1e-5)

    n_in = Nodo("in")
    n_mid = Nodo("mid")
    n_out = Nodo("out")
    gnd = Nodo("0")

    v_src = VSource("Vac", n_in, gnd, idx_u=0)
    d1 = Diodo("D1", n_in, n_mid, Ron=1e-3, Roff=1e6, Vf=0.7)
    r_series = Resistor("Rs", n_mid, n_out, R=0.5)
    c1 = Capacitor("C_filt", n_out, gnd, C=1e-3)
    r_load = Resistor("R_load", n_out, gnd, R=50.0)

    mna = m.add(SubredMNA("mna_rect", 
        nodos=[n_in, n_mid, n_out, gnd], 
        componentes=[v_src, d1, r_series, c1, r_load], 
        dt=m.dt, 
        mediciones_v=[(n_in, gnd), (n_out, gnd)], 
        mediciones_i=[],
        precomputar=True
    ))

    ac = m.add(FuenteSeno("ac_source", amplitud=311, frecuencia=50))
    m.conectar(ac.salida, mna.entrada)

    res = m.run(t_fin=0.1, registrar=[mna])

    if not HAS_MATPLOTLIB:
        print("Prueba completada.")
        return

    plt.figure(figsize=(8, 4))
    v_in = res.get("mna_rect")
    v_out = res.get("mna_rect")
    if v_in is not None and v_in.ndim == 2:
        v_in = v_in[:, 0]
        v_out = v_out[:, 1]
    plt.plot(res.t, v_in, label="Vin (AC)", alpha=0.5)
    plt.plot(res.t, v_out, label="Vout (Filtrada)", linewidth=2)
    plt.title("Rectificador de Media Onda (Diodo MNA)")
    plt.xlabel("Tiempo [s]")
    plt.ylabel("Tensión [V]")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def test_boost():
    print("Ejecutando Test 3: Convertidor Boost...")
    m = Modelo(dt=1e-6)

    n_in = Nodo("in")
    n_sw = Nodo("sw")
    n_out = Nodo("out")
    gnd = Nodo("0")

    v_dc = VSource("Vdc", n_in, gnd, idx_u=0)
    l1 = Inductor("L1", n_in, n_sw, L=1e-3)

    sw = Switch("S1", n_sw, gnd, idx_ctrl=1, Ron=1e-3, Roff=1e6)
    d1 = Diodo("D1", n_sw, n_out, Ron=1e-3, Roff=1e6, Vf=0.7)
    c1 = Capacitor("C1", n_out, gnd, C=10e-6)
    r1 = Resistor("R_load", n_out, gnd, R=20.0)

    mna = m.add(SubredMNA("mna_boost", 
        nodos=[n_in, n_sw, n_out, gnd], 
        componentes=[v_dc, l1, sw, d1, c1, r1], 
        dt=m.dt, 
        mediciones_v=[(n_out, gnd)], 
        mediciones_i=[v_dc],
        precomputar=False
    ))

    vin = m.add(FuenteConstante("vin", valor=12.0))
    pwm = m.add(PulsoRectangular("pwm", amplitud=1.0, periodo=1/50000.0, duty=0.5))

    mux = m.add(Multiplexor("mux", n_canales=2))
    m.conectar(vin.salida, mux.entradas[0])
    m.conectar(pwm.salida, mux.entradas[1])

    m.conectar(mux.salida, mna.entrada)

    res = m.run(t_fin=0.01, registrar=[mna])

    if not HAS_MATPLOTLIB:
        print("Prueba completada.")
        return

    plt.figure(figsize=(8, 4))
    vout = res.get("mna_boost")
    if vout is not None and vout.ndim == 2:
        vout = vout[:, 0]
    plt.plot(res.t, vout, label="Vout (Boost)")

    plt.axhline(24.0 - 0.7, color='r', linestyle='--', label="Teórico (ideal - Vf)")

    plt.title("Convertidor Boost (MNA con Switch Controlado por Signal Flow)")
    plt.xlabel("Tiempo [s]")
    plt.ylabel("Tensión [V]")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    test_rc()
    test_rectificador()
    test_boost()
