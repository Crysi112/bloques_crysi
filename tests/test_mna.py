"""
Tests del bloque SubredMNA (Modified Nodal Analysis).
Circuito de prueba: RLC serie excitado por fuente de voltaje senoidal.
  vs(t) = Vm*sin(wt)  --- R --- L --- C ---
La respuesta en régimen permanente (ω²LC=1, resonancia) permite
verificar valores analíticos.
"""
import pytest
import numpy as np
from bloques_crysi import (Modelo, FuenteSeno, GeneradorCodigo,
                           SubredMNA, Nodo, Resistor, Capacitor, Inductor,
                           VSource, ISource, Switch)


def _hacer_circuito_rlc(R=10.0, L=1e-3, C=1e-6, dt=1e-6):
    """
    Crea un circuito RLC serie: vs(t) -> R -> L -> C -> GND.
    Nodos: N1 (entre VS y R), N2 (entre R y L), N3 (entre L y C).
    GND = nodo de referencia.
    Mide: V_C = voltaje en N3 (entre capacitor y GND).
    Entrada 0: vs(t)
    """
    gnd = Nodo("GND")   # Referencia
    n1 = Nodo("N1")
    n2 = Nodo("N2")
    n3 = Nodo("N3")

    vs = VSource("Vs", n1, gnd, idx_u=0)
    r  = Resistor("R", n1, n2, R)
    l  = Inductor("L", n2, n3, L)
    cap = Capacitor("C", n3, gnd, C)

    red = SubredMNA(
        "rlc",
        nodos=[gnd, n1, n2, n3],
        componentes=[vs, r, l, cap],
        dt=dt,
        mediciones_v=[(n3, gnd)],   # salida 0 = V_C
        mediciones_i=[vs],           # salida 1 = I_vs
    )
    return red


def test_mna_instancia():
    """El bloque se crea sin error y tiene n_state > 0."""
    red = _hacer_circuito_rlc()
    assert red.n_state > 0
    assert red.n_out == 2   # V_C e I_vs


def test_mna_rlc_serie():
    """
    Simula RLC serie en resonancia y verifica que la corriente crece
    (respuesta transitoria de energia creciente) antes de que la
    resistencia la amortigüe.  En resonancia ω = 1/sqrt(LC).
    """
    R, L, C = 10.0, 1e-3, 1e-6
    dt = 1e-7
    w0 = 1.0 / np.sqrt(L * C)          # ~31623 rad/s, frec ~5 kHz
    f0 = w0 / (2 * np.pi)

    m = Modelo(dt=dt)
    with m:
        src = m.add(FuenteSeno("vs", amplitud=100.0, frecuencia=f0))
        red = m.add(_hacer_circuito_rlc(R=R, L=L, C=C, dt=dt))
        m.conectar(src.salida, red.entrada[0:1])

    t_fin = 5.0 / f0   # 5 ciclos completos
    res = m.run(t_fin=t_fin, registrar=["rlc"])

    datos = res["rlc"]       # shape (n_steps, n_out)
    vc = datos[:, 0]         # Voltaje en capacitor
    i_vs = datos[:, 1]       # Corriente de la fuente

    # En resonancia, la corriente de estado estacionario debería ser Vm/R = 10 A
    # Esperamos que llegue a al menos la mitad de ese valor en 5 ciclos
    assert np.max(np.abs(vc)) > 1.0, "Voltaje en capacitor debe subir"
    assert np.max(np.abs(i_vs)) > 0.5, "Corriente debe ser > 0.5 A"


def test_mna_switch():
    """
    Circuito RC con interruptor: switch abierto -> V_C = 0.
    Switch cerrado -> V_C carga hacia V_source.
    """
    gnd = Nodo("GND")
    n1  = Nodo("N1")
    n2  = Nodo("N2")

    # Fuente V (entrada 0), interruptor controlado (entrada 1), RC carga
    vs  = VSource("Vs", n1, gnd, idx_u=0)
    sw  = Switch("S1", n1, n2, idx_ctrl=1, Ron=0.01, Roff=1e7)
    r   = Resistor("R", n2, gnd, 100.0)
    cap = Capacitor("C", n2, gnd, 1e-5)

    red = SubredMNA(
        "rc_sw",
        nodos=[gnd, n1, n2],
        componentes=[vs, sw, r, cap],
        dt=1e-4,
        mediciones_v=[(n2, gnd)],
        mediciones_i=[vs],
    )

    from bloques_crysi.bloques import FuenteConstante, FuenteEscalon

    m = Modelo(dt=1e-4)
    with m:
        v_src  = m.add(FuenteConstante("vdc", 100.0))    # 100 V fuente
        ctrl   = m.add(FuenteEscalon("ctrl", valor_final=1.0, t_paso=0.05))  # cierra a 50 ms
        bloque = m.add(red)
        m.conectar(v_src.salida, bloque.entrada[0:1])
        m.conectar(ctrl.salida,  bloque.entrada[1:2])

    res = m.run(t_fin=0.2, registrar=["rc_sw"])
    vc = res["rc_sw"][:, 0]

    n = len(vc)
    # La señal de control sube a t=0.05s -> paso en la mitad del total (t_fin=0.2s)
    cuarto = n // 4  # Primeros 50ms: switch abierto

    # Antes del cierre: con Roff=1e7 hay una pequeña corriente de fuga
    # V_C(t) = Vs * (1 - exp(-t/RC_off)) con RC_off = 1e7 * C = 100 s
    # A t=0.05s: V_C ~ 100 * 0.05/100 ~ 0.05 V  (¡muy pequeño!)
    assert np.max(np.abs(vc[:cuarto])) < 1.0, "Switch abierto -> V_C debe ser << 100V"
    # Después del cierre: voltaje debe haber subido significativamente (>>1V)
    assert vc[-1] > 40.0, "Switch cerrado -> V_C debe subir hacia 100V"
