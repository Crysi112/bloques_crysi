import pytest
import numpy as np
from bloques_crysi import Modelo
from bloques_crysi.bloques import FuenteSeno
from bloques_crysi.mna import SubredMNA, Nodo, Resistor, VSource, Diodo

def test_mna_diodo():
    """
    Circuito rectificador de media onda:
    V_ac -> Diodo -> Resistor -> GND
    """
    gnd = Nodo("GND")
    n1  = Nodo("N1")
    n2  = Nodo("N2")

    vs  = VSource("Vs", n1, gnd, idx_u=0)
    d1  = Diodo("D1", n1, n2, Vf=0.7, Ron=1e-3, Roff=1e5)
    r   = Resistor("R", n2, gnd, 10.0)

    red = SubredMNA(
        "rectificador",
        nodos=[gnd, n1, n2],
        componentes=[vs, d1, r],
        dt=1e-4,
        mediciones_v=[(n2, gnd)], # Salida 0: V_R (voltaje rectificado)
        mediciones_i=[vs]
    )

    m = Modelo(dt=1e-4)
    with m:
        vac = m.add(FuenteSeno("vac", amplitud=10.0, frecuencia=50.0))
        bloque = m.add(red)
        m.conectar(vac.salida, bloque.entrada[0:1])

    res = m.run(t_fin=0.04, registrar=["vac", "rectificador"])
    
    v_in_raw = np.asarray(res["vac"])
    v_in = v_in_raw[:,0] if v_in_raw.ndim==2 else v_in_raw
    v_out = np.asarray(res["rectificador"])[:, 0]
    
    # El diodo debe rectificar: V_out debe ser positivo o cero
    assert np.all(v_out >= -0.1), "El diodo no debe conducir en inversa"
    
    # Cuando V_in > 0.7, V_out deberia ser aprox V_in - 0.7
    idx_pico = np.argmax(v_in)
    assert v_out[idx_pico] > 8.0, "El diodo debe conducir en directa"
