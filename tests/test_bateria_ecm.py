import pytest
from bloques_crysi import Modelo, FuenteConstante
from bloques_crysi.bloques import BateriaECM

def test_bateria_ecm():
    m = Modelo(dt=0.1)
    with m:
        load = m.add(FuenteConstante("load", 10.0))
        bat = m.add(BateriaECM("bat", soc_init=0.9))
        m.conectar(load.salida, bat.entrada)
    
    res = m.run(t_fin=10.0, registrar=["bat"])
    
    assert res["bat"][:, 1][-1] < 0.9 # SOC
    v = res["bat"][:, 0][-1]          # V_term
    assert v > 300.0 and v < 450.0

    print("Test de BateriaECM superado!")
