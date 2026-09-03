import numpy as np
from bloques_crysi import Modelo
from bloques_crysi.bloques import FuenteConstante, Integrador, FuncionTransferencia

def test_multirate_tf():
    # Creamos un modelo a 1ms
    m = Modelo(dt=1e-3)
    
    # Fuente
    f = m.add(FuenteConstante("f", 1.0))
    
    # Agregamos una Función de Transferencia 1/(s+1) y le decimos que corra cada 10ms (Ts=0.01)
    tf = m.add(FuncionTransferencia("tf", num=[1], den=[1, 1], Ts=0.01))
    
    m.conectar(f.salida, tf.entrada)
    
    res = m.run(0.1, registrar=[tf])
    
    # Como la TF solo se evalúa cada 10ms, su salida debería verse como escalones (ZOH) en pasos de 10 puntos (ya que dt=1ms)
    y = res["tf"]
    
    # Verificamos que se mantenga constante entre pasos
    # Entre t=0.001 y t=0.009 (índices 1 a 9), el valor debe ser constante
    assert np.allclose(y[1:10], y[1])
    
    # Entre t=0.010 y t=0.019 (índices 10 a 19), debe ser constante y mayor que antes
    assert np.allclose(y[10:20], y[10])
    assert y[10] > y[1]
