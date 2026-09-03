"""
Demostración de Co-Simulación Norton (Fase 4)
---------------------------------------------
Muestra cómo un modelo EMT (alta frecuencia) interactúa con un
flujo de potencia estacionario utilizando el equivalente Norton
sin necesidad de iteraciones ni rebobinados de estado.
"""
import numpy as np
import matplotlib.pyplot as plt

from bloques_crysi import Modelo
from bloques_crysi.bloques import FuenteSeno, MedidorPotencia
from bloques_crysi.red.backend_base import BackendRedMock, TensionBus
from bloques_crysi.red.norton import CoSimNorton, PuertoNorton

def main():
    print("=== Configurando Modelo EMT ===")
    m = Modelo(dt=1e-4)
    
    with m:
        # La red se representa como una fuente de voltaje en el modelo EMT
        red = m.add(FuenteSeno("Red_AC", amplitud=311.12, frecuencia=50.0))
        
        # Medidor para Fase 1
        medidor = m.add(MedidorPotencia("PCC", fases=1))
        m.conectar(red.salida, medidor.entrada)
        
        # Corriente de carga
        corriente = m.add(FuenteSeno("Carga_I", amplitud=10.0, frecuencia=50.0, fase=-np.pi/4))
        m.conectar(corriente.salida, medidor.corrientes)

    print("\n=== Configurando Red de Flujo de Potencia (Mock) ===")
    backend = BackendRedMock(v_slack_v=220.0)
    
    # Inicializar el modelo para obtener IDs de las señales
    m.iniciar()
    idx_v = red.salida.indices()[0]
    idx_i = corriente.salida.indices()[0]
    
    puerto_1 = PuertoNorton(
        bus_pcc=1, 
        idx_v_pcc=idx_v, 
        idx_ia=idx_i, 
        idx_ib=idx_i, 
        idx_ic=idx_i,
        v_nominal_ln=220.0
    )
    
    cosim = CoSimNorton(
        modelo=m,
        backend=backend,
        puertos=[puerto_1],
        dt_red=0.02, # Comunicación cada 20 ms (1 ciclo)
        alpha=0.8
    )
    
    print("\n--- Ejecutando Co-Simulacion Fuerte (Norton) ---")
    res = cosim.run(t_fin=0.2)
    print("-> Co-Simulacion completada.")
    
    t = res[1]["t"]
    v_pcc = res[1]["V_pcc"]
    p_pcc = res[1]["P"]
    
    plt.figure(figsize=(10, 6))
    
    plt.subplot(2, 1, 1)
    plt.plot(t, v_pcc, marker='o', label="V_pcc (Magnitud)")
    plt.ylabel("Voltaje [V]")
    plt.title("Co-simulacion Norton: Tension y Potencia Activa en el PCC")
    plt.grid(True)
    plt.legend()
    
    plt.subplot(2, 1, 2)
    plt.plot(t, p_pcc / 1000.0, marker='x', color='red', label="P_pcc (Inyectada)")
    plt.xlabel("Tiempo [s]")
    plt.ylabel("Potencia Activa [kW]")
    plt.grid(True)
    plt.legend()
    
    plt.tight_layout()
    plt.savefig("demo_cosim_plot.png")
    print("-> Grafica guardada como 'demo_cosim_plot.png'")

if __name__ == "__main__":
    main()
