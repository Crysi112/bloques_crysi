"""
Demostración Completa de las Nuevas Características (Fases 1, 2 y 4)
------------------------------------------------------------------
Este script demuestra:
1. El uso del solver MNA (Modified Nodal Analysis) para circuitos.
2. El uso del modelo avanzado de Batería (BateriaECM).
3. La generación automática de código C (Codegen) de la topología.
"""
import numpy as np
import matplotlib.pyplot as plt

from bloques_crysi import Modelo, GeneradorCodigo
from bloques_crysi.bloques import FuenteConstante, FuenteEscalon, BateriaECM
from bloques_crysi.mna import SubredMNA, Nodo, Resistor, Capacitor, Inductor, VSource, Switch

def main():
    print("=== Configurando Modelo con MNA y Bateria ECM ===")
    
    gnd = Nodo("GND")
    n_in = Nodo("IN")
    n_mid = Nodo("MID")
    n_out = Nodo("OUT")
    
    v_in = VSource("Vin", n_in, gnd, idx_u=0)
    s_high = Switch("S_High", n_in, n_mid, idx_ctrl=1, Ron=0.01, Roff=1e6)
    s_low  = Switch("S_Low", gnd, n_mid, idx_ctrl=2, Ron=0.01, Roff=1e6)
    
    ind = Inductor("L", n_mid, n_out, L=1e-3)
    cap = Capacitor("C", n_out, gnd, C=1e-4)
    r_load = Resistor("Rload", n_out, gnd, R=10.0)
    
    mna_buck = SubredMNA(
        "Buck_MNA",
        nodos=[gnd, n_in, n_mid, n_out],
        componentes=[v_in, s_high, s_low, ind, cap, r_load],
        dt=1e-6,
        mediciones_v=[(n_out, gnd), (n_mid, gnd)],
        mediciones_i=[v_in]
    )

    m = Modelo(dt=1e-6)
    
    with m:
        bat = m.add(BateriaECM("BateriaMain", soc_init=0.85))
        
        pwm_high = m.add(FuenteEscalon("PWM_H", valor_final=1.0, t_paso=0.0))
        pwm_low  = m.add(FuenteConstante("PWM_L", valor=0.0))
        
        load_bat = m.add(FuenteConstante("LoadBat", valor=15.0))
        m.conectar(load_bat.salida, bat.entrada)
        
        buck = m.add(mna_buck)
        m.conectar(bat.salida[0], buck.entrada[0:1])
        m.conectar(pwm_high.salida, buck.entrada[1:2])
        m.conectar(pwm_low.salida, buck.entrada[2:3])
        
    print("\n--- Generando Codigo C (Codegen) ---")
    gen = GeneradorCodigo(m)
    codigo_c = gen.generar("mi_buck_bateria")
    with open("buck_generado.c", "w", encoding="utf-8") as f:
        f.write(codigo_c)
    print("-> Codigo C hiper-optimizado guardado en 'buck_generado.c'")
    
    print("\n--- Ejecutando Simulacion ---")
    res = m.run(t_fin=0.05, registrar=["BateriaMain", "Buck_MNA"])
    print("-> Simulacion completada.")
    
    t = res["BateriaMain"]["t"] if isinstance(res["BateriaMain"], dict) and "t" in res["BateriaMain"] else np.arange(len(res["BateriaMain"][:, 0])) * 1e-6
    
    plt.figure(figsize=(10, 8))
    
    plt.subplot(3, 1, 1)
    plt.plot(t, res["BateriaMain"][:, 0], label="Voltaje Bateria (V_term)", color="blue")
    plt.ylabel("Voltaje [V]")
    plt.title("Bateria ECM alimentando circuito MNA")
    plt.legend()
    plt.grid(True)
    
    plt.subplot(3, 1, 2)
    plt.plot(t, res["Buck_MNA"][:, 0], label="Voltaje Salida Buck (MNA)", color="orange")
    plt.ylabel("Voltaje [V]")
    plt.legend()
    plt.grid(True)
    
    plt.subplot(3, 1, 3)
    plt.plot(t, res["BateriaMain"][:, 1] * 100, label="SOC Bateria (%)", color="green")
    plt.xlabel("Tiempo [s]")
    plt.ylabel("SOC [%]")
    plt.legend()
    plt.grid(True)
    
    plt.tight_layout()
    plt.savefig("demo_completa_plot.png")
    print("-> Grafica guardada como 'demo_completa_plot.png'")

if __name__ == "__main__":
    main()
