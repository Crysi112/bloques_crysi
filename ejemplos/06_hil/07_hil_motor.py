import numpy as np
import matplotlib.pyplot as plt
from bloques_crysi import (
    Modelo, MotorHardwareCH32, FuenteSeno, LimitadorRapidez, Scope,
)

COM_PORT = 3

m = Modelo(dt=0.002)

ref_duty = m.add(FuenteSeno("Ref_Duty", amplitud=1, frecuencia=0.2))

limit = m.add(LimitadorRapidez("Limit_Duty", subida=2.0, bajada=2.0))
m.conectar(ref_duty.salida, limit.entrada)

motor = m.add(MotorHardwareCH32("MotorFisico", puerto_com=COM_PORT, baudrate=115200))

m.conectar(limit.salida, motor.duty)

print(f"Iniciando HIL en COM{COM_PORT}. ¡Asegúrate de que el motor esté conectado!")
print("La simulación avanzará en tiempo real...")
try:
    res = m.run(t_fin=10.0, registrar=[
        motor.sensorAngulo(),
        motor.sensorVelocidad()
    ])
    print("Simulación terminada.")

    t = res.t
    ang = res["MotorFisico_ang"]
    rpm = res["MotorFisico_rpm"]

    plt.figure(figsize=(10, 6))

    plt.subplot(2, 1, 1)
    plt.plot(t, ang, label='Ángulo (raw)')
    plt.ylabel("Posición")
    plt.legend()
    plt.grid(True)

    plt.subplot(2, 1, 2)
    plt.plot(t, rpm, label='Velocidad (RPM)', color='orange')
    plt.xlabel("Tiempo [s]")
    plt.ylabel("RPM")
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show()

except Exception as e:
    print("Error durante HIL:", e)
