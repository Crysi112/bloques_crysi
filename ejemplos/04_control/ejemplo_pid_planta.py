import numpy as np

from bloques_crysi import (
    Modelo, FuenteEscalon, Suma, PID, FuncionTransferencia, Scope,
    LimitadorRapidez, FiltroPasoBajo, Puerto,
)

DT = 1e-3

def par(bloque, tipo, offset, n):
    return Puerto(bloque, tipo, offset, n)

def np_abs_max(a):
    return max(abs(float(x)) for x in a)

m = Modelo(dt=DT)
ref = m.add(FuenteEscalon("ref", valor_final=1.0, t_paso=1e-5))
planta = m.add(FuncionTransferencia("planta", num=[1.0],
                                    den=[1.0, 2.0, 1.0]))
m.conectar(ref.salida, planta.entrada)
res = m.run(t_fin=10.0, registrar=[ref, planta])
print("Lazo abierto: G(s) = 1/(s^2 + 2s + 1) ante escalon unitario")
print("=" * 58)
print(f"salida final : {res['planta'][-1]:.4f}  (esperado 1.0)")
print(f"polos de G(s): -1 y -1 (doble) -> respuesta sin oscilar")
print(f"valor a 1 s  : {res['planta'][int(1.0/DT)]:.4f}  "
      f"(1 - e^-1·(1 + 1) = {1 - np.exp(-1) * 2:.4f})")

m = Modelo(dt=DT)
ref = m.add(FuenteEscalon("ref", valor_final=1.0, t_paso=1e-5))
err = m.add(Suma("err", (1.0, -1.0)))
pid = m.add(PID("pid", Kp=2.0, Ki=0.8, Kd=0.3, Tf=0.01,
                u_min=-10.0, u_max=10.0))
actuador = m.add(LimitadorRapidez("actuador", subida=4.0, bajada=4.0))
planta = m.add(FuncionTransferencia("planta", num=[1.0],
                                    den=[1.0, 2.0, 1.0]))
sensor = m.add(FiltroPasoBajo("sensor", fc=10.0, orden=1))
m.conectar(ref.salida, par(err, "ent", 0, 1))
m.conectar(sensor.salida, par(err, "ent", 1, 1))
m.conectar(err.salida, pid.entrada)
m.conectar(pid.salida, actuador.entrada)
m.conectar(actuador.salida, planta.entrada)
m.conectar(planta.salida, sensor.entrada)

sc = m.add(Scope("scope", max_canales=3,
                 guiones=["referencia", "salida", "salida filtrada"]))
m.conectar(ref.salida, sc.canales[0])
m.conectar(planta.salida, sc.canales[1])
m.conectar(sensor.salida, sc.canales[2])

res = m.run(t_fin=100.0, registrar=[sc, ref, err, pid, actuador,
                                    planta, sensor])
salida = res["planta"]
print("Lazo cerrado con PID (Kp=2, Ki=0.8, Kd=0.3) y actuador limitado")
print("=" * 58)
print("t(s)   ref       salida   accion PID")
paso = 5000
for i in range(0, len(salida), paso):
    print(f"{res.t[i]:5.2f}  {res['ref'][i]:8.3f}  {salida[i]:8.3f}"
          f"  {res['pid'][i]:9.3f}")
print(f"\nsalida final: {salida[-1]:.4f}  (objetivo: 1.0)")
if abs(salida[-1] - 1.0) < 0.05:
    print("El PID sigue la referencia en regimen.")
else:
    print("OJO: la salida no convergio a 1.0.")
print(f"pendiente maxima de la accion: "
      f"{np_abs_max(np.diff(np.asarray(res['actuador'])))/DT:.3f} 1/s "
      f"(limitada a 4.0 1/s por el actuador)")
