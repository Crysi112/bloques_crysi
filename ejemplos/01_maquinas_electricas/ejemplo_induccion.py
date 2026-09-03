from bloques_crysi import (
    Modelo, FuenteTrifasica, FuenteConstante, MaquinaInduccion,
    LimitadorRapidez, Scope,
)

DT = 1e-4
WS = 2 * 3.1416 * 50 / (4 / 2)

def np_abs_max(a):
    return max(abs(float(x)) for x in a)

m = Modelo(dt=DT)
red = m.add(FuenteTrifasica("red", amplitud=310.0, frecuencia=50.0))
tl = m.add(FuenteConstante("tl", 0.0))
maq = m.add(MaquinaInduccion(
    "mi", rs=0.5, rr=0.4, Lm=0.1, Lls=0.005, Llr=0.005, P=4, J=0.5,
    Bm=0.01))
m.conectar(red.salida, maq.terminales)
m.conectar(tl.salida, maq.T_L)
sc = m.add(Scope("scope", anchos=[1, 1, 3],
                 guiones=["velocidad (rad/s)", "par (N.m)",
                          "ia (A)", "ib (A)", "ic (A)"]))
m.conectar(maq.sensorVelocidad(), sc.canales[0])
m.conectar(maq.sensorPar(), sc.canales[1])
m.conectar(maq.sensor3I(), sc.canales[2])
res = m.run(t_fin=2.0, registrar=[
    sc, maq.sensor3I(), maq.sensorVelocidad(), maq.sensorPar()])
wm, Te = res["wm"], res["Te"]
print("Maquina de induccion jaula de ardilla (arranque directo, sin carga)")
print("=" * 58)
print(f"velocidad sincrona : {WS:.2f} rad/s")
print(f"velocidad final   : {wm[-1]:9.3f} rad/s  "
      f"({100*wm[-1]/WS:.1f}% del sincronismo)")
print(f"par maximo en arranque : {np_abs_max(Te):9.3f} N.m")
print(f"corriente pico de fase  : {abs(res['I']).max():9.3f} A")

m = Modelo(dt=DT)
red = m.add(FuenteTrifasica("red", amplitud=310.0, frecuencia=50.0))
tl = m.add(FuenteConstante("tl", 10.0))
tl_ramp = m.add(LimitadorRapidez("tl_ramp", subida=10.0, bajada=10.0))
maq = m.add(MaquinaInduccion(
    "mi", rs=0.5, rr=0.4, Lm=0.1, Lls=0.005, Llr=0.005, P=4, J=0.5,
    Bm=0.01))
m.conectar(red.salida, maq.terminales)
m.conectar(tl.salida, tl_ramp.entrada)
m.conectar(tl_ramp.salida, maq.T_L)
sc = m.add(Scope("scope2", anchos=[1, 1, 3],
                 guiones=["velocidad (rad/s)", "par (N.m)",
                          "ia (A)", "ib (A)", "ic (A)"]))
m.conectar(maq.sensorVelocidad(), sc.canales[0])
m.conectar(maq.sensorPar(), sc.canales[1])
m.conectar(maq.sensor3I(), sc.canales[2])
res = m.run(t_fin=2.5, registrar=[
    sc, maq.sensor3I(), maq.sensorVelocidad(), maq.sensorPar()])
wm, Te = res["wm"], res["Te"]
print("Con carga 10 N.m (rampa de 10 N.m/s)")
print("-" * 58)
print(f"velocidad final   : {wm[-1]:9.3f} rad/s  "
      f"({100*wm[-1]/WS:.1f}% del sincronismo {WS:.2f} rad/s)")
print(f"par final en regimen    : {Te[-1]:9.3f} N.m  (carga {10.0} N.m)")
print(f"corriente pico de fase  : {abs(res['I']).max():9.3f} A")
