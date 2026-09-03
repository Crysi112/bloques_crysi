from bloques_crysi import (
    Modelo, FuenteTrifasica, FuenteConstante, MaquinaImanesPermanentes,
    LimitadorRapidez, Scope,
)

DT = 1e-6

def np_abs_max(a):
    return max(abs(float(x)) for x in a)

m = Modelo(dt=DT)
red = m.add(FuenteTrifasica("red", amplitud=310.0, frecuencia=50.0))
tl = m.add(FuenteConstante("tl", 0.0))
maq = m.add(MaquinaImanesPermanentes(
    "pmsm", rs=0.1, Ld=1e-3, Lq=1e-3, lam_m=0.1, P=6, J=0.01, Bm=0.001))
m.conectar(red.salida, maq.terminales)
m.conectar(tl.salida, maq.T_L)
res = m.run(t_fin=0.5, registrar=[
    maq.sensor3V(), maq.sensor3I(), maq.sensorVelocidad(),
    maq.sensorPosicion(), maq.sensorPosicionElectrica(), maq.sensorPar(),
])
wm, Te = res["wm"], res["Te"]
print("PMAC (imanes permanentes) conectada a red trifasica sin control")
print("=" * 72)
print(f"velocidad final  : {wm[-1]:9.3f} rad/s")
print(f"|velocidad| max  : {np_abs_max(wm):9.3f} rad/s")
print(f"par final        : {Te[-1]:9.3f} N.m")
print(f"|par| max        : {np_abs_max(Te):9.3f} N.m")
print(f"posicion mecc    : {res['th_rm'][-1]:9.3f} rad")
print(f"posicion elec    : {res['th_e'][-1]:9.3f} rad")
print(f"|V| pico-fase    : {res['V'][:, 0].max():9.3f} V")
print(f"|I| pico         : {abs(res['I']).max():9.3f} A")
print()
print("Sin devanados amortiguadores ni control FOC, la PMAC no se")
print("sincroniza sola con la red: la velocidad oscila.")

m = Modelo(dt=DT)
red = m.add(FuenteTrifasica("red", amplitud=310.0, frecuencia=50.0))
tl = m.add(FuenteConstante("tl", 5.0))
tl_ramp = m.add(LimitadorRapidez("tl_ramp", subida=10.0, bajada=10.0))
maq = m.add(MaquinaImanesPermanentes(
    "pmsm", rs=0.1, Ld=1e-3, Lq=1e-3, lam_m=0.1, P=6, J=0.01, Bm=0.001))
m.conectar(red.salida, maq.terminales)
m.conectar(tl.salida, tl_ramp.entrada)
m.conectar(tl_ramp.salida, maq.T_L)
sc = m.add(Scope("scope", anchos=[1, 1, 3],
                 guiones=["velocidad (rad/s)", "par (N.m)",
                          "ia (A)", "ib (A)", "ic (A)"]))
m.conectar(maq.sensorVelocidad(), sc.canales[0])
m.conectar(maq.sensorPar(), sc.canales[1])
m.conectar(maq.sensor3I(), sc.canales[2])
res = m.run(t_fin=1, registrar=[
    sc,
    maq.sensor3V(), maq.sensor3I(), maq.sensorVelocidad(),
    maq.sensorPosicion(), maq.sensorPosicionElectrica(), maq.sensorPar(),
])
wm, Te = res["wm"], res["Te"]
print("Con carga rampa de 10 N.m/s hasta 5.0 N.m")
print("-" * 72)
print(f"velocidad final  : {wm[-1]:9.3f} rad/s")
print(f"|velocidad| max  : {np_abs_max(wm):9.3f} rad/s")
print(f"par final        : {Te[-1]:9.3f} N.m  (carga: 5.0 N.m)")
print(f"|par| max        : {np_abs_max(Te):9.3f} N.m")
print(f"|I| pico         : {abs(res['I']).max():9.3f} A")
print()
print("Sigue sin sincronizarse: hace falta un lazo cerrado (FOC).")
