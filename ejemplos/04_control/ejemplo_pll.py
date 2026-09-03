import matplotlib.pyplot as plt

from bloques_crysi import Modelo, FuenteTrifasica, PLLTrifasico, FiltroPasoBajo
from bloques_crysi import Puerto

DT = 1e-4

def Puerto_pll(pll, k):
    return Puerto(pll, "sal", k, 1)

def _corre():
    m = Modelo(dt=DT)
    red = m.add(FuenteTrifasica("red", amplitud=311.0, frecuencia=50.0))
    pll = m.add(PLLTrifasico("pll", Kp=10.0, Ki=100.0, f_ff=50.0))
    m.conectar(red.salida, pll.entrada)
    w_f = m.add(FiltroPasoBajo("w_f", fc=100.0, orden=1))
    m.conectar(Puerto_pll(pll, 0), w_f.entrada)
    res = m.run(t_fin=0.2, registrar=[pll.salida, w_f])
    return res

res = _corre()
t = res.t
w = res["pll"][:, 0]
th = res["pll"][:, 1]
w_filt = res["w_f"]
print("PLL trifásico sobre red de 50 Hz")
print("=" * 46)
print(f"w final     : {w[-1]:.4f} rad/s  "
      f"(esperado 2*pi*50 = {2*3.14159265358979*50:.4f})")
print(f"th final    : {th[-1]:.4f} rad  "
      f"(esperado w*t - pi/2 = {w[-1]*t[-1]-3.14159265358979/2:.4f})")
print(f"error de f  : {abs(w[-1] - 2*3.14159265358979*50):.2e} rad/s")
print(f"w filtrada  : {w_filt[-1]:.4f} rad/s  "
      f"(FiltroPasoBajo fc=100 Hz)")

res = _corre()
t = res.t
w = res["pll"][:, 0]
th = res["pll"][:, 1]
w_filt = res["w_f"]
plt.figure(figsize=(10, 4))
plt.plot(t, w, label="w (rad/s)")
plt.plot(t, th, label="th (rad)")
plt.plot(t, w_filt, label="w filtrada (rad/s)", lw=2)
plt.xlabel("tiempo [s]")
plt.grid(True, alpha=0.3)
plt.legend()
plt.title("PLL trifásico: enganche de frecuencia y fase")
plt.show()
