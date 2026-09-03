import numpy as np
from . import opcodes as ops
from .bloques import Bloque
from .puertos import Puerto
GREF = 1000.0
TREF = 25.0
def _resolver_panel(V, p, G=GREF, T=TREF):
    Ns, Np, Voc, Isc, ki = p[0], p[1], p[2], p[3], p[6]
    Rs, Rsh, n = p[7], p[8], p[9]
    k_voc = p[10] if len(p) > 10 else 0.0
    Vt = 0.02585 * (T + 273.15) / 298.15
    a = Ns * n * Vt
    Iph = Np * (Isc + ki * (T - TREF)) * (G / GREF)
    Voc_T = Voc * (1.0 + k_voc * (T - TREF))
    I0 = (Iph - Voc_T / Rsh) / (np.exp(Voc_T / a) - 1.0)
    I = Iph - V / Rsh
    for _ in range(10):
        u = (V + I * Rs) / a
        if u > 700.0:
            u = 700.0
        e = np.exp(u)
        g = Iph - I0 * (e - 1.0) - (V + I * Rs) / Rsh - I
        gp = -I0 * (Rs / a) * e - Rs / Rsh - 1.0
        dI = g / gp
        lam, g1 = 1.0, abs(g)
        for _ in range(4):
            In = I - lam * dI
            u2 = (V + In * Rs) / a
            if u2 > 700.0:
                u2 = 700.0
            g2 = abs(Iph - I0 * (np.exp(u2) - 1.0)
                     - (V + In * Rs) / Rsh - In)
            if g2 <= g1 or lam <= 0.125:
                break
            lam *= 0.5
        I -= lam * dI
    return I
class PanelSolar(Bloque):
    op = ops.OP_PANEL_SOLAR
    n_in = 3
    n_out = 1
    n_state = 0
    etiqueta = "PanelSolar"
    NOMBRES = ["I"]
    def __init__(self, nombre, Ns=60, Np=1, Voc=37.6, Isc=9.12, Vmp=29.9,
                 Imp=8.63, ki=0.004, Rs=0.2, Rsh=400.0, n=1.3, k_voc=-0.003):
        super().__init__(nombre)
        if Ns <= 0 or Np <= 0 or Voc <= 0 or Isc <= 0 or Rsh <= 0 or n <= 0:
            raise ValueError("Ns, Np, Voc, Isc, Rsh y n deben ser > 0.")
        if ki < 0:
            raise ValueError("ki no puede ser negativo.")
        self.param = [float(Ns), float(Np), float(Voc), float(Isc),
                      float(Vmp), float(Imp), float(ki), float(Rs),
                      float(Rsh), float(n), float(k_voc)]
        self._datos = dict(Vmp=float(Vmp), Imp=float(Imp))
        self.etiqueta = (f"PanelSolar (Voc={Voc:g} V, Isc={Isc:g} A, "
                         f"{Ns:g}x{Np:g} celdas)")
        self.entrada = Puerto(self, "ent", 0, 3)
        self.salida = Puerto(self, "sal", 0, 1)
    def curvaIV(self, G=GREF, T=TREF, n=400):
        Vs = np.linspace(0.0, 1.15 * self.param[2], n)
        Is = np.array([_resolver_panel(float(v), self.param, G, T)
                       for v in Vs])
        return Vs, Is, Vs * Is
    def graficar_curvas(self, Gs=(200.0, 400.0, 600.0, 800.0, 1000.0),
                        T=25.0, n=300):
        import matplotlib.pyplot as plt
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
        for G in Gs:
            Vs, Is, Ps = self.curvaIV(G=G, T=T, n=n)
            ax1.plot(Vs, Is, lw=2, label=f"G = {G:g} W/m²")
            ax2.plot(Vs, Ps, lw=2)
        ax1.set_xlabel("V [V]")
        ax1.set_ylabel("I [A]")
        ax1.set_title(f"Panel {self.etiqueta} a {T:g} °C")
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        ax2.set_xlabel("V [V]")
        ax2.set_ylabel("P [W]")
        ax2.set_title("Potencia")
        ax2.grid(True, alpha=0.3)
        fig.tight_layout()
        plt.show()
