import numpy as np
from . import opcodes as ops
from .bloques import Bloque
from .puertos import Puerto
TIPOS = {
    "plomo_acido": dict(
        alias=("plomo-acido", "plomo-ácido", "plomo_ácido", "lead_acid",
               "lead-acid", "la"),
        Vnom=12.0, Qrated=100.0, f_qmax=1.09,
        f_vfull=0.15, f_a=0.05, f_qexp=0.10, f_qnom=0.90, rho=0.02),
    "litio": dict(
        alias=("lithium", "li", "litio_ion", "li_ion", "li-ion"),
        Vnom=200.0, Qrated=100.0, f_qmax=1.04,
        f_vfull=0.07, f_a=0.03, f_qexp=0.05, f_qnom=0.90, rho=0.012),
    "niquel_cadmio": dict(
        alias=("niquel-cadmio", "níquel-cadmio", "nicd", "ni-cd",
               "nickel_cadmium", "nickel-cadmium"),
        Vnom=12.0, Qrated=100.0, f_qmax=1.04,
        f_vfull=0.20, f_a=0.08, f_qexp=0.10, f_qnom=0.90, rho=0.02),
    "niquel_metal_hidruro": dict(
        alias=("niquel-metal-hidruro", "níquel-metal-hidruro", "nimh",
               "ni-mh", "nickel_metal_hydride", "nickel-metal-hydride"),
        Vnom=200.0, Qrated=100.0, f_qmax=1.04,
        f_vfull=0.18, f_a=0.07, f_qexp=0.10, f_qnom=0.90, rho=0.015),
}
def _tipo(tipo):
    t = str(tipo).strip().lower().replace(" ", "_")
    if t in TIPOS:
        return t, TIPOS[t]
    for nombre, config in TIPOS.items():
        if t in config["alias"]:
            return nombre, config
    raise ValueError(
        f"Tipo de bateria desconocido: {tipo!r}. Validos: "
        + ", ".join(TIPOS))
class Bateria(Bloque):
    op = ops.OP_BATERIA
    n_in = 1
    n_out = 3
    n_state = 4
    etiqueta = "Bateria"
    NOMBRES = ["vbat", "SOC", "T"]
    def __init__(self, nombre, tipo="plomo_acido", Vnom=None, Qrated=None,
                 SOCinit=100.0, Qmax=None, tau=30.0, R=None, eta_c=1.0,
                 histeresis=False, termica=True, C_th=30_000.0,
                 R_th=0.3, alpha=0.0035, T_amb=25.0):
        super().__init__(nombre)
        tipo, cfg = _tipo(tipo)
        self.tipo = tipo
        Vnom = float(Vnom if Vnom is not None else cfg["Vnom"])
        Qrated = float(Qrated if Qrated is not None else cfg["Qrated"])
        if Vnom <= 0 or Qrated <= 0:
            raise ValueError("Vnom y Qrated deben ser > 0.")
        if not (10.0 < SOCinit <= 100.0):
            raise ValueError("SOCinit debe estar en (10, 100] %.")
        tau = float(tau)
        if tau <= 0:
            raise ValueError("tau debe ser > 0.")
        eta_c = float(eta_c)
        if not (0.0 < eta_c <= 1.0):
            raise ValueError("eta_c debe estar en (0, 1].")
        Q = float(Qmax if Qmax is not None else cfg["f_qmax"] * Qrated)
        if Q < Qrated:
            raise ValueError("Qmax no puede ser menor que Qrated.")
        Qexp = cfg["f_qexp"] * Q
        Qnom = cfg["f_qnom"] * Q
        Vfull = Vnom * (1.0 + cfg["f_vfull"])
        A = cfg["f_a"] * Vnom
        Vexp = Vfull - A
        B = 3.0 / Qexp
        i_c = Qrated / 3600.0
        K = ((Vfull - Vnom - A * (1.0 - np.exp(-B * Qnom)))
             * (Q - Qnom) / (Q * (Qnom + i_c)))
        R = float(R if R is not None else cfg["rho"] * Vnom / Qrated)
        E0 = Vfull + (R + K) * i_c - A
        Vcap = 1.25 * Vfull
        C_th = float(C_th if termica else 0.0)
        if C_th < 0.0:
            raise ValueError("C_th debe ser >= 0.")
        if R_th <= 0.0:
            raise ValueError("R_th debe ser > 0.")
        T_amb = float(T_amb)
        self.param = [float(E0), float(K), float(Q), float(A), float(B),
                      float(R), float(tau), float(Vcap), float(eta_c),
                      float(histeresis), float(C_th), float(R_th),
                      float(alpha), float(T_amb)]
        it0 = (1.0 - SOCinit / 100.0) * Q
        self.estados_iniciales = [it0, 0.0, A * np.exp(-B * it0), T_amb]
        self.SOCinit = float(SOCinit)
        self._curva = dict(Vfull=Vfull, Vexp=Vexp, Qexp=Qexp,
                           Qnom=Qnom, Vnom=Vnom, Q=Q, i_c=i_c)
        self.etiqueta = (f"Bateria {tipo} ({Vnom} V, {Qrated} Ah, "
                         f"SOC0={SOCinit:g}%)")
        self.entrada = Puerto(self, "ent", 0, 1)
        self.salida = Puerto(self, "sal", 0, 3)
    def curva(self, corriente=None, n=200, carga=False):
        E0, K, Q, A, B, R = self.param[:6]
        i = float(corriente) if corriente is not None else self._curva["i_c"]
        if not carga and i <= 0:
            i = self._curva["i_c"]
        if carga and i >= 0:
            i = -self._curva["i_c"]
        qmin = 0.1 * Q
        qmax = 0.9 * Q
        vcap = self.param[7]
        its = np.linspace(0.0, qmax, n)
        v = np.zeros(n)
        for k, it in enumerate(its):
            if carga:
                E = (E0 - K * Q * i / (it - qmin) - K * Q * it / (Q - it)
                     + A * np.exp(-B * it))
                E = min(E, vcap)
            else:
                E = (E0 - K * Q * (it + i) / (Q - it) + A * np.exp(-B * it))
            v[k] = E - R * i
        return its, v
    def graficar_curvas(self, corriente=None, n=200):
        import matplotlib.pyplot as plt
        its_d, v_d = self.curva(corriente=corriente, n=n, carga=False)
        its_c, v_c = self.curva(corriente=corriente, n=n, carga=True)
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.plot(its_d, v_d, label="descarga", lw=2)
        ax.plot(its_c, v_c, label="carga", lw=2)
        c = self._curva
        ax.axvline(c["Qexp"], ls="--", color="gray", lw=1)
        ax.axvline(c["Qnom"], ls="--", color="gray", lw=1)
        ax.plot([0.0, c["Qnom"]], [c["Vfull"], c["Vnom"]], "o", ms=4,
                label="puntos de ajuste (Vfull, Qexp, Qnom)")
        ax.set_title(f"Bateria {self.tipo}: curvas de carga/descarga "
                     f"({self.param[2]:g} Ah)")
        ax.set_xlabel("carga extraida it [Ah]")
        ax.set_ylabel("tension [V]")
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.tight_layout()
        plt.show()
    def sensorVoltaje(self):
        from .puertos import Sensor
        return Sensor("Vbat", self, "sal", 0, 1, canales=["Vbat"])
    def sensorSOC(self):
        from .puertos import Sensor
        return Sensor("SOC", self, "sal", 1, 1, canales=["SOC"])
    def sensorTemperatura(self):
        from .puertos import Sensor
        return Sensor("T", self, "sal", 2, 1, canales=["T"])
