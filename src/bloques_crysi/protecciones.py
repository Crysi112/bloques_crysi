"""Módulo de Protecciones Eléctricas según Estándar ANSI/IEEE C37.2."""

from __future__ import annotations
import math
import cmath
from typing import Optional, Sequence, Tuple, Union, Dict, Any
import numpy as np

from .bloques import Bloque
from .puertos import Puerto, Sensor
from . import opcodes as ops

CURVAS_51 = {
    "IEEE_MODERATELY_INVERSE": {"A": 0.0515, "B": 0.1140, "p": 0.02},
    "IEEE_VERY_INVERSE":       {"A": 19.61,  "B": 0.4910, "p": 2.00},
    "IEEE_EXTREMELY_INVERSE":  {"A": 28.20,  "B": 0.1217, "p": 2.00},
    "IEC_STANDARD_INVERSE":    {"A": 0.1400, "B": 0.0000, "p": 0.02},
    "IEC_VERY_INVERSE":        {"A": 13.500, "B": 0.0000, "p": 1.00},
    "IEC_EXTREMELY_INVERSE":   {"A": 80.000, "B": 0.0000, "p": 2.00},
    "IEC_LONG_TIME_INVERSE":   {"A": 120.00, "B": 0.0000, "p": 1.00},
}


class ReleProteccion(Bloque):
    def __init__(self, nombre: str, codigo_ansi: str, n_in: int, n_out: int = 3, n_state: int = 4):
        super().__init__(nombre)
        self.codigo_ansi = codigo_ansi
        self.n_in = n_in
        self.n_out = n_out
        self.n_state = n_state
        self.estados_iniciales = [0.0] * n_state
        self.etiqueta = f"ANSI {codigo_ansi} ({nombre})"
        self.entrada = Puerto(self, "ent", 0, self.n_in)
        self.salida = Puerto(self, "sal", 0, self.n_out, canales=["Trip", "PickUp", "Info"])


class Rele50(ReleProteccion):
    op = ops.OP_RELE_50
    def __init__(self, nombre: str, I_pickup: float, t_retardo: float = 0.0, sufijo: str = "P"):
        super().__init__(nombre, f"50{sufijo}", n_in=1, n_out=3, n_state=2)
        if I_pickup <= 0:
            raise ValueError("I_pickup debe ser > 0")
        self.I_pickup = float(I_pickup)
        self.t_retardo = max(0.0, float(t_retardo))
        self.param = [self.I_pickup, self.t_retardo]
        self.estados_iniciales = [0.0, 0.0]

    def paso(self, i_medida: float, dt: float) -> Tuple[float, float, float]:
        i_mag = abs(i_medida)
        pickup = 1.0 if i_mag >= self.I_pickup else 0.0
        t_acum = self.estados_iniciales[0]
        trip = self.estados_iniciales[1]

        if pickup:
            t_acum += dt
            if t_acum >= self.t_retardo:
                trip = 1.0
        else:
            t_acum = 0.0

        self.estados_iniciales[0] = t_acum
        self.estados_iniciales[1] = trip
        return trip, pickup, t_acum


class Rele51(ReleProteccion):
    op = ops.OP_RELE_51
    def __init__(self, nombre: str, I_pickup: float, dial_tds: float = 1.0,
                 tipo_curva: str = "IEC_STANDARD_INVERSE", sufijo: str = "P",
                 emulacion_disco: bool = True):
        super().__init__(nombre, f"51{sufijo}", n_in=1, n_out=3, n_state=3)
        if tipo_curva not in CURVAS_51:
            raise ValueError(f"Curva desconocida: {tipo_curva}. Válidas: {list(CURVAS_51.keys())}")
        c = CURVAS_51[tipo_curva]
        self.I_pickup = float(I_pickup)
        self.tds = float(dial_tds)
        self.curva = c
        self.emulacion_disco = emulacion_disco
        self.param = [self.I_pickup, self.tds, c["A"], c["B"], c["p"], 1.0 if emulacion_disco else 0.0]
        self.estados_iniciales = [0.0, 0.0, 0.0]

    def calcular_tiempo_trip(self, i_mag: float) -> float:
        m = i_mag / self.I_pickup
        if m <= 1.001:
            return float("inf")
        c = self.curva
        return self.tds * (c["A"] / ((m ** c["p"]) - 1.0) + c["B"])

    def paso(self, i_medida: float, dt: float) -> Tuple[float, float, float]:
        i_mag = abs(i_medida)
        pickup = 1.0 if i_mag >= self.I_pickup else 0.0
        disco = self.estados_iniciales[0]
        trip = self.estados_iniciales[1]

        if pickup:
            t_trip = self.calcular_tiempo_trip(i_mag)
            disco += dt / max(1e-4, t_trip)
            if disco >= 1.0:
                disco = 1.0
                trip = 1.0
        else:
            if self.emulacion_disco:
                disco = max(0.0, disco - (dt / (self.tds * 5.0)))
            else:
                disco = 0.0

        self.estados_iniciales[0] = disco
        self.estados_iniciales[1] = trip
        return trip, pickup, disco


class Rele67(ReleProteccion):
    op = ops.OP_RELE_67
    def __init__(self, nombre: str, I_pickup: float, rca_deg: float = 45.0,
                 arco_disparo_deg: float = 85.0, t_retardo: float = 0.05):
        super().__init__(nombre, "67", n_in=2, n_out=3, n_state=2)
        self.I_pickup = float(I_pickup)
        self.rca_rad = math.radians(rca_deg)
        self.arco_rad = math.radians(arco_disparo_deg)
        self.t_retardo = float(t_retardo)
        self.param = [self.I_pickup, self.rca_rad, self.arco_rad, self.t_retardo]
        self.estados_iniciales = [0.0, 0.0]

    def evaluar_zona(self, v_complex: complex, i_complex: complex) -> bool:
        if abs(i_complex) < self.I_pickup or abs(v_complex) < 1e-3:
            return False
        ang_v = cmath.phase(v_complex)
        ang_i = cmath.phase(i_complex)
        delta = (ang_i - ang_v) - self.rca_rad
        delta = (delta + math.pi) % (2 * math.pi) - math.pi
        return abs(delta) <= self.arco_rad


class Rele46(ReleProteccion):
    op = ops.OP_RELE_46
    def __init__(self, nombre: str, I2_pickup: float, t_retardo: float = 3.0):
        super().__init__(nombre, "46", n_in=1, n_out=3, n_state=2)
        if I2_pickup <= 0:
            raise ValueError("I2_pickup debe ser > 0")
        self.I2_pickup = float(I2_pickup)
        self.t_retardo = max(0.0, float(t_retardo))
        self.param = [self.I2_pickup, self.t_retardo]
        self.estados_iniciales = [0.0, 0.0]

    def paso(self, i2_medida: float, dt: float) -> Tuple[float, float, float]:
        i_mag = abs(i2_medida)
        pickup = 1.0 if i_mag >= self.I2_pickup else 0.0
        t_acum = self.estados_iniciales[0]
        trip = self.estados_iniciales[1]

        if pickup:
            t_acum += dt
            if t_acum >= self.t_retardo:
                trip = 1.0
        else:
            t_acum = 0.0

        self.estados_iniciales[0] = t_acum
        self.estados_iniciales[1] = trip
        return trip, pickup, t_acum


class Rele62(ReleProteccion):
    op = ops.OP_RELE_62
    def __init__(self, nombre: str, t_retardo: float = 300.0, umbral: float = 0.5):
        super().__init__(nombre, "62", n_in=1, n_out=3, n_state=2)
        self.t_retardo = max(0.0, float(t_retardo))
        self.umbral = float(umbral)
        self.param = [self.t_retardo, self.umbral]
        self.estados_iniciales = [0.0, 0.0]

    def paso(self, senal: float, dt: float) -> Tuple[float, float, float]:
        pickup = 1.0 if senal >= self.umbral else 0.0
        t_acum = self.estados_iniciales[0]
        trip = self.estados_iniciales[1]

        if pickup:
            t_acum += dt
            if t_acum >= self.t_retardo:
                trip = 1.0
        else:
            t_acum = 0.0

        self.estados_iniciales[0] = t_acum
        self.estados_iniciales[1] = trip
        return trip, pickup, t_acum


class Rele63(ReleProteccion):
    op = ops.OP_RELE_63
    def __init__(self, nombre: str, umbral: float = 0.5, t_retardo: float = 0.0):
        super().__init__(nombre, "63", n_in=1, n_out=3, n_state=2)
        self.umbral = float(umbral)
        self.t_retardo = max(0.0, float(t_retardo))
        self.param = [self.umbral, self.t_retardo]
        self.estados_iniciales = [0.0, 0.0]

    def paso(self, senal: float, dt: float) -> Tuple[float, float, float]:
        pickup = 1.0 if senal >= self.umbral else 0.0
        t_acum = self.estados_iniciales[0]
        trip = self.estados_iniciales[1]

        if pickup:
            t_acum += dt
            if t_acum >= self.t_retardo:
                trip = 1.0
        else:
            t_acum = 0.0

        self.estados_iniciales[0] = t_acum
        self.estados_iniciales[1] = trip
        return trip, pickup, t_acum


class Rele49(ReleProteccion):
    op = ops.OP_RELE_49
    def __init__(self, nombre: str, I_nominal: float, tau_segundos: float = 600.0,
                 alarma_pct: float = 90.0, disparo_pct: float = 100.0):
        super().__init__(nombre, "49", n_in=1, n_out=3, n_state=2)
        self.I_nom = float(I_nominal)
        self.tau = float(tau_segundos)
        self.theta_alarm = float(alarma_pct) / 100.0
        self.theta_trip = float(disparo_pct) / 100.0
        self.param = [self.I_nom, self.tau, self.theta_alarm, self.theta_trip]
        self.estados_iniciales = [0.0, 0.0]


class Rele27(ReleProteccion):
    op = ops.OP_RELE_27
    def __init__(self, nombre: str, V_pickup: float, t_retardo: float = 0.5):
        super().__init__(nombre, "27", n_in=1, n_out=3, n_state=2)
        self.V_pickup = float(V_pickup)
        self.t_retardo = float(t_retardo)
        self.param = [self.V_pickup, self.t_retardo]
        self.estados_iniciales = [0.0, 0.0]

    def paso(self, v_medida: float, dt: float) -> Tuple[float, float, float]:
        v_mag = abs(v_medida)
        pickup = 1.0 if v_mag <= self.V_pickup else 0.0
        t_acum = self.estados_iniciales[0]
        trip = self.estados_iniciales[1]

        if pickup:
            t_acum += dt
            if t_acum >= self.t_retardo:
                trip = 1.0
        else:
            t_acum = 0.0
        self.estados_iniciales[0] = t_acum
        self.estados_iniciales[1] = trip
        return trip, pickup, t_acum


class Rele59(ReleProteccion):
    op = ops.OP_RELE_59
    def __init__(self, nombre: str, V_pickup: float, t_retardo: float = 0.2, sufijo: str = "P"):
        super().__init__(nombre, f"59{sufijo}", n_in=1, n_out=3, n_state=2)
        self.V_pickup = float(V_pickup)
        self.t_retardo = float(t_retardo)
        self.param = [self.V_pickup, self.t_retardo]
        self.estados_iniciales = [0.0, 0.0]

    def paso(self, v_medida: float, dt: float) -> Tuple[float, float, float]:
        v_mag = abs(v_medida)
        pickup = 1.0 if v_mag >= self.V_pickup else 0.0
        t_acum = self.estados_iniciales[0]
        trip = self.estados_iniciales[1]

        if pickup:
            t_acum += dt
            if t_acum >= self.t_retardo:
                trip = 1.0
        else:
            t_acum = 0.0
        self.estados_iniciales[0] = t_acum
        self.estados_iniciales[1] = trip
        return trip, pickup, t_acum


class Rele81(ReleProteccion):
    op = ops.OP_RELE_81
    def __init__(self, nombre: str, f_pickup: float, modo: str = "U", t_retardo: float = 0.1):
        codigo = "81U" if modo.upper() == "U" else "81O"
        super().__init__(nombre, codigo, n_in=1, n_out=3, n_state=2)
        self.f_pickup = float(f_pickup)
        self.modo = modo.upper()
        self.t_retardo = float(t_retardo)
        self.param = [self.f_pickup, 0.0 if self.modo == "U" else 1.0, self.t_retardo]
        self.estados_iniciales = [0.0, 0.0]

    def paso(self, f_medida: float, dt: float) -> Tuple[float, float, float]:
        pickup = 1.0 if (f_medida <= self.f_pickup if self.modo == "U" else f_medida >= self.f_pickup) else 0.0
        t_acum = self.estados_iniciales[0]
        trip = self.estados_iniciales[1]

        if pickup:
            t_acum += dt
            if t_acum >= self.t_retardo:
                trip = 1.0
        else:
            t_acum = 0.0
        self.estados_iniciales[0] = t_acum
        self.estados_iniciales[1] = trip
        return trip, pickup, t_acum


class Rele24(ReleProteccion):
    op = ops.OP_RELE_24
    def __init__(self, nombre: str, v_f_ratio_pickup: float = 1.10, t_retardo: float = 2.0):
        super().__init__(nombre, "24", n_in=2, n_out=3, n_state=2)
        self.ratio_pickup = float(v_f_ratio_pickup)
        self.t_retardo = float(t_retardo)
        self.param = [self.ratio_pickup, self.t_retardo]
        self.estados_iniciales = [0.0, 0.0]


class Rele87(ReleProteccion):
    op = ops.OP_RELE_87
    def __init__(self, nombre: str, I_min: float = 0.2, slope1: float = 0.25,
                 slope2: float = 0.60, I_break: float = 2.0, sufijo: str = "T",
                 bloqueo_armonicos: bool = True):
        super().__init__(nombre, f"87{sufijo}", n_in=2, n_out=4, n_state=2)
        self.I_min = float(I_min)
        self.s1 = float(slope1)
        self.s2 = float(slope2)
        self.I_bp = float(I_break)
        self.bloqueo_armonicos = bloqueo_armonicos
        self.param = [self.I_min, self.s1, self.s2, self.I_bp, 1.0 if bloqueo_armonicos else 0.0]
        self.salida = Puerto(self, "sal", 0, 4, canales=["Trip", "I_diff", "I_rest", "HarmonicBlock"])

    def evaluar_disparo(self, i1: complex, i2: complex, pct_2do_armonico: float = 0.0) -> Tuple[bool, float, float]:
        i_diff = abs(i1 + i2)
        i_rest = (abs(i1) + abs(i2)) / 2.0

        if self.bloqueo_armonicos and pct_2do_armonico > 0.15:
            return False, i_diff, i_rest

        if i_rest <= self.I_bp:
            umbral = self.I_min + self.s1 * i_rest
        else:
            umbral = self.I_min + self.s1 * self.I_bp + self.s2 * (i_rest - self.I_bp)

        return (i_diff > umbral), i_diff, i_rest


class Rele21(ReleProteccion):
    op = ops.OP_RELE_21
    def __init__(self, nombre: str, z_alcance_ohm: complex, angulo_mho_deg: float = 75.0,
                 t_retardo: float = 0.0, zona: int = 1):
        super().__init__(nombre, f"21_Z{zona}", n_in=2, n_out=4, n_state=2)
        self.z_reach = complex(z_alcance_ohm)
        self.t_retardo = float(t_retardo)
        self.centro = self.z_reach / 2.0
        self.radio = abs(self.z_reach) / 2.0
        self.param = [self.centro.real, self.centro.imag, self.radio, self.t_retardo]
        self.salida = Puerto(self, "sal", 0, 4, canales=["Trip", "R_aparente", "X_aparente", "EnZona"])

    def paso(self, v_c: complex, i_c: complex, dt: float) -> Tuple[float, float, float, float]:
        if abs(i_c) < 1e-4:
            return 0.0, 999.0, 999.0, 0.0
        z_aparente = v_c / i_c
        en_zona = abs(z_aparente - self.centro) <= self.radio
        t_acum = self.estados_iniciales[0]
        trip = self.estados_iniciales[1]

        if en_zona:
            t_acum += dt
            if t_acum >= self.t_retardo:
                trip = 1.0
        else:
            t_acum = 0.0

        self.estados_iniciales[0] = t_acum
        self.estados_iniciales[1] = trip
        return trip, z_aparente.real, z_aparente.imag, 1.0 if en_zona else 0.0


class Disyuntor52(Bloque):
    op = ops.OP_DISYUNTOR_52
    n_in = 2
    n_out = 3
    n_state = 2
    etiqueta = "Disyuntor 52"

    def __init__(self, nombre: str, t_apertura_mecanica: float = 0.05, cerrado_inicial: bool = True):
        super().__init__(nombre)
        self.codigo_ansi = "52"
        self.t_apertura = float(t_apertura_mecanica)
        self.param = [self.t_apertura]
        self.estados_iniciales = [1.0 if cerrado_inicial else 0.0, -1.0]
        self.entrada_trip = Puerto(self, "ent", 0, 1)
        self.entrada_i = Puerto(self, "ent", 1, 1)
        self.salida = Puerto(self, "sal", 0, 3, canales=["Estado", "52a", "52b"])

    def paso(self, orden_trip: float, i_circulante: float, dt: float) -> Tuple[float, float, float]:
        estado = self.estados_iniciales[0]
        timer = self.estados_iniciales[1]

        if estado == 1.0:
            if orden_trip > 0.5 and timer < 0:
                timer = 0.0

            if timer >= 0.0:
                timer += dt
                if timer >= self.t_apertura:
                    if abs(i_circulante) < 1e-1:
                        estado = 0.0
                        timer = -1.0
        else:
            if orden_trip < 0.5:
                pass

        self.estados_iniciales[0] = estado
        self.estados_iniciales[1] = timer
        return estado, estado, 1.0 - estado


class Rele86(Bloque):
    op = ops.OP_RELE_86
    n_in = 2
    n_out = 2
    n_state = 1
    etiqueta = "Lockout 86"

    def __init__(self, nombre: str):
        super().__init__(nombre)
        self.codigo_ansi = "86"
        self.estados_iniciales = [0.0]
        self.entrada_trip = Puerto(self, "ent", 0, 1)
        self.entrada_reset = Puerto(self, "ent", 1, 1)
        self.salida = Puerto(self, "sal", 0, 2, canales=["Bloqueo", "Listo"])

    def paso(self, trip_in: float, reset_in: float) -> Tuple[float, float]:
        bloqueado = self.estados_iniciales[0]
        if trip_in > 0.5:
            bloqueado = 1.0
        elif reset_in > 0.5:
            bloqueado = 0.0

        self.estados_iniciales[0] = bloqueado
        return bloqueado, 1.0 - bloqueado


class Rele79(Bloque):
    op = ops.OP_RELE_79
    n_in = 2
    n_out = 2
    n_state = 4
    etiqueta = "Reconectador 79"

    def __init__(self, nombre: str, tiempos_muertos: Sequence[float] = (0.5, 15.0, 30.0),
                 t_reinicio_secuencia: float = 60.0):
        super().__init__(nombre)
        self.codigo_ansi = "79"
        self.tiempos_muertos = list(tiempos_muertos)
        self.max_intentos = len(self.tiempos_muertos)
        self.t_reinicio = float(t_reinicio_secuencia)
        self.estados_iniciales = [0.0, 0.0, 0.0, 0.0]
        self.salida = Puerto(self, "sal", 0, 2, canales=["PulsoRecierre", "Lockout"])


class Rele25(Bloque):
    op = ops.OP_RELE_25
    n_in = 6
    n_out = 2
    etiqueta = "Sincrocheck 25"

    def __init__(self, nombre: str, delta_v_max: float = 0.10,
                 delta_f_max: float = 0.10, delta_ang_deg: float = 15.0):
        super().__init__(nombre)
        self.codigo_ansi = "25"
        self.dv_max = float(delta_v_max)
        self.df_max = float(delta_f_max)
        self.dang_rad = math.radians(delta_ang_deg)
        self.param = [self.dv_max, self.df_max, self.dang_rad]
        self.salida = Puerto(self, "sal", 0, 2, canales=["PermisoCierre", "DeltaAngDeg"])


def crear_proteccion(codigo_ansi: str, nombre: str, **kwargs) -> Bloque:
    cod = codigo_ansi.upper().strip()

    if cod.startswith("50"):
        sufijo = cod[2:] or "P"
        return Rele50(nombre, I_pickup=kwargs.get("I_pickup", 100.0),
                      t_retardo=kwargs.get("t_retardo", 0.0), sufijo=sufijo)

    elif cod.startswith("51"):
        sufijo = cod[2:] or "P"
        return Rele51(nombre, I_pickup=kwargs.get("I_pickup", 100.0),
                      dial_tds=kwargs.get("dial_tds", kwargs.get("TMS", 1.0)),
                      tipo_curva=kwargs.get("tipo_curva", "IEC_STANDARD_INVERSE"),
                      sufijo=sufijo)

    elif cod.startswith("27"):
        return Rele27(nombre, V_pickup=kwargs.get("V_pickup", 0.85),
                      t_retardo=kwargs.get("t_retardo", 0.5))

    elif cod.startswith("59"):
        sufijo = cod[2:] or "P"
        return Rele59(nombre, V_pickup=kwargs.get("V_pickup", 1.15),
                      t_retardo=kwargs.get("t_retardo", 0.2), sufijo=sufijo)

    elif cod.startswith("81"):
        modo = "O" if "O" in cod else "U"
        return Rele81(nombre, f_pickup=kwargs.get("f_pickup", 49.0 if modo == "U" else 51.0),
                      modo=modo, t_retardo=kwargs.get("t_retardo", 0.1))

    elif cod.startswith("87"):
        sufijo = cod[2:] or "T"
        return Rele87(nombre, I_min=kwargs.get("I_min", 0.2),
                      slope1=kwargs.get("slope1", 0.25),
                      slope2=kwargs.get("slope2", 0.60),
                      I_break=kwargs.get("I_break", 2.0),
                      sufijo=sufijo)

    elif cod.startswith("21"):
        zona = 1
        if "_Z" in cod:
            zona = int(cod.split("_Z")[-1])
        return Rele21(nombre, z_alcance_ohm=kwargs.get("z_alcance_ohm", 5.0 + 10.0j),
                      t_retardo=kwargs.get("t_retardo", 0.0), zona=zona)

    elif cod == "52":
        return Disyuntor52(nombre, t_apertura_mecanica=kwargs.get("t_apertura_mecanica", 0.05),
                           cerrado_inicial=kwargs.get("cerrado_inicial", True))

    elif cod == "86":
        return Rele86(nombre)

    elif cod == "79":
        return Rele79(nombre, tiempos_muertos=kwargs.get("tiempos_muertos", (0.5, 15.0, 30.0)))

    elif cod == "25":
        return Rele25(nombre)

    elif cod.startswith("49"):
        return Rele49(nombre, I_nominal=kwargs.get("I_nominal", 100.0))

    elif cod.startswith("46"):
        return Rele46(nombre, I2_pickup=kwargs.get("I2_pickup", kwargs.get("I_pickup", 10.0)),
                      t_retardo=kwargs.get("t_retardo", 3.0))

    elif cod.startswith("62"):
        return Rele62(nombre, t_retardo=kwargs.get("t_retardo", 300.0),
                      umbral=kwargs.get("umbral", 0.5))

    elif cod.startswith("63"):
        return Rele63(nombre, umbral=kwargs.get("umbral", 0.5),
                      t_retardo=kwargs.get("t_retardo", 0.0))

    else:
        raise NotImplementedError(
            f"El código ANSI {codigo_ansi} no tiene plantilla automática aún. "
            f"Usa un bloque base o compón uno con comparadores y temporizadores."
        )