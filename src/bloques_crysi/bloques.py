import math
from typing import List, Optional, Sequence, Tuple, Union
import numpy as np
from . import opcodes as ops
from .puertos import Puerto
class Bloque:
    op: Optional[int] = None
    n_in: int = 0
    n_out: int = 1
    n_state: int = 0
    etiqueta: str = "Base"
    def __init__(self, nombre: str, Ts: float = 0.0) -> None:
        self.nombre: str = nombre
        self.Ts: float = float(Ts)
        self.param: List[float] = []
        self.estados_iniciales: List[float] = [0.0] * self.n_state
        self.entrada: Optional[Puerto] = Puerto(self, "ent", 0, self.n_in) if self.n_in else None
        self.salida: Puerto = Puerto(self, "sal", 0, self.n_out)
    def __repr__(self) -> str:
        return f"<{self.etiqueta} {self.nombre!r}>"
class FuenteConstante(Bloque):
    op = ops.OP_SRC_CONST
    etiqueta = "FuenteConstante"
    def __init__(self, nombre: str, valor: float = 0.0) -> None:
        super().__init__(nombre)
        self.param = [float(valor)]
        self.etiqueta = f"FuenteConst({valor})"
class FuenteEscalon(Bloque):
    op = ops.OP_SRC_STEP
    etiqueta = "FuenteEscalon"
    def __init__(self, nombre: str, valor_final: float = 1.0,
                 t_paso: float = 0.0, valor_inicial: float = 0.0) -> None:
        super().__init__(nombre)
        self.param = [float(valor_final), float(t_paso), float(valor_inicial)]
class FuenteRampa(Bloque):
    op = ops.OP_SRC_RAMP
    etiqueta = "FuenteRampa"
    def __init__(self, nombre: str, pendiente: float = 1.0,
                 t_inicio: float = 0.0, offset: float = 0.0) -> None:
        super().__init__(nombre)
        self.param = [float(pendiente), float(t_inicio), float(offset)]
class FuenteSeno(Bloque):
    op = ops.OP_SRC_SIN
    etiqueta = "FuenteSeno"
    def __init__(self, nombre: str, amplitud: float = 1.0,
                 frecuencia: float = 1.0, fase: float = 0.0,
                 offset: float = 0.0) -> None:
        super().__init__(nombre)
        self.param = [float(amplitud), float(frecuencia), float(fase), float(offset)]
class FuenteTrifasica(Bloque):
    op = ops.OP_SRC_TRIF
    n_out = 3
    etiqueta = "FuenteTrifasica"
    def __init__(self, nombre: str, amplitud: float = 1.0,
                 frecuencia: float = 50.0, fase: float = 0.0) -> None:
        super().__init__(nombre)
        self.param = [float(amplitud), float(frecuencia), float(fase)]
        self.etiqueta = f"Fuente3F({amplitud} V, {frecuencia} Hz)"
class Ganancia(Bloque):
    op = ops.OP_GAIN
    n_in = 1
    etiqueta = "Ganancia"
    def __init__(self, nombre: str, valor: float = 1.0) -> None:
        super().__init__(nombre)
        self.param = [float(valor)]
class Suma(Bloque):
    op = ops.OP_SUM
    etiqueta = "Suma"
    def __init__(self, nombre: str, signos: Sequence[float] = (1.0, -1.0)) -> None:
        super().__init__(nombre)
        self.n_in = len(signos)
        self.param = [float(s) for s in signos]
        self.etiqueta = f"Suma({signos})"
        self.entrada = Puerto(self, "ent", 0, self.n_in)
class Integrador(Bloque):
    op = ops.OP_INTEGRADOR
    n_in = 1
    n_state = 1
    etiqueta = "Integrador"
    def __init__(self, nombre: str, valor_inicial: float = 0.0, Ts: float = 0.0) -> None:
        super().__init__(nombre, Ts=Ts)
        self.estados_iniciales = [float(valor_inicial)]
class FuncionTransferencia(Bloque):
    op = ops.OP_TF
    n_in = 1
    etiqueta = "FuncionTransferencia"
    def __init__(self, nombre: str, num: Sequence[float],
                 den: Sequence[float], Ts: float = 0.0) -> None:
        super().__init__(nombre, Ts=Ts)
        self.orden = max(len(num), len(den)) - 1
        self.param = self._tustin(num, den)
        self.n_state = 2 * self.orden
        self.estados_iniciales = [0.0] * self.n_state

    def _tustin(self, num, den, dt=0.0):
        self._num = [float(c) for c in num]
        self._den = [float(c) for c in den]
        return [0.0]
    def _discretizar(self, dt):
        from scipy.signal import bilinear
        import numpy as np
        num, den = self._num, self._den
        bd, ad = bilinear(num, den, fs=1.0 / dt)
        n = len(ad) - 1
        bd_padded = list(bd) + [0.0] * (n + 1 - len(bd))
        bd_norm = (np.asarray(bd_padded) / ad[0]).tolist()
        ad_norm = (np.asarray(ad) / ad[0]).tolist()
        return [float(n)] + bd_norm + ad_norm[1:]
class PID(Bloque):
    op = ops.OP_PID
    n_in = 1
    n_state = 3
    etiqueta = "PID"
    def __init__(self, nombre: str, Kp: float = 1.0, Ki: float = 0.0,
                 Kd: float = 0.0, Tf: float = 0.0,
                 u_min: float = -1e12, u_max: float = 1e12,
                 Ts: float = 0.0) -> None:
        super().__init__(nombre, Ts=Ts)
        self.param = [float(Kp), float(Ki), float(Kd), float(Tf),
                      float(u_min), float(u_max)]
class Clarke(Bloque):
    op = ops.OP_CLARKE
    n_in = 3
    n_out = 2
    etiqueta = "Clarke"
    def __init__(self, nombre):
        super().__init__(nombre)
class InvClarke(Bloque):
    op = ops.OP_INV_CLARKE
    n_in = 2
    n_out = 3
    etiqueta = "InvClarke"
    def __init__(self, nombre):
        super().__init__(nombre)
class Park(Bloque):
    op = ops.OP_PARK
    n_in = 3
    n_out = 2
    etiqueta = "Park"
    def __init__(self, nombre):
        super().__init__(nombre)
class InvPark(Bloque):
    op = ops.OP_INV_PARK
    n_in = 3
    n_out = 2
    etiqueta = "InvPark"
    def __init__(self, nombre):
        super().__init__(nombre)
class TransformadaQD(Bloque):
    op = ops.OP_QD
    n_in = 7
    n_out = 4
    etiqueta = "TransformadaQD"
    def __init__(self, nombre):
        super().__init__(nombre)
        self.vabc = Puerto(self, "ent", 0, 3, canales=["va", "vb", "vc"])
        self.iabc = Puerto(self, "ent", 3, 3, canales=["ia", "ib", "ic"])
        self.th = Puerto(self, "ent", 6, 1, canales=["th"])
        self.salida = Puerto(self, "sal", 0, 4,
                             canales=["vqs", "vds", "iqs", "ids"])
class Saturar(Bloque):
    op = ops.OP_SATURAR
    n_in = 1
    etiqueta = "Saturar"
    def __init__(self, nombre, u_min=0.0, u_max=1.0):
        super().__init__(nombre)
        if u_min > u_max:
            raise ValueError("u_min no puede ser mayor que u_max.")
        self.param = [float(u_min), float(u_max)]
class Relay(Bloque):
    op = ops.OP_RELAY
    n_in = 1
    n_out = 1
    n_state = 1
    etiqueta = "Relay"
    def __init__(self, nombre, umbral_on=1.0, umbral_off=0.0,
                 salida_on=1.0, salida_off=0.0):
        super().__init__(nombre)
        if umbral_off > umbral_on:
            raise ValueError("umbral_off no puede ser mayor que umbral_on.")
        self.param = [float(umbral_on), float(umbral_off),
                      float(salida_on), float(salida_off)]
        self.estados_iniciales = [float(salida_off)]
class PulsoRectangular(Bloque):
    op = ops.OP_PULSO_RECT
    n_out = 1
    etiqueta = "PulsoRectangular"
    def __init__(self, nombre, amplitud=1.0, periodo=1.0, duty=0.5,
                 fase=0.0, offset=0.0):
        super().__init__(nombre)
        if periodo <= 0:
            raise ValueError("periodo debe ser > 0.")
        if not (0.0 < duty <= 1.0):
            raise ValueError("duty debe estar en (0, 1].")
        self.param = [float(amplitud), float(periodo), float(duty),
                      float(fase), float(offset)]
class Display(Bloque):
    op = ops.OP_GAIN
    n_in = 1
    n_out = 1
    etiqueta = "Display"
    def __init__(self, nombre, formato="%.4g"):
        super().__init__(nombre)
        self.formato = formato
        self.param = [1.0]
class PLLTrifasico(Bloque):
    op = ops.OP_PLL
    n_in = 3
    n_out = 2
    n_state = 2
    etiqueta = "PLLTrifasico"
    def __init__(self, nombre, Kp=10.0, Ki=100.0, f_ff=50.0, theta0=0.0):
        super().__init__(nombre)
        if f_ff <= 0:
            raise ValueError("f_ff debe ser > 0.")
        self.param = [float(Kp), float(Ki), 2.0 * math.pi * f_ff, float(theta0)]
        self.estados_iniciales = [float(theta0), 0.0]
class EjeMecanico(Bloque):
    op = ops.OP_EJE_MECANICO
    n_out = 2
    n_state = 2
    etiqueta = "EjeMecanico"
    def __init__(self, nombre, n_maquinas=2, J_eq=0.1, Bm_eq=0.0):
        super().__init__(nombre)
        self.param = [float(J_eq), float(Bm_eq)]
        self.estados_iniciales = [0.0, 0.0]
        self.n_in = n_maquinas + 1
        self.entradas = Puerto(self, "ent", 0, self.n_in)
        self.salida = Puerto(self, "sal", 0, 2)
class MedidorPotencia(Bloque):
    op = ops.OP_MEDIDOR_POTENCIA
    etiqueta = "MedidorPotencia"
    def __init__(self, nombre: str, fases: int = 3, con_Q: bool = True):
        super().__init__(nombre)
        if fases not in (1, 3):
            raise ValueError("fases debe ser 1 (DC) o 3 (AC).")
        self.fases = fases
        self.n_out = 3 if (fases == 3 and con_Q) else 2
        self.n_in = 2 * fases + 2
        self.NOMBRES = (["P_e", "Q_e", "P_m"] if self.n_out == 3
                        else ["P_e", "P_m"])
        self.entrada = Puerto(self, "ent", 0, fases)
        self.corrientes = Puerto(self, "ent", fases, fases)
        self.mecanica = Puerto(self, "ent", 2 * fases, 2)
        self.salida = Puerto(self, "sal", 0, self.n_out,
                             canales=list(self.NOMBRES))
        self.opcionales = 2
class InterruptorIdeal(Bloque):
    op = ops.OP_INTERRUPTOR
    n_in = 3
    n_out = 2
    etiqueta = "InterruptorIdeal"
    def __init__(self, nombre: str, R_on: float = 1e-3,
                 R_off: float = 1e6):
        super().__init__(nombre)
        self.param = [float(R_on), float(R_off), 0.0]
        self.NOMBRES = ["I_sw", "V_sw"]
        self.control = Puerto(self, "ent", 0, 1)
        self.terminales = Puerto(self, "ent", 1, 2)
        self.salida = Puerto(self, "sal", 0, 2)
class DiodoIdeal(Bloque):
    op = ops.OP_DIODO
    n_in = 2
    n_out = 2
    n_state = 1
    etiqueta = "DiodoIdeal"
    def __init__(self, nombre: str, R_on: float = 1e-3,
                 R_off: float = 1e6, V_f: float = 0.0,
                 histeresis: float = 1e-3):
        super().__init__(nombre)
        self.param = [float(R_on), float(R_off), float(V_f),
                      float(histeresis)]
        self.estados_iniciales = [0.0]
        self.NOMBRES = ["I_sw", "V_sw"]
        self.terminales = Puerto(self, "ent", 0, 2)
        self.salida = Puerto(self, "sal", 0, 2)
class PuenteInversorTrifasico(Bloque):
    op = ops.OP_PUENTE_INV_3F
    n_in = 4
    n_out = 3
    etiqueta = "PuenteInversorTrifasico"
    def __init__(self, nombre: str, promediado: bool = False):
        super().__init__(nombre)
        self.param = [1.0 if promediado else 0.0]
        self.NOMBRES = ["va", "vb", "vc"]
        self.entrada = Puerto(self, "ent", 0, 4)
        self.salida = Puerto(self, "sal", 0, 3)
class PuenteInversorMonofasico(Bloque):
    op = ops.OP_PUENTE_INV_1F
    n_in = 3
    n_out = 1
    etiqueta = "PuenteInversorMonofasico"
    def __init__(self, nombre: str, promediado: bool = False):
        super().__init__(nombre)
        self.param = [1.0 if promediado else 0.0]
        self.NOMBRES = ["Vout"]
        self.entrada = Puerto(self, "ent", 0, 3)
        self.salida = Puerto(self, "sal", 0, 1)
class Multiplexor(Bloque):
    op = ops.OP_MUX
    etiqueta = "Multiplexor"
    def __init__(self, nombre, n_canales=2):
        super().__init__(nombre)
        if n_canales < 1:
            raise ValueError("n_canales debe ser >= 1.")
        self.n_in = n_canales
        self.n_out = n_canales
        self.entradas = [Puerto(self, "ent", k, 1) for k in range(n_canales)]
        self.salida = Puerto(self, "sal", 0, n_canales)
class Demultiplexor(Bloque):
    op = ops.OP_DEMUX
    etiqueta = "Demultiplexor"
    def __init__(self, nombre, n_canales=2, entrada=None):
        super().__init__(nombre)
        if entrada is not None:
            n_canales = getattr(entrada, "n", None) or n_canales
        if n_canales < 1:
            raise ValueError("n_canales debe ser >= 1.")
        self.n_in = n_canales
        self.n_out = n_canales
        self.entrada = Puerto(self, "ent", 0, n_canales)
        self.salidas = [Puerto(self, "sal", k, 1) for k in range(n_canales)]
        if entrada is not None:
            e = entrada.puerto if hasattr(entrada, "puerto") else entrada
            self._pendientes_conectar = [(e, self.entrada)]
    def __iter__(self):
        return iter(self.salidas)
    def __getitem__(self, item):
        return self.salidas[item]
class Tabla1D(Bloque):
    op = ops.OP_LUT1D
    n_in = 1
    etiqueta = "Tabla1D"
    def __init__(self, nombre, puntos_x, valores_y):
        super().__init__(nombre)
        puntos_x = [float(v) for v in puntos_x]
        valores_y = [float(v) for v in valores_y]
        if len(puntos_x) != len(valores_y):
            raise ValueError("puntos_x y valores_y deben tener igual longitud.")
        if len(puntos_x) < 2:
            raise ValueError("Se necesitan al menos 2 puntos.")
        if any(puntos_x[i + 1] <= puntos_x[i] for i in range(len(puntos_x) - 1)):
            raise ValueError("puntos_x debe ser estrictamente creciente.")
        self.param = [float(len(puntos_x))] + puntos_x + valores_y
class Tabla2D(Bloque):
    op = ops.OP_LUT2D
    n_in = 2
    etiqueta = "Tabla2D"
    def __init__(self, nombre, puntos_x, puntos_y, tabla):
        super().__init__(nombre)
        puntos_x = [float(v) for v in puntos_x]
        puntos_y = [float(v) for v in puntos_y]
        if len(puntos_x) < 2 or len(puntos_y) < 2:
            raise ValueError("Se necesitan al menos 2 puntos por eje.")
        if any(puntos_x[i + 1] <= puntos_x[i] for i in range(len(puntos_x) - 1)):
            raise ValueError("puntos_x debe ser estrictamente creciente.")
        if any(puntos_y[i + 1] <= puntos_y[i] for i in range(len(puntos_y) - 1)):
            raise ValueError("puntos_y debe ser estrictamente creciente.")
        filas = [list(f) for f in tabla]
        if len(filas) != len(puntos_y) or \
           any(len(f) != len(puntos_x) for f in filas):
            raise ValueError(
                "tabla debe tener len(puntos_y) filas de len(puntos_x) "
                "columnas (fila i = puntos_y[i], columna j = puntos_x[j])."
            )
        z = [float(v) for f in filas for v in f]
        self.param = ([float(len(puntos_x)), float(len(puntos_y))]
                      + puntos_x + puntos_y + z)
        self.entrada1 = Puerto(self, "ent", 0, 1)
        self.entrada2 = Puerto(self, "ent", 1, 1)
class Tabla3D(Bloque):
    op = ops.OP_LUT3D
    n_in = 3
    etiqueta = "Tabla3D"
    def __init__(self, nombre, puntos_x, puntos_y, puntos_z, tabla):
        super().__init__(nombre)
        puntos_x = [float(v) for v in puntos_x]
        puntos_y = [float(v) for v in puntos_y]
        puntos_z = [float(v) for v in puntos_z]
        if len(puntos_x) < 2 or len(puntos_y) < 2 or len(puntos_z) < 2:
            raise ValueError("Se necesitan al menos 2 puntos por eje.")
        for nombre_eje, pts in (("puntos_x", puntos_x),
                                ("puntos_y", puntos_y),
                                ("puntos_z", puntos_z)):
            if any(pts[i + 1] <= pts[i] for i in range(len(pts) - 1)):
                raise ValueError(f"{nombre_eje} debe ser estrictamente creciente.")
        capas = [list(c) for c in tabla]
        if (not capas or len(capas) != len(puntos_z)
                or any(len(c) != len(puntos_y) for c in capas)
                or any(len(f) != len(puntos_x) for f in capas[0])):
            raise ValueError(
                "tabla debe tener len(puntos_z) capas de len(puntos_y) "
                "filas de len(puntos_x) columnas."
            )
        z = [float(v) for c in capas for f in c for v in f]
        self.param = ([float(len(puntos_x)), float(len(puntos_y)),
                       float(len(puntos_z))]
                      + puntos_x + puntos_y + puntos_z + z)
        self.entrada1 = Puerto(self, "ent", 0, 1)
        self.entrada2 = Puerto(self, "ent", 1, 1)
        self.entrada3 = Puerto(self, "ent", 2, 1)
class Logico(Bloque):
    op = ops.OP_LOGICO
    etiqueta = "Logico"
    OPS = {"AND": 0, "OR": 1, "NAND": 2, "NOR": 3,
           "XOR": 4, "XNOR": 5, "NOT": 6}
    def __init__(self, nombre, opcion="AND", n_entradas=2, umbral=0.5):
        super().__init__(nombre)
        if opcion not in self.OPS:
            raise ValueError(f"opcion debe ser una de {list(self.OPS)}.")
        if opcion == "NOT":
            n_entradas = 1
        if n_entradas < 1:
            raise ValueError("n_entradas debe ser >= 1.")
        self.n_in = n_entradas
        self.param = [float(self.OPS[opcion]), float(umbral)]
        self.entradas = [Puerto(self, "ent", k, 1) for k in range(n_entradas)]
class Relacional(Bloque):
    op = ops.OP_RELACIONAL
    n_in = 2
    etiqueta = "Relacional"
    OPS = {"==": 0, "!=": 1, "<": 2, "<=": 3, ">": 4, ">=": 5}
    def __init__(self, nombre, opcion="<", tol=0.0):
        super().__init__(nombre)
        if opcion not in self.OPS:
            raise ValueError(f"opcion debe ser una de {list(self.OPS)}.")
        if tol < 0.0:
            raise ValueError("tol debe ser >= 0.")
        self.param = [float(self.OPS[opcion]), float(tol)]
        self.a = Puerto(self, "ent", 0, 1)
        self.b = Puerto(self, "ent", 1, 1)
class LimitadorRapidez(Bloque):
    op = ops.OP_LIM_RAPIDEZ
    n_in = 1
    n_state = 1
    etiqueta = "LimitadorRapidez"
    def __init__(self, nombre, subida=1.0, bajada=1.0, valor_inicial=0.0):
        super().__init__(nombre)
        if subida < 0.0 or bajada < 0.0:
            raise ValueError("subida y bajada deben ser >= 0.")
        self.param = [float(subida), float(bajada)]
        self.estados_iniciales = [float(valor_inicial)]
class RetenedorDisparado(Bloque):
    op = ops.OP_RETENEDOR
    n_in = 2
    n_state = 2
    etiqueta = "RetenedorDisparado"
    def __init__(self, nombre, umbral=0.5, valor_inicial=0.0):
        super().__init__(nombre)
        self.param = [float(umbral)]
        self.estados_iniciales = [float(valor_inicial), -1.0]
        self.senal = Puerto(self, "ent", 0, 1)
        self.trigger = Puerto(self, "ent", 1, 1)
class MaquinaEstados(Bloque):
    op = ops.OP_MAQ_ESTADOS
    n_state = 1
    etiqueta = "MaquinaEstados"
    OPS = {"<": 0, "<=": 1, ">": 2, ">=": 3, "==": 4, "!=": 5}
    def __init__(self, nombre, n_estados, n_entradas, transiciones,
                 estado_inicial=0):
        super().__init__(nombre)
        if n_estados < 1:
            raise ValueError("n_estados debe ser >= 1.")
        if n_entradas < 1:
            raise ValueError("n_entradas debe ser >= 1.")
        if not (0 <= estado_inicial < n_estados):
            raise ValueError("estado_inicial fuera de rango.")
        self.n_in = n_entradas
        self.entradas = [Puerto(self, "ent", k, 1) for k in range(n_entradas)]
        plano = []
        for tr in transiciones:
            desde, hacia, sig, cond, umb = tr
            if not (0 <= desde < n_estados and 0 <= hacia < n_estados):
                raise ValueError("Transicion con estado fuera de rango.")
            if not (0 <= sig < n_entradas):
                raise ValueError("indice_senal fuera de rango.")
            if cond not in self.OPS:
                raise ValueError(f"condicion debe ser una de {list(self.OPS)}.")
            plano += [float(desde), float(hacia), float(sig),
                      float(self.OPS[cond]), float(umb)]
        self.param = [float(n_estados), float(len(transiciones))] + plano
        self.estados_iniciales = [float(estado_inicial)]
class FiltroPasoBajo(FuncionTransferencia):
    etiqueta = "FiltroPasoBajo"
    def __init__(self, nombre, fc, zeta=0.7071067811865476, orden=2):
        if fc <= 0.0:
            raise ValueError("fc debe ser > 0.")
        w = 2.0 * math.pi * fc
        if orden == 1:
            num, den = [w], [1.0, w]
        elif orden == 2:
            num, den = [w * w], [1.0, 2.0 * zeta * w, w * w]
        else:
            raise ValueError("orden debe ser 1 o 2.")
        super().__init__(nombre, num, den)
class FiltroPasoAlto(FuncionTransferencia):
    etiqueta = "FiltroPasoAlto"
    def __init__(self, nombre, fc, zeta=0.7071067811865476, orden=2):
        if fc <= 0.0:
            raise ValueError("fc debe ser > 0.")
        w = 2.0 * math.pi * fc
        if orden == 1:
            num, den = [1.0, 0.0], [1.0, w]
        elif orden == 2:
            num, den = [1.0, 0.0, 0.0], [1.0, 2.0 * zeta * w, w * w]
        else:
            raise ValueError("orden debe ser 1 o 2.")
        super().__init__(nombre, num, den)
class FiltroNotch(FuncionTransferencia):
    etiqueta = "FiltroNotch"
    def __init__(self, nombre, fn, zeta=0.3):
        if fn <= 0.0:
            raise ValueError("fn debe ser > 0.")
        w = 2.0 * math.pi * fn
        num, den = [1.0, 0.0, w * w], [1.0, 2.0 * zeta * w, w * w]
        super().__init__(nombre, num, den)
class MasaTermica(Bloque):
    op = ops.OP_MASA_TERMICA
    n_out = 1
    n_state = 1
    etiqueta = "MasaTermica"
    def __init__(self, nombre, C_th, *entradas, T_inicial=0.0, n_entradas=None,
                 T_amb=0.0, R_amb=0.0):
        super().__init__(nombre)
        if C_th <= 0.0:
            raise ValueError("C_th debe ser > 0.")
        if entradas:
            n_in = len(entradas)
        else:
            n_in = n_entradas if n_entradas is not None else 1
        if n_in < 1:
            raise ValueError("n_entradas debe ser >= 1.")
        self.n_in = n_in
        self.param = [float(C_th), float(T_amb), float(R_amb)]
        self.estados_iniciales = [float(T_inicial)]
        self.entradas = [Puerto(self, "ent", k, 1) for k in range(n_in)]
        if entradas:
            self._pendientes_conectar = [
                (e.puerto if hasattr(e, "puerto") else e, self.entradas[k])
                for k, e in enumerate(entradas)
            ]
class ResistenciaTermica(Bloque):
    op = ops.OP_RES_TERMICA
    n_in = 2
    n_out = 1
    etiqueta = "ResistenciaTermica"
    def __init__(self, nombre, R):
        super().__init__(nombre)
        if R <= 0.0:
            raise ValueError("R debe ser > 0.")
        self.param = [float(R)]
        self.T1 = Puerto(self, "ent", 0, 1)
        self.T2 = Puerto(self, "ent", 1, 1)
class Engranaje(Bloque):
    op = ops.OP_ENGRANAJE
    n_in = 2
    n_out = 2
    etiqueta = "Engranaje"
    def __init__(self, nombre, relacion=1.0):
        super().__init__(nombre)
        if relacion <= 0.0:
            raise ValueError("relacion debe ser > 0.")
        self.param = [float(relacion)]
        self.w1 = Puerto(self, "ent", 0, 1)
        self.T1 = Puerto(self, "ent", 1, 1)
class EjeFlexible(Bloque):
    op = ops.OP_EJE_FLEXIBLE
    n_in = 2
    n_out = 1
    n_state = 2
    etiqueta = "EjeFlexible"
    def __init__(self, nombre, K, B=0.0, theta1_0=0.0, theta2_0=0.0):
        super().__init__(nombre)
        if K < 0.0 or B < 0.0:
            raise ValueError("K y B deben ser >= 0.")
        self.param = [float(K), float(B)]
        self.estados_iniciales = [float(theta1_0), float(theta2_0)]
        self.w1 = Puerto(self, "ent", 0, 1)
        self.w2 = Puerto(self, "ent", 1, 1)
class Embrague(Bloque):
    op = ops.OP_EMBRAGUE
    n_in = 2
    n_out = 1
    etiqueta = "Embrague"
    def __init__(self, nombre, T_max, umbral=0.5):
        super().__init__(nombre)
        if T_max < 0.0:
            raise ValueError("T_max debe ser >= 0.")
        self.param = [float(T_max), float(umbral)]
        self.entrada = Puerto(self, "ent", 0, 1)
        self.control = Puerto(self, "ent", 1, 1)
def _leer_csv(archivo, columna_t, columna_y):
    import csv
    with open(archivo, newline="", encoding="utf-8-sig") as f:
        primera = f.readline()
        f.seek(0)
        delim = ";" if ";" in primera else ","
        pts = []
        for fila in csv.reader(f, delimiter=delim):
            if len(fila) <= max(columna_t, columna_y):
                continue
            try:
                t = float(fila[columna_t].strip())
                y = float(fila[columna_y].strip())
            except ValueError:
                continue
            pts.append((t, y))
    return pts
class FuenteCSV(Bloque):
    op = ops.OP_SRC_CSV
    n_out = 1
    etiqueta = "FuenteCSV"
    def __init__(self, nombre, archivo, columna_t=0, columna_y=1,
                 interpolar=True):
        super().__init__(nombre)
        pts = _leer_csv(archivo, columna_t, columna_y)
        if len(pts) < 2:
            raise ValueError(
                "El CSV debe tener al menos 2 puntos validos (t, y).")
        pts.sort(key=lambda p: p[0])
        for (t0, _), (t1, _) in zip(pts, pts[1:]):
            if t1 <= t0:
                raise ValueError("Los tiempos del CSV deben ser estrictamente "
                                 "crecientes.")
        ts = [t for t, _ in pts]
        ys = [y for _, y in pts]
        self.param = ([float(len(pts)), float(1.0 if interpolar else 0.0)]
                      + [float(t) for t in ts] + [float(y) for y in ys])
class FalloProgramado(Bloque):
    op = ops.OP_FALLO_PROG
    n_in = 1
    n_out = 1
    etiqueta = "FalloProgramado"
    def __init__(self, nombre, t_fallo, valor, modo=0):
        super().__init__(nombre)
        if t_fallo < 0.0:
            raise ValueError("t_fallo debe ser >= 0.")
        if modo not in (0, 1):
            raise ValueError("modo debe ser 0 (reemplazar) o 1 (sumar).")
        self.param = [float(t_fallo), float(valor), float(modo)]
class FalloEvento(Bloque):
    op = ops.OP_FALLO_EVENTO
    n_in = 2
    n_out = 1
    etiqueta = "FalloEvento"
    def __init__(self, nombre, umbral, valor, modo=0):
        super().__init__(nombre)
        if modo not in (0, 1):
            raise ValueError("modo debe ser 0 (reemplazar) o 1 (sumar).")
        self.param = [float(umbral), float(valor), float(modo)]
        self.senal = Puerto(self, "ent", 0, 1)
        self.disparo = Puerto(self, "ent", 1, 1)
class Multiplicador(Bloque):
    op = ops.OP_MULTIPLICADOR
    n_in = 2
    n_out = 1
    etiqueta = "Multiplicador"
    NOMBRES = ["y"]
    def __init__(self, nombre):
        super().__init__(nombre)
        self.entrada = Puerto(self, "ent", 0, 2)
        self.salida = Puerto(self, "sal", 0, 1)
class SaturarVectorial(Bloque):
    op = ops.OP_SAT_VECTORIAL
    n_in = 2
    n_out = 2
    n_state = 0
    etiqueta = "SaturarVectorial"
    NOMBRES = ["Vd_sat", "Vq_sat"]
    def __init__(self, nombre, Vmax=1.0):
        super().__init__(nombre)
        self.param = [float(Vmax)]
        self.entrada = Puerto(self, "ent", 0, 2)
        self.salida = Puerto(self, "sal", 0, 2)
class FuenteTabla(Bloque):
    op = ops.OP_SRC_TABLE
    n_in = 0
    n_out = 1
    etiqueta = "FuenteTabla"
    NOMBRES = ["y"]
    def __init__(self, nombre, puntos, interpolar=True):
        super().__init__(nombre)
        if len(puntos) < 2:
            raise ValueError("Necesita al menos 2 puntos (t, y).")
        pts = sorted(puntos, key=lambda p: p[0])
        for (t0, _), (t1, _) in zip(pts, pts[1:]):
            if t1 <= t0:
                raise ValueError("Tiempos deben ser estrictamente crecientes.")
        ts = [float(t) for t, _ in pts]
        ys = [float(y) for _, y in pts]
        n_pts = len(pts)
        self.param = [float(n_pts), float(1.0 if interpolar else 0.0)] + ts + ys
        self.etiqueta = f"FuenteTabla({n_pts} pts)"
class BateriaECM(Bloque):
    op = ops.OP_BATERIA_ECM
    n_in = 1
    n_out = 6
    n_state = 4
    etiqueta = "BateriaECM"
    NOMBRES = ["V_term", "SOC", "T_cell", "P_loss", "I_chg_lim", "I_dch_lim"]
    def __init__(self, nombre,
                 Q_nom=60.0, V_nom=3.65, V_min_cell=2.5, V_max_cell=4.2,
                 R0=0.8e-3, R1=0.5e-3, C1=8000.0, R2=0.3e-3, C2=80000.0,
                 N_series=96, N_parallel=2,
                 ocv_soc=None,
                 I_chg_cont=150.0, I_dch_cont=400.0,
                 T_min_chg=0.0, T_max_chg=45.0, T_min_dch=-20.0, T_max_dch=60.0,
                 R_th_pack=0.3, C_th_pack=10000.0, T_amb=25.0,
                 deg_per_cycle=5e-5,
                 soc_init=0.9):
        super().__init__(nombre)
        if ocv_soc is None:
            ocv_soc = [
                [0.00, 2.80], [0.05, 3.25], [0.10, 3.45], [0.15, 3.52],
                [0.20, 3.58], [0.30, 3.65], [0.40, 3.71], [0.50, 3.75],
                [0.60, 3.80], [0.70, 3.87], [0.80, 3.96], [0.90, 4.05],
                [0.95, 4.12], [1.00, 4.20]
            ]
        n_ocv = len(ocv_soc)
        ocv_flat = [float(v) for pt in ocv_soc for v in pt]
        self.param = [
            float(Q_nom), float(V_nom), float(V_min_cell), float(V_max_cell),
            float(R0), float(R1), float(C1), float(R2), float(C2),
            float(N_series), float(N_parallel),
            float(n_ocv)] + ocv_flat + [
            float(I_chg_cont), float(I_dch_cont),
            float(T_min_chg), float(T_max_chg), float(T_min_dch), float(T_max_dch),
            float(R_th_pack), float(C_th_pack), float(T_amb),
            float(deg_per_cycle)
        ]
        self.estados_iniciales = [float(soc_init), 0.0, 0.0, float(T_amb)]
        self.entrada = Puerto(self, "ent", 0, 1)
        self.salida = Puerto(self, "sal", 0, 6)
class CalculoIdc(Bloque):
    op = ops.OP_CALCULO_IDC
    n_in = 5
    n_out = 1
    etiqueta = "CalculoIdc"
    NOMBRES = ["I_dc"]
    def __init__(self, nombre, eff=1.0):
        super().__init__(nombre)
        self.param = [float(eff)]
        self.entrada = Puerto(self, "ent", 0, 5)
        self.salida = Puerto(self, "sal", 0, 1)
class Vehiculo(Bloque):
    op = ops.OP_VEHICULO
    n_in = 3
    n_out = 7
    n_state = 1
    estados_iniciales = [0.0]
    opcionales = 1
    etiqueta = "Vehiculo"
    NOMBRES = ["T_load", "v_veh", "v_kmh", "omega_ref", "T_ff", "grade",
               "omega_m_real"]
    def __init__(self, nombre,
                 mass=1800.0, Cd=0.28, A_frontal=2.3, Crr=0.01,
                 rho_air=1.225, g=9.81,
                 gear_ratio_total=8.0, wheel_radius=0.33,
                 gear_efficiency=0.95, regen_eff=0.7,
                 ref_timescale=10.0):
        super().__init__(nombre)
        self.param = [
            float(mass), float(Cd), float(A_frontal), float(Crr),
            float(rho_air), float(g),
            float(gear_ratio_total), float(wheel_radius),
            float(gear_efficiency), float(regen_eff),
            float(ref_timescale)
        ]
        self.entrada = Puerto(self, "ent", 0, 3)
        self.salida = Puerto(self, "sal", 0, 7)
