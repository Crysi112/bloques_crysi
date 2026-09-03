import numpy as np
from . import opcodes as ops
from .bloques import Bloque
from .puertos import Puerto
class ConvertidorDC_DC(Bloque):
    n_in = 2
    n_out = 2
    n_state = 2
    NOMBRES = ["vout", "iL"]
    def __init__(self, nombre, L, C, R):
        super().__init__(nombre)
        self.param = [float(L), float(C), float(R)]
        self.etiqueta = f"{self.etiqueta}(L={L}, C={C}, R={R})"
        self.entrada = Puerto(self, "ent", 0, 1)
        self.d = Puerto(self, "ent", 1, 1)
        self.salida = Puerto(self, "sal", 0, 2)
class ConvertidorBuck(ConvertidorDC_DC):
    op = ops.OP_POT_BUCK
    etiqueta = "ConvertidorBuck"
class ConvertidorBoost(ConvertidorDC_DC):
    op = ops.OP_POT_BOOST
    etiqueta = "ConvertidorBoost"
class ConvertidorBuckBoost(ConvertidorDC_DC):
    op = ops.OP_POT_BUCKBOOST
    etiqueta = "ConvertidorBuckBoost"
class RectificadorTrifasico(Bloque):
    op = ops.OP_POT_RECT_3F
    n_in = 3
    n_out = 2
    n_state = 1
    etiqueta = "Rectificador3F"
    NOMBRES = ["vdc", "idc"]
    def __init__(self, nombre, C, R, Rint=1e-3):
        super().__init__(nombre)
        self.param = [float(C), float(R), float(Rint)]
        self.etiqueta = f"Rect3F(C={C}, R={R})"
        self.entrada = Puerto(self, "ent", 0, 3)
        self.salida = Puerto(self, "sal", 0, 2)
class InversorTrifasico(Bloque):
    op = ops.OP_POT_INV_3F
    n_in = 1
    n_out = 6
    n_state = 6
    etiqueta = "Inversor3F"
    NOMBRES = ["vCa", "vCb", "vCc", "iLa", "iLb", "iLc"]
    def __init__(self, nombre, f_out=60.0, fsw=5000.0,
                 m_start=0.4, m_end=1.0, t_ramp=0.0,
                 Lf=5e-3, Cf=20e-6, R=50.0, conmutada=False):
        super().__init__(nombre)
        self.param = [float(f_out), float(fsw), float(m_start), float(m_end),
                      float(t_ramp), float(Lf), float(Cf), float(R),
                      float(1 if conmutada else 0)]
        self.etiqueta = (f"Inv3F(f_out={f_out}, fsw={fsw}, "
                         f"m={m_start}->{m_end})")
        self.entrada = Puerto(self, "ent", 0, 1)
        self.salida = Puerto(self, "sal", 0, 6)
    def sensorVoltajesSalida(self):
        from .puertos import Sensor
        return Sensor("V_f", self, "sal", 0, 3, canales=["vCa", "vCb", "vCc"])
    def sensorCorrientesFase(self):
        from .puertos import Sensor
        return Sensor("I_f", self, "sal", 3, 3, canales=["iLa", "iLb", "iLc"])
class InversorMonofasico(Bloque):
    op = ops.OP_POT_INV_1F
    n_in = 1
    n_out = 2
    n_state = 2
    etiqueta = "Inversor1F"
    NOMBRES = ["vC", "iL"]
    def __init__(self, nombre, f_out=50.0, fsw=10000.0,
                 m_start=0.4, m_end=1.0, t_ramp=0.0,
                 Lf=5e-3, Cf=20e-6, R=50.0, conmutada=False):
        super().__init__(nombre)
        self.param = [float(f_out), float(fsw), float(m_start), float(m_end),
                      float(t_ramp), float(Lf), float(Cf), float(R),
                      float(1 if conmutada else 0)]
        self.etiqueta = (f"Inv1F(f_out={f_out}, fsw={fsw}, "
                         f"m={m_start}->{m_end})")
class CargaRLTrifasica(Bloque):
    op = ops.OP_CARGA_RL_3F
    n_in = 3
    n_out = 3
    n_state = 2
    etiqueta = "CargaRL3F"
    NOMBRES = ["ia", "ib", "ic"]
    def __init__(self, nombre, R=1.0, L=1e-3):
        super().__init__(nombre)
        self.param = [float(R), float(L)]
    @classmethod
    def desde_pq(cls, nombre: str, p_w: float, q_var: float, v_ll: float = 400.0, f: float = 50.0):
        v_ln = v_ll / np.sqrt(3.0)
        p = float(p_w)
        q = float(q_var)
        den = p * p + q * q
        if den <= 0:
            raise ValueError("P y Q no pueden ser simultáneamente cero.")
        r = 3.0 * (v_ln ** 2) * p / den
        x = 3.0 * (v_ln ** 2) * q / den
        l = x / (2.0 * np.pi * f) if f > 0 else 0.0
        return cls(nombre, R=r, L=l)
class CargaPQTrifasica(Bloque):
    op = ops.OP_CARGA_PQ_3F
    n_in = 3
    n_out = 3
    n_state = 0
    etiqueta = "CargaPQ3F"
    NOMBRES = ["ia", "ib", "ic"]
    def __init__(self, nombre, p_w: float = 1000.0, q_var: float = 500.0):
        super().__init__(nombre)
        self.param = [float(p_w), float(q_var)]
class CargaPQMonofasica(Bloque):
    op = ops.OP_CARGA_PQ_1F
    n_in = 1
    n_out = 1
    n_state = 0
    etiqueta = "CargaPQ1F"
    NOMBRES = ["i"]
    def __init__(self, nombre, p_w: float = 1000.0, q_var: float = 500.0):
        super().__init__(nombre)
        self.param = [float(p_w), float(q_var)]
class Transformador(Bloque):
    op = ops.OP_TRANSFORMADOR
    etiqueta = "Transformador"
    def __init__(self, nombre, a=1.0, fases=1):
        super().__init__(nombre)
        if a <= 0:
            raise ValueError("a debe ser > 0.")
        self.n_in = 2 * fases
        self.n_out = 2 * fases
        self.param = [float(a)]
        self.entrada = Puerto(self, "ent", 0, self.n_in)
        self.salida = Puerto(self, "sal", 0, self.n_out)
        if fases == 3:
            self.NOMBRES = ["va2", "vb2", "vc2", "ia1", "ib1", "ic1"]
        else:
            self.NOMBRES = ["v2", "i1"]
class Resistencia(Bloque):
    op = ops.OP_RESISTENCIA
    n_in = 1
    n_out = 1
    n_state = 0
    etiqueta = "Resistencia"
    NOMBRES = ["i"]
    def __init__(self, nombre, R=1.0):
        super().__init__(nombre)
        if R <= 0:
            raise ValueError("R debe ser > 0.")
        self.param = [float(R)]
        self.etiqueta = f"Resistencia(R={R})"
class Inductor(Bloque):
    op = ops.OP_INDUCTOR
    n_in = 1
    n_out = 1
    n_state = 1
    etiqueta = "Inductor"
    NOMBRES = ["i"]
    def __init__(self, nombre, L=1e-3):
        super().__init__(nombre)
        if L <= 0:
            raise ValueError("L debe ser > 0.")
        self.param = [float(L)]
        self.etiqueta = f"Inductor(L={L})"
class Capacitor(Bloque):
    op = ops.OP_CAPACITOR
    n_in = 1
    n_out = 1
    n_state = 1
    etiqueta = "Capacitor"
    NOMBRES = ["v"]
    def __init__(self, nombre, C=1e-3):
        super().__init__(nombre)
        if C <= 0:
            raise ValueError("C debe ser > 0.")
        self.param = [float(C)]
        self.etiqueta = f"Capacitor(C={C})"
