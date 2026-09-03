from . import opcodes as ops
from .bloques import Bloque
from .puertos import Puerto, Sensor
class GeneradorPWM(Bloque):
    op = ops.OP_PWM_1F
    n_in = 1
    n_out = 1
    n_state = 0
    etiqueta = "PWM"
    def __init__(self, nombre, fsw=10_000.0, dead_time=0.0):
        super().__init__(nombre)
        self.param = [float(fsw), float(dead_time)]
        self.etiqueta = f"PWM(fsw={fsw})"
        self.d = Puerto(self, "ent", 0, 1)
        self.salida = Puerto(self, "sal", 0, 1)
    def sensorDisparo(self):
        return Sensor("S_pwm", self, "sal", 0, 1, canales=["S"])
class GeneradorSPWM(Bloque):
    op = ops.OP_PWM_SPWM
    n_in = 1
    n_out = 3
    n_state = 0
    etiqueta = "SPWM"
    NOMBRES = ["Sa", "Sb", "Sc"]
    def __init__(self, nombre, f_out=50.0, fsw=10_000.0,
                 fase_ini=0.0, dead_time=0.0):
        super().__init__(nombre)
        self.param = [float(f_out), float(fsw), float(fase_ini), float(dead_time)]
        self.etiqueta = f"SPWM(f={f_out},fsw={fsw})"
        self.m = Puerto(self, "ent", 0, 1)
        self.salida = Puerto(self, "sal", 0, 3)
    def sensorDisparos(self):
        return Sensor("S_spwm", self, "sal", 0, 3, canales=["Sa", "Sb", "Sc"])
class GeneradorSVPWM(Bloque):
    op = ops.OP_PWM_SVPWM
    n_in = 2
    n_out = 3
    n_state = 0
    etiqueta = "SVPWM"
    NOMBRES = ["Sa", "Sb", "Sc"]
    def __init__(self, nombre, Vdc=600.0, fsw=10_000.0, dead_time=0.0):
        super().__init__(nombre)
        self.param = [float(Vdc), float(fsw), float(dead_time)]
        self.etiqueta = f"SVPWM(Vdc={Vdc},fsw={fsw})"
        self.v_alpha = Puerto(self, "ent", 0, 1)
        self.v_beta = Puerto(self, "ent", 1, 1)
        self.salida = Puerto(self, "sal", 0, 3)
    def sensorDisparos(self):
        return Sensor("S_svpwm", self, "sal", 0, 3, canales=["Sa", "Sb", "Sc"])
