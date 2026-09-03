from . import opcodes as ops
from .bloques import Bloque
from .puertos import Puerto, Sensor
class MotorHardwareCH32(Bloque):
    op = ops.OP_HW_SERIAL
    n_in = 1
    n_out = 2
    n_state = 0
    etiqueta = "HW_MOTOR"
    NOMBRES = ["angulo", "rpm"]
    def __init__(self, nombre, puerto_com, baudrate=115200):
        super().__init__(nombre)
        self.param = [float(puerto_com), float(baudrate)]
        self.etiqueta = f"CH32(COM{puerto_com})"
        self.duty = Puerto(self, "ent", 0, 1)
        self.salida_angulo = Puerto(self, "sal", 0, 1)
        self.salida_rpm = Puerto(self, "sal", 1, 1)
    def sensorAngulo(self):
        return Sensor(f"{self.nombre}_ang", self, "sal", 0, 1)
    def sensorVelocidad(self):
        return Sensor(f"{self.nombre}_rpm", self, "sal", 1, 1)
