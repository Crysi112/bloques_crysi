import numpy as np
from . import opcodes as ops
from .bloques import Bloque
from .puertos import Puerto, Sensor
class Maquina(Bloque):
    n_in = 0
    n_out = 10
    n_state = 0
    NOMBRES = ["ia", "ib", "ic", "iqs", "ids", "wm", "th_rm", "th_e", "Te", "P_cu_s"]
    def __init__(self, nombre):
        super().__init__(nombre)
        self.terminales = Puerto(self, "ent", 0, 3)
        self.salida = Puerto(self, "sal", 0, self.n_out)
    def sensor3V(self):
        return Sensor("V", self, "ent", 0, 3,
                      canales=["va", "vb", "vc"])
    def sensor3I(self):
        return Sensor("I", self, "sal", 0, 3,
                      canales=["ia", "ib", "ic"])
    def sensorVelocidad(self):
        return Sensor("wm", self, "sal", 5, 1)
    def sensorPosicion(self):
        return Sensor("th_rm", self, "sal", 6, 1)
    def sensorPosicionElectrica(self):
        return Sensor("th_e", self, "sal", 7, 1)
    def sensorPar(self):
        return Sensor("Te", self, "sal", 8, 1)
    def sensorCorrienteD(self):
        return Sensor("ids", self, "sal", 4, 1)
    def sensorCorrienteQ(self):
        return Sensor("iqs", self, "sal", 3, 1)
    def sensorPerdidasEstator(self):
        return Sensor("P_cu_s", self, "sal", 9, 1)
    def resumen(self, res: "Resultado") -> str:
        import numpy as np
        out = []
        out.append("=" * 60)
        out.append(f" RESUMEN DINÁMICO: {self.nombre} ({self.etiqueta})")
        out.append("=" * 60)
        try:
            Te_pico = res.pico("Te", t_max=0.1)
            out.append(f"  • Par electromagnético pico (Te) : {Te_pico:.2f} N·m")
        except Exception:
            pass
        try:
            I_pico = res.pico("I", t_max=0.1)
            out.append(f"  • Corriente de estator pico (Ia) : {I_pico:.2f} A")
        except Exception:
            pass
        try:
            ws = self.velocidad_sincronica
            wm_final = res.final("wm")
            out.append(f"  • Velocidad final (wm) : {wm_final:.2f} rad/s ({wm_final*60/(2*np.pi):.1f} rpm)")
            out.append(f"  • Deslizamiento nominal (s) : {100*(ws-wm_final)/ws:.2f} %")
        except Exception:
            pass
        try:
            ts = res.tiempo_establecimiento("wm", 0.02)
            out.append(f"  • Tiempo de arranque (ts, 98%) : {ts:.3f} s")
        except Exception:
            pass
        try:
            Pcu_final = res.final("P_cu_s")
            out.append(f"  • Pérdidas Joule estator (Pcu) : {Pcu_final:.1f} W")
        except Exception:
            pass
        out.append("=" * 60)
        return "\n".join(out)
class MaquinaImanesPermanentes(Maquina):
    op = ops.OP_MAQ_PMAC
    n_in = 4
    n_state = 4
    n_out = 10
    etiqueta = "PMAC"
    def __init__(self, nombre, rs, Ld, Lq, lam_m, P, J, Bm=0.0,
                 th_inicial=0.0, mecanica_interna=True, saturacion=None):
        super().__init__(nombre)
        self.mecanica_interna = bool(mecanica_interna)
        self.param = [float(rs), float(Ld), float(Lq), float(lam_m),
                      float(P), float(J), float(Bm), 0.0 if mecanica_interna else 1.0]
        if saturacion is not None:
            pts = [(float(a), float(b)) for a, b in saturacion]
            if len(pts) < 2:
                raise ValueError("saturacion necesita al menos 2 puntos.")
            ids_lut = [a for a, _ in pts]
            if any(ids_lut[i] >= ids_lut[i + 1]
                   for i in range(len(ids_lut) - 1)):
                raise ValueError("los puntos de saturacion deben tener Id creciente.")
            self.param += [float(len(pts))] + \
                [v for pt in pts for v in pt]
            self.saturacion = pts
        self.estados_iniciales = [0.0, 0.0, 0.0, float(th_inicial)]
        self.etiqueta = f"PMAC(P={P})"
        if mecanica_interna:
            self.n_in = 4
            self.T_L = Puerto(self, "ent", 3, 1)
        else:
            self.n_in = 5
            self.puerto_mecanico = Puerto(self, "ent", 3, 2)
            self.entrada = Puerto(self, "ent", 0, self.n_in)
class MaquinaInduccion(Maquina):
    op = ops.OP_MAQ_INDUCCION
    n_in = 4
    n_state = 6
    n_out = 13
    etiqueta = "MaquinaMI"
    def __init__(self, nombre, rs, rr, Lm, Lls, Llr, P, J, Bm=0.0, w_frame=0.0, mecanica_interna=True):
        super().__init__(nombre)
        self.mecanica_interna = bool(mecanica_interna)
        Ls = Lls + Lm
        Lr = Llr + Lm
        det = Ls * Lr - Lm * Lm
        Li00, Li01, Li11 = Lr / det, -Lm / det, Ls / det
        self.param = [float(rs), float(rr), float(Li00), float(Li01), float(Li11),
                      float(w_frame), float(P), float(J), float(Bm), 0.0 if mecanica_interna else 1.0]
        self.estados_iniciales = [0.0] * 6
        self.etiqueta = f"MI(P={P})"
        if mecanica_interna:
            self.n_in = 4
            self.T_L = Puerto(self, "ent", 3, 1)
        else:
            self.n_in = 5
            self.puerto_mecanico = Puerto(self, "ent", 3, 2)
            self.entrada = Puerto(self, "ent", 0, self.n_in)
    def sensorCorrienteRotor(self):
        return Sensor("I_rotor", self, "sal", 10, 3,
                      canales=["i'ar", "i'br", "i'cr"])
    @property
    def velocidad_sincronica(self):
        P = self.param[6] if len(self.param) > 6 else 4
        f = 60.0
        return 4.0 * np.pi * f / P
class MaquinaSincrona(Maquina):
    op = ops.OP_MAQ_SINCRONA
    n_in = 5
    n_state = 8
    n_out = 10
    etiqueta = "MaquinaSincrona"
    def __init__(self, nombre, rs, rfd, rkq1, rkq2, rkd,
                 Lls, Lmq, Llkq1, Llkq2, Lmd, Llf, Llkd,
                 P, J, Bm=0.0, th_inicial=0.0, mecanica_interna=True):
        super().__init__(nombre)
        self.mecanica_interna = bool(mecanica_interna)
        Lq = np.array([
            [Lls + Lmq, Lmq, Lmq],
            [Lmq, Llkq1 + Lmq, Lmq],
            [Lmq, Lmq, Llkq2 + Lmq],
        ])
        Ld = np.array([
            [Lls + Lmd, Lmd, Lmd],
            [Lmd, Llf + Lmd, Lmd],
            [Lmd, Lmd, Llkd + Lmd],
        ])
        Liq = np.linalg.inv(Lq).ravel()
        Lid = np.linalg.inv(Ld).ravel()
        self.param = ([float(rs), float(rfd), float(rkq1), float(rkq2), float(rkd),
                       float(P), float(J), float(Bm)]
                      + [float(x) for x in Liq] + [float(x) for x in Lid] + [0.0 if mecanica_interna else 1.0])
        self.estados_iniciales = [0.0] * 7 + [float(th_inicial)]
        self.etiqueta = f"MS(P={P})"
        self.vfd = Puerto(self, "ent", 3, 1)
        if mecanica_interna:
            self.n_in = 5
            self.T_L = Puerto(self, "ent", 4, 1)
        else:
            self.n_in = 6
            self.puerto_mecanico = Puerto(self, "ent", 4, 2)
            self.entrada = Puerto(self, "ent", 0, self.n_in)
    def sensorVoltajeCampo(self):
        return Sensor("v_fd", self, "ent", 3, 1)
class MaquinaCorrienteContinua(Maquina):
    op = ops.OP_MAQ_CC
    n_in = 3
    n_out = 8
    n_state = 4
    etiqueta = "MaquinaCC"
    NOMBRES = ["ia", "if", "wm", "th_rm", "Te", "Ea", "V_t", "P_cu"]
    def __init__(self, nombre, r_a, L_a, r_f, L_f, L_AF, J, Bm=0.0, mecanica_interna=True):
        super().__init__(nombre)
        self.mecanica_interna = bool(mecanica_interna)
        self.n_out = 8
        self.param = [float(r_a), float(L_a), float(r_f), float(L_f),
                      float(L_AF), float(J), float(Bm), 0.0 if mecanica_interna else 1.0]
        self.estados_iniciales = [0.0] * 4
        self.entrada = Puerto(self, "ent", 0, 1)
        self.campo = Puerto(self, "ent", 1, 1)
        if mecanica_interna:
            self.n_in = 3
            self.T_L = Puerto(self, "ent", 2, 1)
        else:
            self.n_in = 4
            self.puerto_mecanico = Puerto(self, "ent", 2, 2)
        self.salida = Puerto(self, "sal", 0, 8)
    def sensorCorriente(self):
        return Sensor("ia", self, "sal", 0, 1)
    def sensorCampo(self):
        return Sensor("if", self, "sal", 1, 1)
    def sensorVelocidad(self):
        return Sensor("wm", self, "sal", 2, 1)
    def sensorPosicion(self):
        return Sensor("th_rm", self, "sal", 3, 1)
    def sensorPar(self):
        return Sensor("Te", self, "sal", 4, 1)
    def sensorEa(self):
        return Sensor("Ea", self, "sal", 5, 1)
    def sensorVoltajeTerminal(self):
        return Sensor("V_t", self, "sal", 6, 1)
    def sensorPerdidas(self):
        return Sensor("P_cu", self, "sal", 7, 1)
class MaquinaDCImanesPermanentes(Maquina):
    op = ops.OP_MAQ_DC_PM
    n_in = 2
    n_out = 7
    n_state = 3
    etiqueta = "MaquinaDCPM"
    NOMBRES = ["ia", "wm", "th_rm", "Te", "Ea", "V_t", "P_cu"]
    def __init__(self, nombre, r_a, L_a, Kt, J, Bm=0.0, mecanica_interna=True):
        super().__init__(nombre)
        self.n_out = 7
        self.mecanica_interna = bool(mecanica_interna)
        self.param = [float(r_a), float(L_a), float(Kt),
                      float(J), float(Bm), 0.0 if mecanica_interna else 1.0]
        self.estados_iniciales = [0.0] * 3
        self.entrada = Puerto(self, "ent", 0, 1)
        if mecanica_interna:
            self.n_in = 2
            self.T_L = Puerto(self, "ent", 1, 1)
        else:
            self.n_in = 3
            self.puerto_mecanico = Puerto(self, "ent", 1, 2)
        self.salida = Puerto(self, "sal", 0, 7)
    def sensorCorriente(self):
        return Sensor("ia", self, "sal", 0, 1)
    def sensorVelocidad(self):
        return Sensor("wm", self, "sal", 1, 1)
    def sensorPosicion(self):
        return Sensor("th_rm", self, "sal", 2, 1)
    def sensorPar(self):
        return Sensor("Te", self, "sal", 3, 1)
    def sensorEa(self):
        return Sensor("Ea", self, "sal", 4, 1)
    def sensorVoltajeTerminal(self):
        return Sensor("V_t", self, "sal", 5, 1)
    def sensorPerdidas(self):
        return Sensor("P_cu", self, "sal", 6, 1)
    def sensorVoltajeTerminal(self):
        return Sensor("V_t", self, "sal", 5, 1)
