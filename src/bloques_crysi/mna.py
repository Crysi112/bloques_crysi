import numpy as np
from typing import List, Tuple, Dict, Any, Optional

from .bloques import Bloque
from .puertos import Puerto
from . import opcodes as ops

class Nodo:
    def __init__(self, nombre: str):
        self.nombre = nombre
        self.idx = -1

class Componente:
    def __init__(self, nombre: str, n1: Nodo, n2: Nodo):
        self.nombre = nombre
        self.n1 = n1
        self.n2 = n2

class Resistor(Componente):
    def __init__(self, nombre: str, n1: Nodo, n2: Nodo, R: float):
        super().__init__(nombre, n1, n2)
        self.R = float(R)

class Capacitor(Componente):
    def __init__(self, nombre: str, n1: Nodo, n2: Nodo, C: float):
        super().__init__(nombre, n1, n2)
        self.C = float(C)

class Inductor(Componente):
    def __init__(self, nombre: str, n1: Nodo, n2: Nodo, L: float):
        super().__init__(nombre, n1, n2)
        self.L = float(L)

class VSource(Componente):
    def __init__(self, nombre: str, n1: Nodo, n2: Nodo, idx_u: int):
        super().__init__(nombre, n1, n2)
        self.idx_u = int(idx_u)
        self.idx_i = -1

class ISource(Componente):
    def __init__(self, nombre: str, n1: Nodo, n2: Nodo, idx_u: int):
        super().__init__(nombre, n1, n2)
        self.idx_u = int(idx_u)

class Switch(Componente):
    def __init__(self, nombre: str, n1: Nodo, n2: Nodo, idx_ctrl: int, Ron=1e-3, Roff=1e6):
        super().__init__(nombre, n1, n2)
        self.idx_ctrl = int(idx_ctrl)
        self.Ron = float(Ron)
        self.Roff = float(Roff)

class Diodo(Componente):
    def __init__(self, nombre: str, n1: Nodo, n2: Nodo, Ron=1e-3, Roff=1e6, Vf=0.7):
        super().__init__(nombre, n1, n2)
        self.Ron = float(Ron)
        self.Roff = float(Roff)
        self.Vf = float(Vf)

class VCVS(Componente):
    def __init__(self, nombre: str, n1: Nodo, n2: Nodo, idx_u: int, gain: float = 1.0):
        super().__init__(nombre, n1, n2)
        self.idx_u = int(idx_u)
        self.gain = float(gain)
        self.idx_i = -1

class VCCS(Componente):
    def __init__(self, nombre: str, n1: Nodo, n2: Nodo, idx_u: int, gm: float = 1.0):
        super().__init__(nombre, n1, n2)
        self.idx_u = int(idx_u)
        self.gm = float(gm)

class MutualInductor(Componente):
    def __init__(self, nombre: str, n1: Nodo, n2: Nodo, n3: Nodo, n4: Nodo, L1: float, L2: float, M: float):
        super().__init__(nombre, n1, n2)
        self.n3 = n3
        self.n4 = n4
        self.L1 = float(L1)
        self.L2 = float(L2)
        self.M = float(M)
        self.idx_i1 = -1
        self.idx_i2 = -1

class SubredMNA(Bloque):
    op = 74

    def __init__(self, nombre: str, nodos: List[Nodo], componentes: List[Componente], dt: float,
                 mediciones_v: List[Tuple[Nodo, Nodo]],
                 mediciones_i: List[VSource],
                 modo: str = "lu", precomputar: bool = False, metodo: int = 0):
        super().__init__(nombre)
        ops.OP_MNA = 74
        ops.ES_DINAMICO.add(74)
        self.op = 74
        self.nodos = []
        for n in nodos:
            if str(n.nombre).lower() in ["0", "gnd", "tierra"]:
                n.idx = -1
            else:
                n.idx = len(self.nodos)
                self.nodos.append(n)
        self.num_nodos = len(self.nodos)
        self.componentes = componentes
        self.dt = float(dt)
        self.metodo = int(metodo)
        self.mediciones_v = mediciones_v
        self.mediciones_i = mediciones_i
        self.modo = modo
        self.resistors = [c for c in componentes if isinstance(c, Resistor)]
        self.capacitors = [c for c in componentes if isinstance(c, Capacitor)]
        self.inductors = [c for c in componentes if isinstance(c, Inductor)]
        self.vsources = [c for c in componentes if isinstance(c, VSource)]
        self.isources = [c for c in componentes if isinstance(c, ISource)]
        self.switches = [c for c in componentes if isinstance(c, Switch)]
        self.diodos = [c for c in componentes if isinstance(c, Diodo)]
        self.vcvss = [c for c in componentes if isinstance(c, VCVS)]
        self.vccss = [c for c in componentes if isinstance(c, VCCS)]
        self.mutuals = [c for c in componentes if isinstance(c, MutualInductor)]

        idx_ext = self.num_nodos
        for vs in self.vsources + self.inductors + self.vcvss:
            vs.idx_i = idx_ext
            idx_ext += 1
        for mut in self.mutuals:
            mut.idx_i1 = idx_ext; idx_ext += 1
            mut.idx_i2 = idx_ext; idx_ext += 1

        self.n_x = idx_ext
        idx_u_max = -1
        for src in self.vsources + self.isources:
            if src.idx_u > idx_u_max: idx_u_max = src.idx_u
        for sw in self.switches:
            if sw.idx_ctrl > idx_u_max: idx_u_max = sw.idx_ctrl
        self.n_u = max(idx_u_max + 1, 1)
        self.n_out = len(mediciones_v) + len(mediciones_i)
        self.n_in = self.n_u
        self.n_state = self.n_x
        self.entrada = Puerto(self, "ent", 0, self.n_in)
        self.salida = Puerto(self, "sal", 0, self.n_out)
        self.etiqueta = f"MNA({self.n_x}x{self.n_x},{'tustin' if self.metodo==1 else 'be'})"

        n_sw_ctrl = len(self.switches)
        n_diodos = len(self.diodos)
        self.n_sw = n_sw_ctrl + n_diodos

        self.param = [float(self.n_x), float(self.n_u), float(n_sw_ctrl), float(self.n_out), float(n_diodos), float(len(self.mediciones_v))]
        for n1, n2 in self.mediciones_v:
            self.param.extend([float(n1.idx), float(n2.idx)])
        self.param.append(float(len(self.mediciones_i)))
        for vsrc in self.mediciones_i:
            self.param.append(float(vsrc.idx_i))
        for sw in self.switches:
            self.param.append(float(sw.idx_ctrl))
        for d in self.diodos:
            self.param.extend([float(d.n1.idx), float(d.n2.idx), float(d.Vf), float(d.Ron), float(d.Roff)])

        self.param.append(float(len(self.resistors)))
        for r in self.resistors: self.param.extend([float(r.n1.idx), float(r.n2.idx), float(r.R)])
        self.param.append(float(len(self.capacitors)))
        for c in self.capacitors: self.param.extend([float(c.n1.idx), float(c.n2.idx), float(c.C)])
        self.param.append(float(len(self.inductors)))
        for ind in self.inductors: self.param.extend([float(ind.n1.idx), float(ind.n2.idx), float(ind.idx_i), float(ind.L)])
        self.param.append(float(len(self.vsources)))
        for vs in self.vsources: self.param.extend([float(vs.n1.idx), float(vs.n2.idx), float(vs.idx_u), float(vs.idx_i)])
        self.param.append(float(len(self.isources)))
        for isrc in self.isources: self.param.extend([float(isrc.n1.idx), float(isrc.n2.idx), float(isrc.idx_u)])
        self.param.append(float(len(self.switches)))
        for sw in self.switches: self.param.extend([float(sw.n1.idx), float(sw.n2.idx), float(sw.Ron), float(sw.Roff)])
        self.param.append(float(len(self.vcvss)))
        for vcv in self.vcvss: self.param.extend([float(vcv.n1.idx), float(vcv.n2.idx), float(vcv.idx_u), float(vcv.idx_i), float(vcv.gain)])
        self.param.append(float(len(self.vccss)))
        for vcc in self.vccss: self.param.extend([float(vcc.n1.idx), float(vcc.n2.idx), float(vcc.idx_u), float(vcc.gm)])
        self.param.append(float(len(self.mutuals)))
        for mut in self.mutuals: self.param.extend([float(mut.n1.idx), float(mut.n2.idx), float(mut.n3.idx), float(mut.n4.idx), float(mut.idx_i1), float(mut.idx_i2), float(mut.L1), float(mut.L2), float(mut.M)])
        self.param.append(float(self.metodo))

        self.matrices_Bu = []
        self.matrices_Bx = []
        if precomputar:
            if self.metodo == 1:
                raise ValueError("Precomputado solo BE")
            for estado_idx in range(2**self.n_sw):
                bits = [(estado_idx >> i) & 1 for i in range(self.n_sw)]
                Bu, Bx = self._construir_matrices(bits)
                self.matrices_Bu.append(Bu)
                self.matrices_Bx.append(Bx)
            self.param.extend([1.0, float(2**self.n_sw)])
            for idx in range(2**self.n_sw):
                self.param.extend(self.matrices_Bx[idx].flatten().tolist())
                self.param.extend(self.matrices_Bu[idx].flatten().tolist())
            self.etiqueta = f"MNA({self.n_x}x{self.n_x},pre={2**self.n_sw})"
        else:
            self.param.append(0.0)
            self.etiqueta = f"MNA({self.n_x}x{self.n_x},lu)"

        self.estados_iniciales = [0.0] * self.n_x
        self.NOMBRES = [f"V_{n1.nombre}_{n2.nombre}" for n1, n2 in self.mediciones_v] + [f"I_{v.nombre}" for v in self.mediciones_i]

    def _construir_matrices(self, estados_sw: List[int]):
        G = np.zeros((self.n_x, self.n_x))
        C = np.zeros((self.n_x, self.n_x))
        W = np.zeros((self.n_x, self.n_u))

        def stamp_admitance(mat, n1, n2, val):
            if n1 >= 0: mat[n1, n1] += val
            if n2 >= 0: mat[n2, n2] += val
            if n1 >= 0 and n2 >= 0:
                mat[n1, n2] -= val; mat[n2, n1] -= val
        def stamp_vsrc(matG, n1, n2, idx_i):
            if n1 >= 0: matG[n1, idx_i] += 1.0; matG[idx_i, n1] += 1.0
            if n2 >= 0: matG[n2, idx_i] -= 1.0; matG[idx_i, n2] -= 1.0

        for r in self.resistors: stamp_admitance(G, r.n1.idx, r.n2.idx, 1.0 / r.R)
        for cap in self.capacitors: stamp_admitance(C, cap.n1.idx, cap.n2.idx, cap.C)
        for ind in self.inductors:
            stamp_vsrc(G, ind.n1.idx, ind.n2.idx, ind.idx_i)
            C[ind.idx_i, ind.idx_i] = -ind.L
        for i, sw in enumerate(self.switches):
            stamp_admitance(G, sw.n1.idx, sw.n2.idx, 1.0 / (sw.Ron if estados_sw[i] else sw.Roff))
        for i, d in enumerate(self.diodos):
            stamp_admitance(G, d.n1.idx, d.n2.idx, 1.0 / (d.Ron if estados_sw[len(self.switches)+i] else d.Roff))
        for vs in self.vsources:
            stamp_vsrc(G, vs.n1.idx, vs.n2.idx, vs.idx_i)
            W[vs.idx_i, vs.idx_u] = 1.0
        for isrc in self.isources:
            if isrc.n1.idx >= 0: W[isrc.n1.idx, isrc.idx_u] = -1.0
            if isrc.n2.idx >= 0: W[isrc.n2.idx, isrc.idx_u] = 1.0
        for vcv in self.vcvss:
            stamp_vsrc(G, vcv.n1.idx, vcv.n2.idx, vcv.idx_i)
            W[vcv.idx_i, vcv.idx_u] = vcv.gain
        for vcc in self.vccss:
            if vcc.n1.idx >= 0: W[vcc.n1.idx, vcc.idx_u] -= vcc.gm
            if vcc.n2.idx >= 0: W[vcc.n2.idx, vcc.idx_u] += vcc.gm
        for mut in self.mutuals:
            stamp_vsrc(G, mut.n1.idx, mut.n2.idx, mut.idx_i1)
            stamp_vsrc(G, mut.n3.idx, mut.n4.idx, mut.idx_i2)
            C[mut.idx_i1, mut.idx_i1] = -mut.L1; C[mut.idx_i1, mut.idx_i2] = -mut.M
            C[mut.idx_i2, mut.idx_i1] = -mut.M;  C[mut.idx_i2, mut.idx_i2] = -mut.L2

        A = G + C / self.dt
        try:
            A_inv = np.linalg.inv(A)
        except np.linalg.LinAlgError:
            A += np.eye(self.n_x) * 1e-12
            A_inv = np.linalg.inv(A)
        return A_inv @ W, A_inv @ (C / self.dt)
