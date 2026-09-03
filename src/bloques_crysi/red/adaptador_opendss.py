from typing import Optional
import cmath
import math
import py_dss_interface
from .backend_base import (
    BackendRed,
    TensionBus,
    ErrorFlujoPotencia,
)
from .unidades import (
    watts_a_kw, vars_a_kvar,
    kw_a_watts, kvar_a_vars,
    volts_linea_a_fase,
    kv_a_volts,
)
import numpy as np
class RedOpenDSS:
    def __init__(self, nombre: str = "Red", v_slack_kv_ll: float = 4.16, f_hz: float = 60.0,
                 fuente_ideal: bool = True):
        self.dss = py_dss_interface.DSS()
        self.dss.text("Clear")
        self.dss.text(f'New Circuit.{nombre} phases=3 basekv={v_slack_kv_ll} '
                      f'frequency={f_hz} pu=1.0 bus1=650')
        if fuente_ideal:
            self.dss.text("Vsource.Source.Isc3=1000000000")
            self.dss.text("Vsource.Source.Isc1=1000000000")
        self.v_base_kv_ll = v_slack_kv_ll
    def _formatear_matriz(self, matriz) -> str:
        if isinstance(matriz, (list, tuple)) and any(isinstance(r, (list, tuple)) and len(r) != len(matriz) for r in matriz):
            filas = [" ".join(f"{float(x):.6f}" for x in row) for row in matriz]
            return "[" + " | ".join(filas) + "]"
        m = np.asarray(matriz)
        filas = []
        for i in range(m.shape[0]):
            filas.append(" ".join(f"{float(x):.6f}" for x in m[i, :i+1]))
        return "[" + " | ".join(filas) + "]"
    def definir_linea_matricial(self, nombre: str, fases: int, r_matrix, x_matrix, c_matrix=None):
        r_str = self._formatear_matriz(r_matrix)
        x_str = self._formatear_matriz(x_matrix)
        cmd = f"New LineCode.{nombre} nphases={fases} rmatrix={r_str} xmatrix={x_str}"
        if c_matrix is not None:
            c_str = self._formatear_matriz(c_matrix)
            cmd += f" cmatrix={c_str}"
        self.dss.text(cmd)
    def definir_linea_simetrica(self, nombre: str, r1: float, x1: float, r0: float, x0: float):
        self.dss.text(f"New LineCode.{nombre} nphases=3 r1={r1} x1={x1} r0={r0} x0={x0}")
    def agregar_linea(self, nombre: str, bus1: str, bus2: str, longitud_km: float = None,
                      linecode: str = None, fases_n: int = 3, conexiones: str = "", longitud_ft: float = None, **kwargs):
        if longitud_ft is not None and longitud_km is None:
            longitud_km = float(longitud_ft) * 0.0003048
        if linecode is None:
            linecode = kwargs.get("linecode", kwargs.get("LineCode"))
        self.dss.text(f"New Line.{nombre} bus1={bus1}{conexiones} bus2={bus2}{conexiones} "
                      f"phases={fases_n} linecode={linecode} length={longitud_km} units=km")
    def agregar_transformador(self, nombre: str, bus_hv: str = None, bus_lv: str = None,
                              kv_hv: float = None, kv_lv: float = None, kva: float = None,
                              conn_hv: str = "wye", conn_lv: str = "wye",
                              x_hl: float = 2.0, r_perc: float = 1.1, **kwargs):
        if kv_hv is None and "kvs" in kwargs:
            kvs = kwargs["kvs"]; kv_hv, kv_lv = float(kvs[0]), float(kvs[1])
        if kva is None and "kvas" in kwargs:
            kvas = kwargs["kvas"]; kva = float(kvas[0])
        if "conns" in kwargs:
            conns = kwargs["conns"]; conn_hv, conn_lv = str(conns[0]), str(conns[1])
        if "xhl" in kwargs: x_hl = float(kwargs["xhl"])
        if "r_loadloss" in kwargs: r_perc = float(kwargs["r_loadloss"])
        if bus_hv is None: bus_hv = kwargs.get("bus1", kwargs.get("bus"))
        if bus_lv is None: bus_lv = kwargs.get("bus2")
        self.dss.text(f"New Transformer.{nombre} phases=3 windings=2 "
                      f"buses=({bus_hv}, {bus_lv}) "
                      f"conns=({conn_hv}, {conn_lv}) "
                      f"kvs=({kv_hv}, {kv_lv}) "
                      f"kvas=({kva}, {kva}) "
                      f"Xhl={x_hl} %LoadLoss={r_perc}")
    def agregar_carga_pq(self, nombre: str, bus: str, kw: float, kvar: float, kv: float, fases: int = 3, conn: str = ""):
        self.dss.text(f"New Load.{nombre} bus1={bus}{conn} phases={fases} kv={kv} kw={kw} kvar={kvar} model=1")
    def agregar_regulador(self, nombre: str, linea_controlada: str = None, bus_medido: str = None,
                          v_set: float = 122.0, pt_ratio: float = 20.0, ct_rating: float = 700.0, band: float = 2.0,
                          r_set: float = 3.0, x_set: float = 9.0, transformer: str = None, **kwargs):
        tr = transformer or linea_controlada
        self.dss.text(f"New RegControl.{nombre} transformer={tr} "
                      f"winding=2 vreg={v_set} ptratio={pt_ratio} CTPrim={ct_rating} band={band} R={r_set} X={x_set}")
    def agregar_banco_reguladores(self, nombre: str, bus_hv: str, bus_lv: str,
                                   kva_fase: float, kv_ln: float = None, kv_ll: float = None,
                                   v_set: float = 122.0, pt_ratio: float = 20.0,
                                   ct_rating: float = 700.0, band: float = 2.0,
                                   r_set: float = 3.0, x_set: float = 9.0,
                                   xhl: float = 0.01, r_perc: float = 0.01,
                                   fases=(1, 2, 3)):
        if kv_ln is None:
            if kv_ll is None:
                raise ValueError("Debes pasar kv_ln o kv_ll.")
            kv_ln = float(kv_ll) / (3 ** 0.5)
        for f in fases:
            tname = f"{nombre}_{f}"
            self.dss.text(
                f"New Transformer.{tname} phases=1 windings=2 "
                f"buses=({bus_hv}.{f}, {bus_lv}.{f}) "
                f"conns=(wye, wye) "
                f"kvs=({kv_ln}, {kv_ln}) "
                f"kvas=({kva_fase}, {kva_fase}) "
                f"Xhl={xhl} %LoadLoss={r_perc}"
            )
            self.dss.text(
                f"New RegControl.{nombre}_{f} transformer={tname} "
                f"winding=2 vreg={v_set} ptratio={pt_ratio} "
                f"CTPrim={ct_rating} band={band} R={r_set} X={x_set}"
            )
    def compilar(self) -> "BackendOpenDSS":
        self.dss.text("CalcVoltageBases")
        self.dss.text("Solve")
        return BackendOpenDSS(self.dss, self.v_base_kv_ll)
    def definir_linecode(self, nombre, fases=3, r_mat=None, x_mat=None, b_mat_us=None, r_matrix=None, x_matrix=None, c_matrix=None, **kwargs):
        if r_mat is not None: r_matrix = r_mat
        if x_mat is not None: x_matrix = x_mat
        if r_matrix is not None:
            r_matrix = [[r/1.60934 for r in row] for row in r_matrix]
        if x_matrix is not None:
            x_matrix = [[x/1.60934 for x in row] for row in x_matrix]
        if b_mat_us is not None:
            import numpy as np
            try:
                c_matrix = [[b*1e-6 / (2*np.pi*60) * 1e9 / 1.60934 for b in row] for row in b_mat_us]
            except Exception:
                c_matrix = b_mat_us
        return self.definir_linea_matricial(nombre, fases, r_matrix, x_matrix, c_matrix)
    def agregar_switch(self, nombre, bus1, bus2, fases=3):
        self.dss.text(f"New Line.{nombre} bus1={bus1} bus2={bus2} phases={fases} r1=1e-4 x1=0 c1=0 length=0.01 units=km")
    def agregar_capacitor(self, nombre, bus, kvar, kv, fases=3):
        ph = fases
        self.dss.text(f"New Capacitor.{nombre} bus1={bus} phases={ph} kv={kv} kvar={kvar}")
    def agregar_carga(self, nombre, bus, kw, kvar, kv, fases=3, conn="", model=1):
        base = kv
        self.dss.text(f"New Load.{nombre} bus1={bus} phases={fases} kv={base} kw={kw} kvar={kvar} model={model}")
    def agregar_transformador_subestacion(self, bus_hv="SourceBus", bus_lv="650", kva=5000, **kwargs):
        self.dss.text("Set VoltageBases=[115, 4.16, 0.48]")
        self.agregar_transformador("Sub", bus_hv=bus_hv, bus_lv=bus_lv, kv_hv=115, kv_lv=4.16, kva=kva, x_hl=8.0, r_perc=1.0)
    def agregar_regulador_monofasico(self, nombre, bus_in, bus_out, kv_ln, kva, **kwargs):
        tr = f"{nombre}_TR"
        self.agregar_transformador(tr, bus_in, bus_out, kv_hv=kv_ln*1.732, kv_lv=kv_ln*1.732, kva=kva, x_hl=0.01, r_perc=0.01)
        self.agregar_regulador(nombre, tr, bus_out, v_set=122.0, pt_ratio=20.0, ct_rating=700.0)
class BackendOpenDSS(BackendRed):
    @classmethod
    def desde_archivo(cls, archivo_dss: str, v_base_kv_ll: float) -> "BackendOpenDSS":
        import py_dss_interface
        dss = py_dss_interface.DSS()
        dss.text('Clear')
        dss.text(f'Redirect "{archivo_dss}"')
        dss.text('Calcv')
        return cls(dss, v_base_kv_ll)
    @classmethod
    def radial(cls, vn_kv: float = 0.4, r_ohm: float = 0.1, x_ohm: float = 0.05,
               longitud_km: float = 0.1) -> "BackendOpenDSS":
        import py_dss_interface
        dss = py_dss_interface.DSS()
        dss.text('Clear')
        dss.text(f'New Circuit.Test phases=3 basekv={vn_kv}')
        dss.text(f'New Source.Source Bus1=Source phases=3 basekv={vn_kv} pu=1.0')
        r_total = r_ohm * longitud_km
        x_total = x_ohm * longitud_km
        dss.text(f'New Line.Line Bus1=Source Bus2=PCC Phases=3 '
                 f'R1={r_total} X1={x_total} C1=0.0 Length={longitud_km} Units=km')
        dss.text('New Load.Load Bus1=PCC phases=3 kv=0.4 kw=0 kvar=0 model=1')
        dss.text('Calcv')
        return cls(dss, v_base_kv_ll=vn_kv)
    def __init__(
        self,
        dss: py_dss_interface.DSS,
        v_base_kv_ll: float,
    ):
        self._dss = dss
        self._convergido = False
        self._v_base_kv_ll = v_base_kv_ll
        self._bus_names = list(self._dss.circuit.buses_names)
        self._bus_to_idx = {name.lower(): i for i, name in enumerate(self._bus_names)}
        self._load_names_by_bus: dict[int, list[str]] = {}
        self._gen_names_by_bus: dict[int, list[str]] = {}
        self._refrescar_cache_elementos()
    def _refrescar_cache_elementos(self) -> None:
        self._load_names_by_bus.clear()
        self._gen_names_by_bus.clear()
        self._dss.loads.first()
        for _ in range(self._dss.loads.count):
            name = self._dss.loads.name
            buses = self._dss.cktelement.bus_names
            if buses:
                bus_name = buses[0].split('.')[0]
                bus_idx = self._bus_to_idx.get(bus_name.lower())
                if bus_idx is not None:
                    self._load_names_by_bus.setdefault(bus_idx, []).append(name)
            self._dss.loads.next()
        self._dss.generators.first()
        for _ in range(self._dss.generators.count):
            name = self._dss.generators.name
            buses = self._dss.cktelement.bus_names
            if buses:
                bus_name = buses[0].split('.')[0]
                bus_idx = self._bus_to_idx.get(bus_name.lower())
                if bus_idx is not None:
                    self._gen_names_by_bus.setdefault(bus_idx, []).append(name)
            self._dss.generators.next()
        self._dss.pvsystems.first()
        for _ in range(self._dss.pvsystems.count):
            name = self._dss.pvsystems.name
            buses = self._dss.cktelement.bus_names
            if buses:
                bus_name = buses[0].split('.')[0]
                bus_idx = self._bus_to_idx.get(bus_name.lower())
                if bus_idx is not None:
                    self._gen_names_by_bus.setdefault(bus_idx, []).append(name)
            self._dss.pvsystems.next()
    def _get_bus_name(self, bus_idx: int) -> str:
        if 0 <= bus_idx < len(self._bus_names):
            return self._bus_names[bus_idx]
        raise ValueError(f"Bus index {bus_idx} fuera de rango (0-{len(self._bus_names)-1})")
    def _activar_bus(self, bus_idx: int) -> None:
        bus_name = self._get_bus_name(bus_idx)
        self._dss.circuit.set_active_bus(bus_name)
    @property
    def nombre_backend(self) -> str:
        return "OpenDSS"
    @property
    def dss(self) -> py_dss_interface.DSS:
        return self._dss
    @property
    def v_base_kv_ll(self) -> float:
        return self._v_base_kv_ll
    def set_carga(self, bus_idx: int, p_w: float, q_var: float) -> None:
        bus_name = self._get_bus_name(bus_idx)
        self._dss.circuit.set_active_bus(bus_name)
        try:
            kV_ll = float(self._dss.bus.kv_base) * 1.7320508
        except Exception:
            kV_ll = self._v_base_kv_ll
        load_name = f"Load_Cosim_Bus{bus_idx}"
        if load_name not in self._load_names_by_bus.get(bus_idx, []):
            cmd = f"New Load.{load_name} bus1={bus_name} phases=3 kV={kV_ll} kW={watts_a_kw(p_w)} kvar={vars_a_kvar(q_var)} model=1"
            self._dss.text(cmd)
            self._load_names_by_bus.setdefault(bus_idx, []).append(load_name)
        else:
            self._dss.loads.name = load_name
            p_kw = watts_a_kw(p_w)
            q_kvar = vars_a_kvar(q_var)
            self._dss.loads.kw = p_kw
            self._dss.loads.kvar = q_kvar
        self._convergido = False
    def set_generacion(self, bus_idx: int, p_w: float, q_var: float) -> None:
        bus_name = self._get_bus_name(bus_idx)
        self._dss.circuit.set_active_bus(bus_name)
        try:
            kV_ll = float(self._dss.bus.kv_base) * 1.7320508
        except Exception:
            kV_ll = self._v_base_kv_ll
        gen_name = f"Gen_Cosim_Bus{bus_idx}"
        if gen_name not in self._gen_names_by_bus.get(bus_idx, []):
            cmd = f"New Generator.{gen_name} bus1={bus_name} phases=3 kV={kV_ll} kW={watts_a_kw(p_w)} kvar={vars_a_kvar(q_var)} model=1"
            self._dss.text(cmd)
            self._gen_names_by_bus.setdefault(bus_idx, []).append(gen_name)
        else:
            self._dss.generators.name = gen_name
            p_kw = watts_a_kw(p_w)
            q_kvar = vars_a_kvar(q_var)
            self._dss.generators.kw = p_kw
            self._dss.generators.kvar = q_kvar
        self._convergido = False
    def runpp(self) -> None:
        try:
            self._dss.solution.mode = 0
            self._dss.solution.solve()
        except Exception as e:
            self._convergido = False
            raise ErrorFlujoPotencia(
                mensaje=f"Error en OpenDSS Solve: {e}",
                backend="OpenDSS",
                detalles={"error": str(e)},
            ) from e
        if not self._dss.solution.converged:
            self._convergido = False
            raise ErrorFlujoPotencia(
                mensaje="Flujo de potencia OpenDSS no convergió",
                backend="OpenDSS",
                detalles={
                    "iterations": self._dss.solution.iterations,
                    "max_iterations": self._dss.solution.max_iterations,
                },
            )
        self._convergido = True
    def get_tension(self, bus_idx: int) -> TensionBus:
        if not self._convergido:
            raise RuntimeError("runpp() no se ha ejecutado o no convergió")
        bus_name = self._get_bus_name(bus_idx)
        self._dss.circuit.set_active_bus(bus_name)
        v_complex_arr = self._dss.bus.voltages
        if len(v_complex_arr) < 2:
            raise RuntimeError(f"No hay datos de tensión para bus {bus_name}")
        v_a_real = v_complex_arr[0]
        v_a_imag = v_complex_arr[1]
        v_a = complex(v_a_real, v_a_imag)
        mag_v = abs(v_a)
        ang_rad = cmath.phase(v_a)
        return TensionBus(magnitud_v=mag_v, angulo_rad=ang_rad, es_linea_linea=False)
    def get_corriente_linea(self, linea_idx: int) -> float:
        if not self._convergido:
            raise RuntimeError("runpp() no se ha ejecutado o no convergió")
        self._dss.lines.first()
        for i in range(self._dss.lines.count):
            if i == linea_idx:
                self._dss.circuit.set_active_element(f"Line.{self._dss.lines.name}")
                currents = self._dss.cktelement.currents
                if len(currents) >= 2:
                    i_a = abs(complex(currents[0], currents[1]))
                    i_b = abs(complex(currents[2], currents[3])) if len(currents) >= 4 else 0.0
                    i_c = abs(complex(currents[4], currents[5])) if len(currents) >= 6 else 0.0
                    n_phases = len(currents) // 2
                    return (i_a + i_b + i_c) / n_phases
                return 0.0
            self._dss.lines.next()
        raise ValueError(f"Línea índice {linea_idx} no encontrada")
    def get_potencias_bus(self, bus_idx: int) -> tuple[float, float]:
        if not self._convergido:
            raise RuntimeError("runpp() no se ha ejecutado o no convergió")
        p_total, q_total = 0.0, 0.0
        for load_name in self._load_names_by_bus.get(bus_idx, []):
            self._dss.loads.name = load_name
            p_total -= kw_a_watts(self._dss.loads.kw)
            q_total -= kvar_a_vars(self._dss.loads.kvar)
        for gen_name in self._gen_names_by_bus.get(bus_idx, []):
            try:
                self._dss.generators.name = gen_name
                p_total += kw_a_watts(self._dss.generators.kw)
                q_total += kvar_a_vars(self._dss.generators.kvar)
            except Exception:
                self._dss.circuit.set_active_element(f"PVSystem.{gen_name}")
                pwrs = self._dss.cktelement.powers
                if pwrs:
                    p_total -= kw_a_watts(sum(pwrs[0::2]))
                    q_total -= kvar_a_vars(sum(pwrs[1::2]))
        return p_total, q_total
    def set_carga_fase(self, bus_idx: int, fase: str, p_w: float, q_var: float) -> None:
        bus_name = self._get_bus_name(bus_idx)
        self._dss.circuit.set_active_bus(bus_name)
        try:
            kV_ll = float(self._dss.bus.kv_base) * 1.7320508
        except Exception:
            kV_ll = self._v_base_kv_ll
        fase_idx = self._fase_a_indice(fase)
        load_name = f"Load_Cosim_Bus{bus_idx}_{fase}"
        kV_ln = kV_ll / 3**0.5
        p_kw = watts_a_kw(p_w)
        q_kvar = vars_a_kvar(q_var)
        bus_fase = f"{bus_name}.{fase_idx + 1}"
        cmd = f"New Load.{load_name} bus1={bus_fase} phases=1 kV={kV_ln} kW={p_kw} kvar={q_kvar} model=1"
        self._dss.text(cmd)
        self._convergido = False
    def set_generacion_fase(self, bus_idx: int, fase: str, p_w: float, q_var: float) -> None:
        bus_name = self._get_bus_name(bus_idx)
        self._dss.circuit.set_active_bus(bus_name)
        try:
            kV_ll = float(self._dss.bus.kv_base) * 1.7320508
        except Exception:
            kV_ll = self._v_base_kv_ll
        fase_idx = self._fase_a_indice(fase)
        gen_name = f"Gen_Cosim_Bus{bus_idx}_{fase}"
        kV_ln = kV_ll / 3**0.5
        p_kw = watts_a_kw(p_w)
        q_kvar = vars_a_kvar(q_var)
        bus_fase = f"{bus_name}.{fase_idx + 1}"
        cmd = f"New Generator.{gen_name} bus1={bus_fase} phases=1 kV={kV_ln} kW={p_kw} kvar={q_kvar} model=1"
        self._dss.text(cmd)
        self._convergido = False
    def _fase_a_indice(self, fase) -> int:
        if isinstance(fase, int):
            if 0 <= fase <= 2:
                return fase
            raise ValueError("Fase entera debe ser 0, 1 o 2")
        fase_str = str(fase).upper()
        if fase_str in ('A', 'B', 'C'):
            return ord(fase_str) - ord('A')
        raise ValueError(f"Fase inválida: {fase}. Use 'A', 'B', 'C' o 0, 1, 2")
    def get_tensiones_fase(self, bus_idx: int) -> list[TensionBus]:
        if not self._convergido:
            raise RuntimeError("runpp() no se ha ejecutado o no convergió")
        bus_name = self._get_bus_name(bus_idx)
        self._dss.circuit.set_active_bus(bus_name)
        v_complex_arr = self._dss.bus.voltages
        nodes = self._dss.bus.nodes
        if len(v_complex_arr) < 2:
            raise RuntimeError(f"Datos de tensión incompletos para bus {bus_name}")
        n_fases = len(v_complex_arr) // 2
        if len(nodes) >= n_fases and n_fases <= 3:
            n_fases = min(len(nodes), n_fases)
        tensiones = []
        for i in range(n_fases):
            if i * 2 + 1 >= len(v_complex_arr):
                break
            v_real = v_complex_arr[i * 2]
            v_imag = v_complex_arr[i * 2 + 1]
            v = complex(v_real, v_imag)
            tensiones.append(TensionBus(magnitud_v=abs(v), angulo_rad=cmath.phase(v), es_linea_linea=False))
        while len(tensiones) < 3 and n_fases < 3:
            tensiones.append(TensionBus(magnitud_v=0.0, angulo_rad=0.0, es_linea_linea=False))
            n_fases += 1
        return tensiones[:3] if len(tensiones) == 3 else tensiones
    def get_corrientes_fase(self, bus_idx: int) -> list[complex]:
        raise NotImplementedError("OpenDSS no provee corrientes inyectadas directas a nivel de bus. Caclúlalas iterando sus elementos conectados.")
    def desactivar_cargas_estaticas(self, bus_idx: int) -> None:
        bus_name = self._get_bus_name(bus_idx)
        load_name = self._load_names_by_bus.get(bus_idx)
        if load_name:
            try:
                self._dss.loads.name = load_name
                self._dss.loads.enabled = False
            except:
                self._dss.text(f"Load.{load_name}.enabled=No")
        gen_name = self._gen_names_by_bus.get(bus_idx)
        if gen_name:
            try:
                self._dss.generators.name = gen_name
                self._dss.generators.enabled = False
            except:
                self._dss.text(f"Generator.{gen_name}.enabled=No")
        self._dss.pvsystems.first()
        for _ in range(self._dss.pvsystems.count):
            pv_name = self._dss.pvsystems.name
            buses = self._dss.cktelement.bus_names
            if buses and buses[0].split('.')[0].lower() == bus_name.lower():
                try:
                    self._dss.pvsystems.enabled = False
                except:
                    self._dss.text(f"PVSystem.{pv_name}.enabled=No")
            self._dss.pvsystems.next()
        self._convergido = False
    def reporte(self) -> str:
        if not self._convergido:
            return "El flujo de potencia no ha sido ejecutado o no convergió."
        lineas_out = []
        lineas_out.append("=" * 80)
        lineas_out.append(f" REPORTE DE FLUJO DE POTENCIA (Backend: {self.nombre_backend})")
        lineas_out.append("=" * 80)
        lineas_out.append("\n--- PERFIL DE TENSIONES POR NODO ---")
        lineas_out.append(f"{'BUS':<12} | {'FASE A (p.u. / deg)':<22} | {'FASE B (p.u. / deg)':<22} | {'FASE C (p.u. / deg)':<22} |")
        lineas_out.append("-" * 86)
        for bus_name in self._dss.circuit.buses_names:
            self._dss.circuit.set_active_bus(bus_name)
            nodes = self._dss.bus.nodes
            v_arr = self._dss.bus.voltages
            v_base_ln = self._dss.bus.kv_base * 1000.0
            fases_str = ["        ---         ", "        ---         ", "        ---         "]
            for idx, f in enumerate(nodes):
                if 1 <= f <= 3 and len(v_arr) >= 2 * (idx + 1):
                    v_c = complex(v_arr[2 * idx], v_arr[2 * idx + 1])
                    mag_pu = abs(v_c) / v_base_ln if v_base_ln > 0 else 0.0
                    ang_deg = math.degrees(cmath.phase(v_c))
                    fases_str[f - 1] = f"{mag_pu:6.4f} pu / {ang_deg:6.2f}°"
            lineas_out.append(f"{bus_name:<12} | {fases_str[0]:<22} | {fases_str[1]:<22} | {fases_str[2]:<22} |")
        lineas_out.append("\n--- FLUJOS DE POTENCIA Y PÉRDIDAS POR ELEMENTO ---")
        lineas_out.append(f"{'ELEMENTO':<20} | {'BUS 1':<10} | {'BUS 2':<10} | {'P (kW)':<12} | {'Q (kVAr)':<12} | {'PÉRDIDAS P (kW)':<15} |")
        lineas_out.append("-" * 92)
        try:
            self._dss.lines.first()
            for _ in range(self._dss.lines.count):
                nombre = self._dss.lines.name
                self._dss.circuit.set_active_element(f"Line.{nombre}")
                buses = self._dss.cktelement.bus_names
                b1 = buses[0].split('.')[0] if buses else "N/A"
                b2 = buses[1].split('.')[0] if len(buses) > 1 else "N/A"
                powers = self._dss.cktelement.powers
                losses = self._dss.cktelement.losses
                p_total = sum(powers[0::2]) if powers else 0.0
                q_total = sum(powers[1::2]) if powers else 0.0
                p_loss_kw = (losses[0] / 1000.0) if losses else 0.0
                lineas_out.append(f"Line.{nombre:<15} | {b1:<10} | {b2:<10} | {p_total:12.2f} | {q_total:12.2f} | {p_loss_kw:15.3f} |")
                self._dss.lines.next()
        except Exception:
            pass
        lineas_out.append("=" * 80)
        return "\n".join(lineas_out)
