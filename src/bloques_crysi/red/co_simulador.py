from dataclasses import dataclass
from typing import Optional, Sequence, List, Tuple
import numpy as np
import warnings
from ..modelo import Modelo
from ..bloques import FuenteTrifasica, MedidorPotencia
from .backend_base import BackendRed, TensionBus, ErrorFlujoPotencia
@dataclass
class PCCConfig:
    bus_pcc: int
    medidor: MedidorPotencia
    fuente_red: FuenteTrifasica
    v_nominal_ln: Optional[float] = None
    es_generacion: bool = False
    fases: str = "ABC"
@dataclass
class ResultadoPCC:
    bus_pcc: int
    t: np.ndarray
    v_pcc_mag: np.ndarray
    v_pcc_ang: np.ndarray
    p_pcc: np.ndarray
    q_pcc: np.ndarray
    convergido: np.ndarray
    iteraciones_ventana: np.ndarray
@dataclass
class ResultadoCoSim:
    t: np.ndarray
    pccs: List[ResultadoPCC]
    def get_pcc(self, bus_pcc: int) -> Optional[ResultadoPCC]:
        for pcc in self.pccs:
            if pcc.bus_pcc == bus_pcc:
                return pcc
        return None
    def plot(self, bus_indices: Optional[Sequence[int]] = None, mostrar: bool = True):
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            warnings.warn("matplotlib no instalado; no se puede graficar.")
            return
        pccs_a_graficar = [p for p in self.pccs if bus_indices is None or p.bus_pcc in bus_indices]
        n = len(pccs_a_graficar)
        if n == 0:
            return
        fig, ejes = plt.subplots(n, 3, sharex=True, figsize=(14, 3.5 * n),
                                 layout="constrained")
        if n == 1:
            ejes = ejes.reshape(1, -1)
        colores = ["blue", "red", "limegreen", "magenta", "darkorange", "purple"]
        for i, pcc in enumerate(pccs_a_graficar):
            t = pcc.t
            ax_v = ejes[i, 0]
            ax_p = ejes[i, 1]
            ax_q = ejes[i, 2]
            ax_v.plot(t, pcc.v_pcc_mag, color=colores[i % len(colores)])
            ax_v.set_ylabel(f"V Bus {pcc.bus_pcc} [V]")
            ax_v.grid(True, linewidth=0.5, alpha=0.5)
            ax_v.set_xlim(t.min(), t.max())
            ax_p.plot(t, pcc.p_pcc / 1000, color=colores[i % len(colores)])
            ax_p.set_ylabel(f"P Bus {pcc.bus_pcc} [kW]")
            ax_p.grid(True, linewidth=0.5, alpha=0.5)
            ax_p.set_xlim(t.min(), t.max())
            ax_q.plot(t, pcc.q_pcc / 1000, color=colores[i % len(colores)])
            ax_q.set_ylabel(f"Q Bus {pcc.bus_pcc} [kVAr]")
            ax_q.grid(True, linewidth=0.5, alpha=0.5)
            ax_q.set_xlim(t.min(), t.max())
            if i == n - 1:
                ax_v.set_xlabel("Tiempo [s]")
                ax_p.set_xlabel("Tiempo [s]")
                ax_q.set_xlabel("Tiempo [s]")
            else:
                ax_v.tick_params(labelbottom=False)
                ax_p.tick_params(labelbottom=False)
                ax_q.tick_params(labelbottom=False)
        if mostrar:
            plt.show()
        return fig
    def guardar_csv(self, ruta: str, separador: str = ",") -> None:
        columnas = ["t"]
        datos = [self.t]
        for pcc in self.pccs:
            bus = pcc.bus_pcc
            columnas.extend([
                f"v_pcc_{bus}", f"p_pcc_{bus}", f"q_pcc_{bus}",
                f"conv_{bus}", f"iter_{bus}"
            ])
            datos.append(pcc.v_pcc_mag)
            datos.append(pcc.p_pcc)
            datos.append(pcc.q_pcc)
            datos.append(pcc.convergido.astype(int))
            datos.append(pcc.iteraciones_ventana)
        data = np.column_stack(datos)
        header = separador.join(columnas)
        np.savetxt(ruta, data, delimiter=separador, fmt="%.17g", header=header, comments="")
class CoSimuladorRed:
    def __init__(
        self,
        modelo: Modelo,
        backend: BackendRed,
        pccs: Sequence[PCCConfig],
        dt_red: float = 0.1,
        tol_convergencia_v: float = 1e-3,
        max_iter_ventana: int = 20,
        relajacion: float = 0.5,
    ):
        self.modelo = modelo
        self.backend = backend
        self.pccs = list(pccs)
        self.dt_red = float(dt_red)
        self.tol_v = float(tol_convergencia_v)
        self.max_iter = int(max_iter_ventana)
        self.alpha = float(relajacion)
        if not (0 < self.alpha <= 1.0):
            raise ValueError("relajacion (alpha) debe estar en (0, 1]")
        if self.dt_red <= 0:
            raise ValueError("dt_red debe ser > 0")
        if self.max_iter < 1:
            raise ValueError("max_iter_ventana debe ser >= 1")
        if self.pccs and len(set(pcc.bus_pcc for pcc in self.pccs)) != len(self.pccs):
            raise ValueError("Los buses PCC deben ser únicos")
        self._iniciado = False
        self._v_pcc_ant: dict = {}
        for pcc in self.pccs:
            v_nom = (pcc.v_nominal_ln if pcc.v_nominal_ln is not None
                    else float(pcc.fuente_red.param[0]) / np.sqrt(3))
            self._v_pcc_ant[pcc.bus_pcc] = TensionBus(magnitud_v=v_nom, angulo_rad=0.0)
    def acoplar_pcc(
        self,
        maquina,
        bus_idx: int,
        v_nominal_ll: float = 400.0,
        es_generacion: bool = False,
        medidor_existente: Optional[MedidorPotencia] = None,
    ) -> PCCConfig:
        if self._iniciado:
            raise RuntimeError(
                "No se pueden agregar PCCs: el modelo ya fue iniciado "
                "(run() ya se llamó al menos una vez)."
            )
        fuente = self.modelo.add(FuenteTrifasica(
            f"red_pcc{bus_idx}", amplitud=v_nominal_ll, frecuencia=50.0
        ))
        if medidor_existente is None:
            medidor = self.modelo.add(MedidorPotencia(f"pcc{bus_idx}", fases=3))
        else:
            medidor = medidor_existente
        self.modelo.conectar(fuente.salida, medidor.entrada)
        pcc = PCCConfig(
            bus_pcc=bus_idx,
            medidor=medidor,
            fuente_red=fuente,
            v_nominal_ln=v_nominal_ll / np.sqrt(3),
            es_generacion=es_generacion,
        )
        self.pccs.append(pcc)
        self._v_pcc_ant[bus_idx] = TensionBus(magnitud_v=pcc.v_nominal_ln, angulo_rad=0.0)
        return pcc
    def _leer_pq(self, valores: dict, medidor: MedidorPotencia) -> Tuple[float, float]:
        arr = np.atleast_1d(valores[medidor.nombre]).astype(float)
        p = float(arr[0])
        q = float(arr[1]) if medidor.fases == 3 and medidor.n_out == 3 else 0.0
        return p, q
    def _iterar_punto_fijo_ventana(
        self, n_pasos: int
    ) -> Tuple[dict, dict, dict, bool, int]:
        foto_inicio = self.modelo.guardar_estado()
        v_pcc = dict(self._v_pcc_ant)
        p_prom, q_prom = {}, {}
        convergido = False
        max_dv = float("inf")
        n_iter = 0
        for n_iter in range(1, self.max_iter + 1):
            self.modelo.restaurar_estado(foto_inicio)
            for pcc in self.pccs:
                v = v_pcc[pcc.bus_pcc]
                if not np.isfinite(v.magnitud_v) or v.magnitud_v <= 0:
                    v = TensionBus(magnitud_v=pcc.v_nominal_ln, angulo_rad=0.0, es_linea_linea=False)
                self.modelo.set_param(pcc.fuente_red, 0, v.magnitud_ll_v)
                if len(pcc.fuente_red.param) > 2:
                    self.modelo.set_param(pcc.fuente_red, 2, v.angulo_rad)
            p_sum = {pcc.bus_pcc: 0.0 for pcc in self.pccs}
            q_sum = {pcc.bus_pcc: 0.0 for pcc in self.pccs}
            for _ in range(n_pasos):
                valores = self.modelo.paso()
                for pcc in self.pccs:
                    p, q = self._leer_pq(valores, pcc.medidor)
                    p_sum[pcc.bus_pcc] += p
                    q_sum[pcc.bus_pcc] += q
            for pcc in self.pccs:
                p_med = p_sum[pcc.bus_pcc] / n_pasos
                q_med = q_sum[pcc.bus_pcc] / n_pasos
                if pcc.es_generacion:
                    p_prom[pcc.bus_pcc] = -p_med
                    q_prom[pcc.bus_pcc] = -q_med
                else:
                    p_prom[pcc.bus_pcc] = p_med
                    q_prom[pcc.bus_pcc] = q_med
            for pcc in self.pccs:
                bus_idx = pcc.bus_pcc
                p = p_prom[bus_idx]
                q = q_prom[bus_idx]
                if pcc.fases == "ABC":
                    if p >= 0:
                        self.backend.set_carga(bus_idx, p, q)
                    else:
                        self.backend.set_generacion(bus_idx, -p, -q)
                else:
                    if p >= 0:
                        self.backend.set_carga_fase(bus_idx, pcc.fases, p, q)
                    else:
                        self.backend.set_generacion_fase(bus_idx, pcc.fases, -p, -q)
            try:
                self.backend.runpp()
            except ErrorFlujoPotencia as e:
                warnings.warn(f"Flujo de potencia no convergió: {e}")
                v_nueva = {}
                for pcc in self.pccs:
                    v_nueva[pcc.bus_pcc] = TensionBus(
                        magnitud_v=pcc.v_nominal_ln,
                        angulo_rad=0.0,
                        es_linea_linea=False
                    )
            else:
                v_nueva = {pcc.bus_pcc: self.backend.get_tension(pcc.bus_pcc)
                          for pcc in self.pccs}
            for bus_idx, v in v_nueva.items():
                if not np.isfinite(v.magnitud_v) or v.magnitud_v <= 0:
                    v_nueva[bus_idx] = TensionBus(
                        magnitud_v=self._v_pcc_ant[bus_idx].magnitud_v,
                        angulo_rad=self._v_pcc_ant[bus_idx].angulo_rad,
                        es_linea_linea=False
                    )
            max_dv = 0.0
            v_relajada = {}
            alpha_eff = self.alpha
            if max_dv > self.tol_v * 10:
                alpha_eff = min(self.alpha, 0.3)
            for bus_idx, v_ant in v_pcc.items():
                vn = v_nueva[bus_idx]
                v_min = v_ant.magnitud_v * 0.5
                v_max = v_ant.magnitud_v * 1.5
                v_clamped = max(v_min, min(v_max, vn.magnitud_v))
                v_relajada[bus_idx] = TensionBus(
                    magnitud_v=self.alpha * v_clamped + (1 - self.alpha) * v_ant.magnitud_v,
                    angulo_rad=self.alpha * vn.angulo_rad + (1 - self.alpha) * v_ant.angulo_rad,
                    es_linea_linea=False,
                )
                max_dv = max(max_dv, abs(v_relajada[bus_idx].magnitud_v - v_ant.magnitud_v))
            v_pcc = v_relajada
            if max_dv < self.tol_v:
                convergido = True
                break
        self.modelo.restaurar_estado(foto_inicio)
        for pcc in self.pccs:
            v = v_pcc[pcc.bus_pcc]
            self.modelo.set_param(pcc.fuente_red, 0, v.magnitud_ll_v)
            if len(pcc.fuente_red.param) > 2:
                self.modelo.set_param(pcc.fuente_red, 2, v.angulo_rad)
        for _ in range(n_pasos):
            self.modelo.paso()
        if not convergido:
            warnings.warn(
                f"Punto fijo no convergió tras {n_iter} iteraciones "
                f"(max_dv={max_dv:.4g} V > tol={self.tol_v:g} V)."
            )
        return v_pcc, p_prom, q_prom, convergido, n_iter
    def run(self, t_fin: float) -> ResultadoCoSim:
        if t_fin <= 0:
            raise ValueError("t_fin debe ser > 0")
        if not self.pccs:
            raise ValueError("Se necesita al menos un PCCConfig (use acoplar_pcc o pase pccs al constructor)")
        n_pasos_ventana = max(1, round(self.dt_red / self.modelo.dt))
        dt_ventana_real = n_pasos_ventana * self.modelo.dt
        if abs(dt_ventana_real - self.dt_red) > 1e-9:
            warnings.warn(
                f"dt_red={self.dt_red:g}s no es múltiplo exacto de "
                f"modelo.dt={self.modelo.dt:g}s; se usa {dt_ventana_real:g}s "
                "como ventana real."
            )
        n_ventanas = max(1, round(t_fin / dt_ventana_real))
        if not self._iniciado:
            self.modelo.iniciar(registrar=[pcc.medidor for pcc in self.pccs])
            self._iniciado = True
        resultados_por_pcc = {
            pcc.bus_pcc: {'t': [], 'v_mag': [], 'v_ang': [],
                         'p': [], 'q': [], 'conv': [], 'iter': []}
            for pcc in self.pccs
        }
        t_actual = 0.0
        for _ in range(n_ventanas):
            t_actual += dt_ventana_real
            v_pcc, p_prom, q_prom, conv, n_iter = self._iterar_punto_fijo_ventana(
                n_pasos_ventana
            )
            self._v_pcc_ant = dict(v_pcc)
            for pcc in self.pccs:
                bus_idx = pcc.bus_pcc
                signo = -1 if pcc.es_generacion else 1
                r = resultados_por_pcc[bus_idx]
                r['t'].append(t_actual)
                r['v_mag'].append(v_pcc[bus_idx].magnitud_v)
                r['v_ang'].append(v_pcc[bus_idx].angulo_rad)
                r['p'].append(signo * p_prom[bus_idx])
                r['q'].append(signo * q_prom[bus_idx])
                r['conv'].append(conv)
                r['iter'].append(n_iter)
        pccs_resultado = []
        for pcc in self.pccs:
            bus_idx = pcc.bus_pcc
            r = resultados_por_pcc[bus_idx]
            pccs_resultado.append(ResultadoPCC(
                bus_pcc=bus_idx,
                t=np.array(r['t']),
                v_pcc_mag=np.array(r['v_mag']),
                v_pcc_ang=np.array(r['v_ang']),
                p_pcc=np.array(r['p']),
                q_pcc=np.array(r['q']),
                convergido=np.array(r['conv']),
                iteraciones_ventana=np.array(r['iter']),
            ))
        t_comun = pccs_resultado[0].t if pccs_resultado else np.array([])
        return ResultadoCoSim(t=t_comun, pccs=pccs_resultado)
    def acoplar_pcc(
        self,
        maquina,
        bus_idx: int,
        v_nominal_ll: float = 400.0,
        es_generacion: bool = False,
        medidor_existente: Optional[MedidorPotencia] = None,
    ) -> PCCConfig:
        if self._iniciado:
            raise RuntimeError(
                "No se pueden agregar PCCs: el modelo ya fue iniciado "
                "(run() ya se llamó al menos una vez)."
            )
        fuente = self.modelo.add(FuenteTrifasica(
            f"red_pcc{bus_idx}", amplitud=v_nominal_ll, frecuencia=50.0
        ))
        if medidor_existente is None:
            medidor = self.modelo.add(MedidorPotencia(f"pcc{bus_idx}", fases=3))
        else:
            medidor = medidor_existente
        self.modelo.conectar(fuente.salida, medidor.entrada)
        pcc = PCCConfig(
            bus_pcc=bus_idx,
            medidor=medidor,
            fuente_red=fuente,
            v_nominal_ln=v_nominal_ll / np.sqrt(3),
            es_generacion=es_generacion,
        )
        self.pccs.append(pcc)
        self._v_pcc_ant[bus_idx] = TensionBus(magnitud_v=pcc.v_nominal_ln, angulo_rad=0.0)
        return pcc
    def _leer_pq(self, valores: dict, medidor: MedidorPotencia) -> Tuple[float, float]:
        arr = np.atleast_1d(valores[medidor.nombre]).astype(float)
        p = float(arr[0])
        q = float(arr[1]) if medidor.fases == 3 and medidor.n_out == 3 else 0.0
        return p, q
    def reporte(self) -> None:
        if hasattr(self.backend, "reporte"):
            print(self.backend.reporte())
        else:
            print("El backend no cuenta con generador de reporte.")
def crear_cosimulador_simple(
    modelo: Modelo,
    backend: BackendRed,
    bus_pcc: int,
    medidor: MedidorPotencia,
    fuente_red: FuenteTrifasica,
    dt_red: float = 0.1,
    tol_convergencia_v: float = 1e-3,
    max_iter_ventana: int = 20,
    relajacion: float = 0.5,
    v_nominal_ln: Optional[float] = None,
) -> CoSimuladorRed:
    pcc = PCCConfig(
        bus_pcc=bus_pcc,
        medidor=medidor,
        fuente_red=fuente_red,
        v_nominal_ln=v_nominal_ln,
        es_generacion=False,
    )
    return CoSimuladorRed(
        modelo=modelo,
        backend=backend,
        pccs=[pcc],
        dt_red=dt_red,
        tol_convergencia_v=tol_convergencia_v,
        max_iter_ventana=max_iter_ventana,
        relajacion=relajacion,
    )
def comparar_backends(
    modelo: Modelo,
    backend_pp,
    backend_odss,
    pccs: Sequence[PCCConfig],
    dt_red: float = 0.1,
    t_fin: float = 1.0,
) -> dict:
    cosim_pp = CoSimuladorRed(modelo, backend_pp, pccs, dt_red)
    res_pp = cosim_pp.run(t_fin)
    cosim_odss = CoSimuladorRed(modelo, backend_odss, pccs, dt_red)
    res_odss = cosim_odss.run(t_fin)
    comparacion = {}
    for pcc_pp, pcc_odss in zip(res_pp.pccs, res_odss.pccs):
        bus = pcc_pp.bus_pcc
        comparacion[bus] = {
            'v_mag_diff': np.max(np.abs(pcc_pp.v_pcc_mag - pcc_odss.v_pcc_mag)),
            'v_ang_diff': np.max(np.abs(pcc_pp.v_pcc_ang - pcc_odss.v_pcc_ang)),
            'p_diff': np.max(np.abs(pcc_pp.p_pcc - pcc_odss.p_pcc)),
            'q_diff': np.max(np.abs(pcc_pp.q_pcc - pcc_odss.q_pcc)),
        }
    return comparacion
