from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Optional
@dataclass(frozen=True)
class TensionBus:
    magnitud_v: float
    angulo_rad: float
    es_linea_linea: bool = False
    @property
    def magnitud_ll_v(self) -> float:
        from .unidades import SQRT3
        return self.magnitud_v * SQRT3 if not self.es_linea_linea else self.magnitud_v
    @property
    def magnitud_ln_v(self) -> float:
        from .unidades import SQRT3
        return self.magnitud_v if not self.es_linea_linea else self.magnitud_v / SQRT3
@dataclass(frozen=True)
class PotenciaInyectada:
    p_w: float
    q_var: float
class ErrorFlujoPotencia(RuntimeError):
    def __init__(self, mensaje: str, backend: str, detalles: Optional[dict] = None):
        super().__init__(f"[{backend}] {mensaje}")
        self.backend = backend
        self.detalles = detalles or {}
class BackendRed(ABC):
    @abstractmethod
    def set_carga(self, bus_idx: int, p_w: float, q_var: float) -> None:
        ...
    @abstractmethod
    def set_generacion(self, bus_idx: int, p_w: float, q_var: float) -> None:
        ...
    @abstractmethod
    def runpp(self) -> None:
        ...
    @abstractmethod
    def get_tension(self, bus_idx: int) -> TensionBus:
        ...
    @abstractmethod
    def get_corriente_linea(self, linea_idx: int) -> float:
        ...
    @abstractmethod
    def desactivar_cargas_estaticas(self, bus_idx: int) -> None:
        ...
    @property
    @abstractmethod
    def nombre_backend(self) -> str:
        ...
    def __enter__(self) -> "BackendRed":
        return self
    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        pass
class BackendRedMock(BackendRed):
    def __init__(
        self,
        v_slack_v: float = 230.0,
        z_linea_ohm: complex = 0.1 + 0.05j,
    ):
        self._v_slack_v = v_slack_v
        self._z_linea = z_linea_ohm
        self._carga_bus1 = PotenciaInyectada(0.0, 0.0)
        self._generacion_bus1 = PotenciaInyectada(0.0, 0.0)
        self._convergido = False
        self._tension_bus1: Optional[TensionBus] = None
    def set_carga(self, bus_idx: int, p_w: float, q_var: float) -> None:
        if bus_idx != 1:
            raise ValueError(f"BackendMock solo soporta bus 1, recibido {bus_idx}")
        self._carga_bus1 = PotenciaInyectada(p_w, q_var)
    def set_generacion(self, bus_idx: int, p_w: float, q_var: float) -> None:
        if bus_idx != 1:
            raise ValueError(f"BackendMock solo soporta bus 1, recibido {bus_idx}")
        self._generacion_bus1 = PotenciaInyectada(p_w, q_var)
    def runpp(self) -> None:
        p_iny = self._generacion_bus1.p_w - self._carga_bus1.p_w
        q_iny = self._generacion_bus1.q_var - self._carga_bus1.q_var
        if abs(p_iny) < 1e-6 and abs(q_iny) < 1e-6:
            v1 = self._v_slack_v
            ang = 0.0
        else:
            v1 = self._v_slack_v
            for _ in range(5):
                p_01 = -p_iny
                q_01 = -q_iny
                i_01 = complex(p_01, -q_01) / v1
                dv = i_01 * self._z_linea
                v1_nuevo = self._v_slack_v - dv.real
                if abs(v1_nuevo - v1) < 1e-6:
                    v1 = v1_nuevo
                    break
                v1 = v1_nuevo
            ang = 0.0
        self._tension_bus1 = TensionBus(magnitud_v=v1, angulo_rad=ang, es_linea_linea=False)
        self._convergido = True
    def get_tension(self, bus_idx: int) -> TensionBus:
        if not self._convergido:
            raise RuntimeError("runpp() no se ha ejecutado o no convergió")
        if bus_idx == 0:
            return TensionBus(magnitud_v=self._v_slack_v, angulo_rad=0.0, es_linea_linea=False)
        if bus_idx == 1:
            if self._tension_bus1 is None:
                raise RuntimeError("runpp() no se ha ejecutado")
            return self._tension_bus1
        raise ValueError(f"Bus {bus_idx} no existe en BackendMock")
    def get_corriente_linea(self, linea_idx: int) -> float:
        if not self._convergido:
            raise RuntimeError("runpp() no se ha ejecutado o no convergió")
        if linea_idx != 0:
            raise ValueError("Solo existe línea 0 en BackendMock")
        if self._tension_bus1 is None:
            return 0.0
        v1 = self._tension_bus1.magnitud_v
        p_net = self._generacion_bus1.p_w - self._carga_bus1.p_w
        q_net = self._generacion_bus1.q_var - self._carga_bus1.q_var
        if v1 > 0:
            return (p_net**2 + q_net**2)**0.5 / v1
        return 0.0
    @property
    def nombre_backend(self) -> str:
        return "Mock"
    def desactivar_cargas_estaticas(self, bus_idx: int) -> None:
        if bus_idx == 1:
            self._carga_bus1 = PotenciaInyectada(0.0, 0.0)
            self._generacion_bus1 = PotenciaInyectada(0.0, 0.0)
            self._convergido = False
