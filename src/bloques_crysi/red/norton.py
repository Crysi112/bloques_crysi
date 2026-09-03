from __future__ import annotations
from dataclasses import dataclass
import warnings
import numpy as np
from typing import List, TYPE_CHECKING
if TYPE_CHECKING:
    from ..modelo import Modelo
    from .backend_base import BackendRed

@dataclass
class PuertoNorton:
    bus_pcc: int
    idx_v_pcc: int
    idx_ia: int
    idx_ib: int
    idx_ic: int
    v_nominal_ln: float
    delta_v: float = 0.01

class CoSimNorton:
    def __init__(
        self,
        modelo: "Modelo",
        backend: "BackendRed",
        puertos: List[PuertoNorton],
        dt_red: float = 0.02,
        tol_v: float = 1e-3,
        max_iter: int = 5,
        alpha: float = 0.8,
    ):
        self.modelo = modelo
        self.backend = backend
        self.puertos = puertos
        self.dt_red = float(dt_red)
        self.tol_v = float(tol_v)
        self.max_iter = int(max_iter)
        self.alpha = float(alpha)
        self._iniciado = False
        self._v_pcc = {p.bus_pcc: np.array([p.v_nominal_ln, p.v_nominal_ln, p.v_nominal_ln])
                       for p in self.puertos}
        self.t_actual = 0.0

    def _avanzar_ventana(self, dt_red: float) -> None:
        n_pasos = max(1, round(dt_red / self.modelo.dt))
        for _ in range(n_pasos):
            self.modelo.paso()
        self.t_actual += dt_red

    def _get_sig_array(self):
        m_c = self.modelo._paso_ctx[0]
        return np.ctypeslib.as_array(m_c.sig, shape=(m_c.n_sig,))

    def _leer_corrientes(self) -> dict:
        sigs = self._get_sig_array()
        result = {}
        for p in self.puertos:
            ia = sigs[p.idx_ia]
            ib = sigs[p.idx_ib]
            ic = sigs[p.idx_ic]
            result[p.bus_pcc] = np.array([ia, ib, ic])
        return result

    def run(self, t_fin: float, callback=None) -> dict:
        if not self._iniciado:
            self.modelo.iniciar(registrar=[self.modelo.bloques[0].nombre])
            self._iniciado = True
            self.t_actual = 0.0

        hist = {p.bus_pcc: {"t": [], "V_pcc": [], "P": [], "Q": [], "I": []}
                for p in self.puertos}
        fallos_consecutivos = 0

        while self.t_actual < t_fin - self.dt_red * 0.5:
            self._avanzar_ventana(self.dt_red)

            I_pccs = self._leer_corrientes()

            for p in self.puertos:
                I_abc = I_pccs[p.bus_pcc]
                V_abc = self._v_pcc[p.bus_pcc]

                P = np.sum(V_abc * np.real(I_abc))
                Q = np.sum(V_abc * np.imag(I_abc)) if np.iscomplexobj(I_abc) else 0.0

                try:
                    if P >= 0:
                        self.backend.set_carga(p.bus_pcc, P, Q)
                    else:
                        self.backend.set_generacion(p.bus_pcc, -P, -Q)
                except Exception as e:
                    warnings.warn(f"Backend error al actualizar bus {p.bus_pcc}: {e}")

            try:
                self.backend.runpp()
                fallos_consecutivos = 0

                for p in self.puertos:
                    try:
                        tb = self.backend.get_tension(p.bus_pcc)
                        V_nuevo = np.array([
                            tb.magnitud_v * np.sqrt(2),
                            tb.magnitud_v * np.sqrt(2),
                            tb.magnitud_v * np.sqrt(2),
                        ])
                        V_ant = self._v_pcc[p.bus_pcc]
                        self._v_pcc[p.bus_pcc] = self.alpha * V_nuevo + (1 - self.alpha) * V_ant

                        sigs = self._get_sig_array()
                        sigs[p.idx_v_pcc] = self._v_pcc[p.bus_pcc][0]

                    except Exception as e:
                        warnings.warn(f"No se pudo leer tension de bus {p.bus_pcc}: {e}")
            except Exception as e:
                fallos_consecutivos += 1
                if fallos_consecutivos >= self.max_iter:
                    raise RuntimeError(
                        f"CoSimNorton abortado: {fallos_consecutivos} divergencias consecutivas en el flujo de potencia."
                    ) from e
                warnings.warn(f"Error resolviendo flujo de potencia ({fallos_consecutivos}/{self.max_iter}): {e}")

            for p in self.puertos:
                I_abc = I_pccs[p.bus_pcc]
                V_abc = self._v_pcc[p.bus_pcc]
                P = float(np.sum(np.real(V_abc * np.conj(I_abc.astype(complex)))) if np.iscomplexobj(I_abc) else np.dot(V_abc, I_abc))
                hist[p.bus_pcc]["t"].append(self.t_actual)
                hist[p.bus_pcc]["V_pcc"].append(np.mean(np.abs(V_abc)))
                hist[p.bus_pcc]["P"].append(P)
                hist[p.bus_pcc]["Q"].append(0.0)
                hist[p.bus_pcc]["I"].append(np.mean(np.abs(I_abc)))

            if callback:
                callback(self.t_actual, self._v_pcc.copy(), I_pccs)

        for bus in hist:
            for k in hist[bus]:
                hist[bus][k] = np.array(hist[bus][k])

        return hist
