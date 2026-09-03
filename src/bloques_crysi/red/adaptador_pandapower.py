from typing import Optional
import pandapower as pp
import pandapower.networks as pp_nets
from .backend_base import (
    BackendRed,
    TensionBus,
    ErrorFlujoPotencia,
)
from .unidades import (
    watts_a_mw, vars_a_mvar,
    tension_pandapower_a_si,
    SistemaUnidades,
)
class BackendPandapower(BackendRed):
    @classmethod
    def radial(cls, vn_kv: float = 0.4, r_ohm: float = 0.1, x_ohm: float = 0.05,
               longitud_km: float = 0.1, s_base_mva: float = 100.0, **runpp_kwargs):
        net = pp.create_empty_network()
        b0 = pp.create_bus(net, vn_kv=vn_kv, name="Slack")
        pp.create_ext_grid(net, bus=b0, vm_pu=1.0, va_degree=0.0)
        b1 = pp.create_bus(net, vn_kv=vn_kv, name="PCC")
        pp.create_line(net, from_bus=b0, to_bus=b1, length_km=longitud_km,
                       std_type="NAYY 4x50 SE", name="Linea_0_1")
        return cls(net, v_base_kv=vn_kv, s_base_mva=s_base_mva, **runpp_kwargs)
    @classmethod
    def desde_red(cls, nombre_red: str, **runpp_kwargs):
        net_func = getattr(pp_nets, f"create_{nombre_red}", None)
        if net_func is None:
            raise ValueError(f"Red estándar '{nombre_red}' no encontrada en pandapower.networks")
        net = net_func()
        if len(net.ext_grid) > 0:
            slack_bus = net.ext_grid.at[0, "bus"]
            v_base = float(net.bus.at[slack_bus, "vn_kv"])
        else:
            v_base = float(net.bus.iloc[0]["vn_kv"])
        return cls(net, v_base_kv=v_base, **runpp_kwargs)
    def __init__(
        self,
        net: pp.pandapowerNet,
        v_base_kv: Optional[float] = None,
        s_base_mva: float = 100.0,
        algoritmo: str = "nr",
        **runpp_kwargs,
    ):
        self._net = net
        self._algoritmo = algoritmo
        self._runpp_kwargs = runpp_kwargs
        self._convergido = False
        if v_base_kv is not None:
            self._v_base_kv = v_base_kv
        else:
            if len(net.ext_grid) > 0:
                slack_bus = net.ext_grid.at[0, "bus"]
                self._v_base_kv = float(net.bus.at[slack_bus, "vn_kv"])

    @property
    def nombre_backend(self) -> str:
        return "pandapower"
    @property
    def net(self) -> pp.pandapowerNet:
        return self._net
    @property
    def v_base_kv(self) -> float:
        return self._v_base_kv
    def _obtener_o_crear_carga(self, bus_idx: int) -> int:
        nombre = f"Load_Cosim_Bus{bus_idx}"
        mask = self._net.load["name"] == nombre
        if mask.any():
            return self._net.load[mask].index[0]
        return pp.create_load(self._net, bus=bus_idx, p_mw=0.0, q_mvar=0.0, name=nombre)
    def _obtener_o_crear_generador(self, bus_idx: int) -> int:
        nombre = f"Gen_Cosim_Bus{bus_idx}"
        mask = self._net.sgen["name"] == nombre
        if mask.any():
            return self._net.sgen[mask].index[0]
        return pp.create_sgen(self._net, bus=bus_idx, p_mw=0.0, q_mvar=0.0, name=nombre)
    def set_carga(self, bus_idx: int, p_w: float, q_var: float) -> None:
        if bus_idx not in self._net.bus.index:
            raise ValueError(f"Bus {bus_idx} no existe en la red pandapower")
        idx = self._obtener_o_crear_carga(bus_idx)
        p_mw = watts_a_mw(p_w)
        q_mvar = vars_a_mvar(q_var)
        self._net.load.at[idx, "p_mw"] = p_mw
        self._net.load.at[idx, "q_mvar"] = q_mvar
        self._net.load.at[idx, "in_service"] = True
        self._convergido = False
    def set_generacion(self, bus_idx: int, p_w: float, q_var: float) -> None:
        if bus_idx not in self._net.bus.index:
            raise ValueError(f"Bus {bus_idx} no existe en la red pandapower")
        idx = self._obtener_o_crear_generador(bus_idx)
        p_mw = watts_a_mw(p_w)
        q_mvar = vars_a_mvar(q_var)
        self._net.sgen.at[idx, "p_mw"] = p_mw
        self._net.sgen.at[idx, "q_mvar"] = q_mvar
        self._net.sgen.at[idx, "in_service"] = True
        self._convergido = False
    def runpp(self) -> None:
        try:
            has_asym = False
            try:
                if hasattr(self._net, "asymmetric_load") and not self._net.asymmetric_load.empty:
                    has_asym = True
                if hasattr(self._net, "asymmetric_sgen") and not self._net.asymmetric_sgen.empty:
                    has_asym = True
            except Exception:
                pass
            if has_asym:
                try:
                    pp.runpp_3ph(self._net, algorithm=self._algoritmo, **self._runpp_kwargs)
                except AttributeError:
                    pp.runpp(self._net, algorithm=self._algoritmo, **self._runpp_kwargs)
            else:
                pp.runpp(
                    self._net,
                    algorithm=self._algoritmo,
                    **self._runpp_kwargs,
                )
        except pp.LoadflowNotConverged as e:
            self._convergido = False
            raise ErrorFlujoPotencia(
                mensaje="Flujo de potencia no convergió",
                backend="pandapower",
                detalles={"algorithm": self._algoritmo, "error": str(e)},
            ) from e
        except Exception as e:
            self._convergido = False
            raise ErrorFlujoPotencia(
                mensaje=f"Error en runpp: {e}",
                backend="pandapower",
                detalles={"error": str(e)},
            ) from e
        self._convergido = True
    def get_tension(self, bus_idx: int) -> TensionBus:
        if not self._convergido:
            raise RuntimeError("runpp() no se ha ejecutado o no convergió")
        if bus_idx not in self._net.bus.index:
            raise ValueError(f"Bus {bus_idx} no existe en la red pandapower")
        if bus_idx not in self._net.res_bus.index:
            raise RuntimeError(f"No hay resultados para bus {bus_idx} (¿runpp falló?)")
        vm_pu = float(self._net.res_bus.at[bus_idx, "vm_pu"])
        va_deg = float(self._net.res_bus.at[bus_idx, "va_degree"])
        vn_kv_bus = float(self._net.bus.at[bus_idx, "vn_kv"])
        v_v, ang_rad = tension_pandapower_a_si(vm_pu, va_deg, vn_kv_bus, es_linea_linea=True)
        return TensionBus(magnitud_v=v_v, angulo_rad=ang_rad, es_linea_linea=False)
    def get_corriente_linea(self, linea_idx: int) -> float:
        if not self._convergido:
            raise RuntimeError("runpp() no se ha ejecutado o no convergió")
        if linea_idx not in self._net.line.index:
            raise ValueError(f"Línea {linea_idx} no existe en la red pandapower")
        if "res_line" not in dir(self._net) or linea_idx not in self._net.res_line.index:
            raise RuntimeError(f"No hay resultados de corriente para línea {linea_idx}")
        i_ka_a = float(self._net.res_line.at[linea_idx, "i_ka"])
        return i_ka_a * 1000.0
    def get_potencias_bus(self, bus_idx: int) -> tuple[float, float]:
        if not self._convergido:
            raise RuntimeError("runpp() no se ha ejecutado o no convergió")
        p_total, q_total = 0.0, 0.0
        loads = self._net.load[self._net.load["bus"] == bus_idx]
        if not loads.empty:
            p_total -= loads["p_mw"].sum() * 1e6
            q_total -= loads["q_mvar"].sum() * 1e6
        sgens = self._net.sgen[self._net.sgen["bus"] == bus_idx]
        if not sgens.empty:
            p_total += sgens["p_mw"].sum() * 1e6
            q_total += sgens["q_mvar"].sum() * 1e6
        gens = self._net.gen[self._net.gen["bus"] == bus_idx]
        if not gens.empty:
            p_total += gens["p_mw"].sum() * 1e6
            if "res_gen" in self._net and not self._net.res_gen.empty:
                q_total += self._net.res_gen.loc[gens.index, "q_mvar"].sum() * 1e6
        return p_total, q_total
    def desactivar_cargas_estaticas(self, bus_idx: int) -> None:
        if bus_idx not in self._net.bus.index:
            return
        if len(self._net.load) > 0:
            self._net.load.loc[self._net.load["bus"] == bus_idx, "in_service"] = False
        if len(self._net.sgen) > 0:
            self._net.sgen.loc[self._net.sgen["bus"] == bus_idx, "in_service"] = False
        if len(self._net.gen) > 0:
            self._net.gen.loc[self._net.gen["bus"] == bus_idx, "in_service"] = False
        self._convergido = False
    def set_carga_fase(self, bus_idx: int, fase: str, p_w: float, q_var: float) -> None:
        fase = str(fase).lower()
        if fase not in ['a', 'b', 'c']:
            raise ValueError("La fase debe ser 'A', 'B' o 'C'")
        p_mw = watts_a_mw(p_w)
        q_mvar = vars_a_mvar(q_var)
        if hasattr(self._net, "asymmetric_load") and not self._net.asymmetric_load.empty:
            mask = self._net.asymmetric_load["bus"] == bus_idx
            if mask.any():
                idx = self._net.asymmetric_load[mask].index[0]
                self._net.asymmetric_load.at[idx, f"p_{fase}_mw"] = p_mw
                self._net.asymmetric_load.at[idx, f"q_{fase}_mvar"] = q_mvar
                self._net.asymmetric_load.at[idx, "in_service"] = True
                self._convergido = False
                return
        import pandapower as pp
        idx = pp.create_asymmetric_load(self._net, bus=bus_idx,
                                        p_a_mw=0, p_b_mw=0, p_c_mw=0,
                                        q_a_mvar=0, q_b_mvar=0, q_c_mvar=0)
        self._net.asymmetric_load.at[idx, f"p_{fase}_mw"] = p_mw
        self._net.asymmetric_load.at[idx, f"q_{fase}_mvar"] = q_mvar
        self._convergido = False
    def set_generacion_fase(self, bus_idx: int, fase: str, p_w: float, q_var: float) -> None:
        self.set_carga_fase(bus_idx, fase, -p_w, -q_var)
    def reporte(self) -> str:
        if not self._convergido:
            return "El flujo de potencia no ha sido ejecutado o no convergió."
        lineas_out = []
        lineas_out.append("=" * 80)
        lineas_out.append(f" REPORTE DE FLUJO DE POTENCIA (Backend: {self.nombre_backend})")
        lineas_out.append("=" * 80)
        lineas_out.append("\n--- PERFIL DE TENSIONES POR NODO ---")
        lineas_out.append(f"{'BUS':<12} | {'VM (p.u.)':<12} | {'ANG (deg)':<12} |")
        lineas_out.append("-" * 40)
        for bus_idx in self._net.bus.index:
            if bus_idx not in self._net.res_bus.index:
                continue
            nombre = self._net.bus.at[bus_idx, 'name']
            vm = self._net.res_bus.at[bus_idx, 'vm_pu']
            va = self._net.res_bus.at[bus_idx, 'va_degree']
            lineas_out.append(f"{nombre:<12} | {vm:<12.4f} | {va:<12.2f} |")
        lineas_out.append("\n--- FLUJOS DE POTENCIA Y PÉRDIDAS POR LÍNEA ---")
        lineas_out.append(f"{'ELEMENTO':<20} | {'BUS 1':<10} | {'BUS 2':<10} | {'P (kW)':<12} | {'Q (kVAr)':<12} | {'PÉRDIDAS (kW)':<15} |")
        lineas_out.append("-" * 92)
        for l_idx, linea in self._net.line.iterrows():
            if l_idx not in self._net.res_line.index:
                continue
            b1 = int(linea.from_bus)
            b2 = int(linea.to_bus)
            n1 = self._net.bus.at[b1, 'name'] if b1 in self._net.bus.index else str(b1)
            n2 = self._net.bus.at[b2, 'name'] if b2 in self._net.bus.index else str(b2)
            p_kw = self._net.res_line.at[l_idx, 'p_from_mw'] * 1000
            q_kvar = self._net.res_line.at[l_idx, 'q_from_mvar'] * 1000
            p_loss_kw = self._net.res_line.at[l_idx, 'pl_mw'] * 1000
            lineas_out.append(f"Line.{l_idx:<15} | {n1:<10} | {n2:<10} | {p_kw:12.2f} | {q_kvar:12.2f} | {p_loss_kw:15.3f} |")
        lineas_out.append("=" * 80)
        return "\n".join(lineas_out)
