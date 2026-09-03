from math import pi
SQRT3 = 3 ** 0.5
def watts_a_mw(p_w: float) -> float:
    return p_w / 1_000_000.0
def mw_a_watts(p_mw: float) -> float:
    return p_mw * 1_000_000.0
def watts_a_kw(p_w: float) -> float:
    return p_w / 1_000.0
def kw_a_watts(p_kw: float) -> float:
    return p_kw * 1_000.0
def vars_a_mvar(q_var: float) -> float:
    return q_var / 1_000_000.0
def mvar_a_vars(q_mvar: float) -> float:
    return q_mvar * 1_000_000.0
def vars_a_kvar(q_var: float) -> float:
    return q_var / 1_000.0
def kvar_a_vars(q_kvar: float) -> float:
    return q_kvar * 1_000.0
def va_a_mva(s_va: float) -> float:
    return s_va / 1_000_000.0
def mva_a_va(s_mva: float) -> float:
    return s_mva * 1_000_000.0
def va_a_kva(s_va: float) -> float:
    return s_va / 1_000.0
def kva_a_va(s_kva: float) -> float:
    return s_kva * 1_000.0
def volts_a_kv(v_v: float) -> float:
    return v_v / 1_000.0
def kv_a_volts(v_kv: float) -> float:
    return v_kv * 1_000.0
def volts_linea_a_fase(v_ll: float) -> float:
    return v_ll / SQRT3
def fase_a_volts_linea(v_ln: float) -> float:
    return v_ln * SQRT3
def rad_a_grados(rad: float) -> float:
    return rad * 180.0 / pi
def grados_a_rad(deg: float) -> float:
    return deg * pi / 180.0
class SistemaUnidades:
    def __init__(
        self,
        v_base_kv: float,
        s_base_mva: float,
    ):
        self.v_base_kv = v_base_kv
        self.s_base_mva = s_base_mva
        self.v_base_ll_v = kv_a_volts(v_base_kv)
        self.v_base_ln_v = volts_linea_a_fase(self.v_base_ll_v)
        self.s_base_va = mva_a_va(s_base_mva)
        self.p_base_w = self.s_base_va
        self.q_base_var = self.s_base_va
        self.z_base_ohm = (self.v_base_ll_v ** 2) / self.s_base_va
        self.i_base_a = self.s_base_va / (SQRT3 * self.v_base_ll_v)
    def v_a_pu(self, v_v: float, es_linea_linea: bool = True) -> float:
        v_base = self.v_base_ll_v if es_linea_linea else self.v_base_ln_v
        return v_v / v_base
    def p_a_pu(self, p_w: float) -> float:
        return p_w / self.p_base_w
    def q_a_pu(self, q_var: float) -> float:
        return q_var / self.q_base_var
    def s_a_pu(self, s_va: float) -> float:
        return s_va / self.s_base_va
    def i_a_pu(self, i_a: float) -> float:
        return i_a / self.i_base_a
    def z_a_pu(self, z_ohm: float) -> float:
        return z_ohm / self.z_base_ohm
    def pu_a_v(self, v_pu: float, es_linea_linea: bool = True) -> float:
        v_base = self.v_base_ll_v if es_linea_linea else self.v_base_ln_v
        return v_pu * v_base
    def pu_a_p(self, p_pu: float) -> float:
        return p_pu * self.p_base_w
    def pu_a_q(self, q_pu: float) -> float:
        return q_pu * self.q_base_var
    def pu_a_s(self, s_pu: float) -> float:
        return s_pu * self.s_base_va
    def pu_a_i(self, i_pu: float) -> float:
        return i_pu * self.i_base_a
    def pu_a_z(self, z_pu: float) -> float:
        return z_pu * self.z_base_ohm
def pq_si_a_pandapower(p_w: float, q_var: float, s_base_mva: float) -> tuple[float, float]:
    return watts_a_mw(p_w), vars_a_mvar(q_var)
def tension_pandapower_a_si(vm_pu: float, va_deg: float, v_base_kv: float, es_linea_linea: bool = True) -> tuple[float, float]:
    s = SistemaUnidades(v_base_kv, 1.0)
    v_ll_v = s.pu_a_v(vm_pu, es_linea_linea=True)
    v_ln_v = volts_linea_a_fase(v_ll_v)
    ang_rad = grados_a_rad(va_deg)
    return v_ln_v, ang_rad
def tension_si_a_opendss(v_v: float, ang_rad: float) -> tuple[float, float]:
    v_kv = volts_a_kv(v_v)
    ang_deg = rad_a_grados(ang_rad)
    return v_kv, ang_deg
def tension_opendss_a_si(v_kv: float, ang_deg: float, es_linea_neutro: bool = True) -> tuple[float, float]:
    v_v = kv_a_volts(v_kv)
    if not es_linea_neutro:
        v_v = fase_a_volts_linea(v_v)
    ang_rad = grados_a_rad(ang_deg)
    return v_v, ang_rad
def pq_si_a_opendss(p_w: float, q_var: float) -> tuple[float, float]:
    return watts_a_kw(p_w), vars_a_kvar(q_var)
def pq_opendss_a_si(p_kw: float, q_kvar: float) -> tuple[float, float]:
    return kw_a_watts(p_kw), kvar_a_vars(q_kvar)
