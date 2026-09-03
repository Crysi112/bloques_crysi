"""
Tests para red.unidades - conversiones de unidades (fuente de verdad única).
"""

import math
import pytest
from bloques_crysi.red.unidades import (
    watts_a_mw, mw_a_watts, watts_a_kw, kw_a_watts,
    vars_a_mvar, mvar_a_vars, vars_a_kvar, kvar_a_vars,
    va_a_mva, mva_a_va, va_a_kva, kva_a_va,
    volts_a_kv, kv_a_volts,
    volts_linea_a_fase, fase_a_volts_linea,
    rad_a_grados, grados_a_rad,
    SistemaUnidades,
    pq_si_a_pandapower,
    tension_pandapower_a_si,
    tension_si_a_opendss,
    tension_opendss_a_si,
    pq_si_a_opendss,
    pq_opendss_a_si,
)


class TestPotencia:
    def test_watts_mw_ida_vuelta(self):
        assert mw_a_watts(watts_a_mw(1_500_000)) == 1_500_000
        assert abs(watts_a_mw(mw_a_watts(2.5)) - 2.5) < 1e-10

    def test_watts_kw_ida_vuelta(self):
        assert kw_a_watts(watts_a_kw(150_000)) == 150_000
        assert abs(watts_a_kw(kw_a_watts(150.0)) - 150.0) < 1e-10

    def test_vars_mvar_ida_vuelta(self):
        assert mvar_a_vars(vars_a_mvar(500_000)) == 500_000

    def test_va_mva_ida_vuelta(self):
        assert mva_a_va(va_a_mva(2_000_000)) == 2_000_000


class TestTension:
    def test_volts_kv_ida_vuelta(self):
        assert kv_a_volts(volts_a_kv(400_000)) == 400_000
        assert abs(volts_a_kv(kv_a_volts(13.8)) - 13.8) < 1e-10

    def test_linea_fase_ida_vuelta(self):
        v_ll = 400.0
        v_ln = volts_linea_a_fase(v_ll)
        assert abs(v_ln - 400.0 / math.sqrt(3)) < 1e-10
        assert abs(fase_a_volts_linea(v_ln) - v_ll) < 1e-10


class TestAngulo:
    def test_rad_grados_ida_vuelta(self):
        assert abs(grados_a_rad(rad_a_grados(math.pi)) - math.pi) < 1e-10
        assert abs(rad_a_grados(grados_a_rad(180.0)) - 180.0) < 1e-10

    def test_valores_conocidos(self):
        assert abs(rad_a_grados(math.pi) - 180.0) < 1e-10
        assert abs(rad_a_grados(math.pi / 2) - 90.0) < 1e-10
        assert abs(grados_a_rad(90.0) - math.pi / 2) < 1e-10


class TestSistemaUnidades:
    def test_bases_calculadas(self):
        s = SistemaUnidades(v_base_kv=13.8, s_base_mva=100.0)
        assert abs(s.v_base_ll_v - 13800.0) < 1e-6
        assert abs(s.v_base_ln_v - 13800.0 / math.sqrt(3)) < 1e-6
        assert abs(s.s_base_va - 100_000_000.0) < 1e-6
        assert abs(s.z_base_ohm - (13800.0 ** 2) / 100_000_000.0) < 1e-6
        assert abs(s.i_base_a - 100_000_000.0 / (math.sqrt(3) * 13800.0)) < 1e-6

    def test_pu_tension(self):
        s = SistemaUnidades(v_base_kv=13.8, s_base_mva=100.0)
        v_pu = s.v_a_pu(13800.0, es_linea_linea=True)
        assert abs(v_pu - 1.0) < 1e-10
        v_pu_ln = s.v_a_pu(13800.0 / math.sqrt(3), es_linea_linea=False)
        assert abs(v_pu_ln - 1.0) < 1e-10

    def test_pu_potencia(self):
        s = SistemaUnidades(v_base_kv=13.8, s_base_mva=100.0)
        assert abs(s.p_a_pu(50_000_000) - 0.5) < 1e-10
        assert abs(s.q_a_pu(30_000_000) - 0.3) < 1e-10

    def test_desde_pu(self):
        s = SistemaUnidades(v_base_kv=13.8, s_base_mva=100.0)
        assert abs(s.pu_a_v(1.0, es_linea_linea=True) - 13800.0) < 1e-6
        assert abs(s.pu_a_p(0.5) - 50_000_000) < 1e-6


class TestConvenienciaAdaptadores:
    def test_pq_si_a_pandapower(self):
        p_mw, q_mvar = pq_si_a_pandapower(100_000_000, 50_000_000, 100.0)
        assert abs(p_mw - 100.0) < 1e-10
        assert abs(q_mvar - 50.0) < 1e-10

    def test_tension_pandapower_a_si(self):
        v_v, ang_rad = tension_pandapower_a_si(1.0, 30.0, 13.8)
        # 1.0 pu de 13.8 kV LL = 13.8/sqrt(3) kV LN ≈ 7967 V LN
        assert abs(v_v - 13800.0 / math.sqrt(3)) < 1e-6
        assert abs(ang_rad - math.radians(30.0)) < 1e-10

    def test_tension_si_a_opendss(self):
        v_kv, ang_deg = tension_si_a_opendss(230.0, math.radians(45.0))
        assert abs(v_kv - 0.23) < 1e-10
        assert abs(ang_deg - 45.0) < 1e-10

    def test_pq_si_a_opendss(self):
        p_kw, q_kvar = pq_si_a_opendss(150_000, 80_000)
        assert abs(p_kw - 150.0) < 1e-10
        assert abs(q_kvar - 80.0) < 1e-10


if __name__ == "__main__":
    pytest.main([__file__, "-v"])