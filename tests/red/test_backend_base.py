"""
Tests para red.backend_base - interfaz base y mock.
"""

import math
import pytest
from bloques_crysi.red.backend_base import (
    BackendRed,
    BackendRedMock,
    TensionBus,
    PotenciaInyectada,
    ErrorFlujoPotencia,
)


class TestTensionBus:
    def test_creacion_fase_neutro(self):
        tb = TensionBus(magnitud_v=230.0, angulo_rad=0.0, es_linea_linea=False)
        assert tb.magnitud_v == 230.0
        assert tb.magnitud_ln_v == 230.0
        assert abs(tb.magnitud_ll_v - 230.0 * math.sqrt(3)) < 1e-10

    def test_creacion_linea_linea(self):
        tb = TensionBus(magnitud_v=400.0, angulo_rad=0.0, es_linea_linea=True)
        assert tb.magnitud_v == 400.0
        assert tb.magnitud_ll_v == 400.0
        assert abs(tb.magnitud_ln_v - 400.0 / math.sqrt(3)) < 1e-10

    def test_inmutable(self):
        tb = TensionBus(magnitud_v=230.0, angulo_rad=0.0)
        with pytest.raises(Exception):
            tb.magnitud_v = 240.0


class TestPotenciaInyectada:
    def test_creacion(self):
        p = PotenciaInyectada(p_w=100_000, q_var=50_000)
        assert p.p_w == 100_000
        assert p.q_var == 50_000


class TestErrorFlujoPotencia:
    def test_atributos(self):
        err = ErrorFlujoPotencia("No convergió", "pandapower", {"iter": 10})
        assert err.backend == "pandapower"
        assert err.detalles == {"iter": 10}
        assert "pandapower" in str(err)


class TestBackendRedMock:
    def test_setup_basico(self):
        mock = BackendRedMock(v_slack_v=230.0)
        assert mock.nombre_backend == "Mock"

    def test_set_carga_bus_invalido(self):
        mock = BackendRedMock()
        with pytest.raises(ValueError):
            mock.set_carga(99, 1000, 500)

    def test_set_carga_y_generacion(self):
        mock = BackendRedMock()
        mock.set_carga(1, 100_000, 50_000)      # 100 kW, 50 kVAr carga
        mock.set_generacion(1, 20_000, 10_000)  # 20 kW, 10 kVAr generación
        mock.runpp()
        v = mock.get_tension(1)
        # Red simple: V1 ≈ V0 - (P - jQ)/V0 * Z
        # P_net = 20k - 100k = -80kW, Q_net = 10k - 50k = -40kVAr
        # I ≈ (80k + j40k) / 230 ≈ 348 + j174 A
        # dV ≈ I * (0.1 + j0.05) ≈ 34.8 + j17.4 + j17.4 - 8.7 ≈ 26.1 + j34.8
        # V1 ≈ 230 - 26.1 ≈ 203.9 V (muy aproximado)
        assert v.magnitud_v > 0
        assert v.es_linea_linea is False

    def test_get_tension_sin_runpp_falla(self):
        mock = BackendRedMock()
        with pytest.raises(RuntimeError):
            mock.get_tension(1)

    def test_get_tension_bus_slack(self):
        mock = BackendRedMock(v_slack_v=230.0)
        mock.runpp()
        v0 = mock.get_tension(0)
        assert abs(v0.magnitud_v - 230.0) < 1e-6
        assert v0.angulo_rad == 0.0

    def test_get_corriente_linea(self):
        mock = BackendRedMock()
        mock.set_carga(1, 100_000, 0.0)
        mock.runpp()
        i = mock.get_corriente_linea(0)
        # I = P / V ≈ 100k / 230 ≈ 435 A (aproximado)
        # Con el mock corregido, la tensión cae a ~219V, así que I ≈ 100k/219 ≈ 456A
        # Aceptamos rango amplio
        assert i > 400 and i < 600

    def test_context_manager(self):
        with BackendRedMock() as mock:
            mock.set_carga(1, 1000, 500)
            mock.runpp()
            v = mock.get_tension(1)
            assert v.magnitud_v > 0


class TestBackendRedABC:
    def test_no_se_puede_instanciar_directamente(self):
        with pytest.raises(TypeError):
            BackendRed()

    def test_subclase_incompleta_falla(self):
        class Incompleto(BackendRed):
            pass
        with pytest.raises(TypeError):
            Incompleto()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])