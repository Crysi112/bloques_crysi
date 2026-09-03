"""
Tests placeholder para adaptador pandapower (Fase 1).
Se activan cuando pandapower esté instalado.
"""

import pytest

# Intentar importar pandapower y el adaptador
try:
    import pandapower as pp
    from bloques_crysi.red.adaptador_pandapower import BackendPandapower
    PANDAPOWER_DISPONIBLE = True
except ImportError:
    PANDAPOWER_DISPONIBLE = False
    BackendPandapower = None  # type: ignore

pytestmark = pytest.mark.skipif(not PANDAPOWER_DISPONIBLE, reason="pandapower no instalado")


def crear_red_referencia():
    """Red de 2 buses: Slack (0) -- Línea -- Carga (1)"""
    net = pp.create_empty_network()
    b0 = pp.create_bus(net, vn_kv=0.4, name="Bus 0 Slack")
    b1 = pp.create_bus(net, vn_kv=0.4, name="Bus 1 PCC")
    pp.create_ext_grid(net, bus=b0, vm_pu=1.0, va_degree=0.0)
    pp.create_line(net, from_bus=b0, to_bus=b1, length_km=0.1, std_type="NAYY 4x50 SE")
    return net, b0, b1


class TestBackendPandapower:
    def test_inicializacion(self):
        net, _, _ = crear_red_referencia()
        backend = BackendPandapower(net)
        assert backend.nombre_backend == "pandapower"

    def test_set_carga_y_runpp(self):
        net, b0, b1 = crear_red_referencia()
        backend = BackendPandapower(net)
        backend.set_carga(b1, 100_000, 50_000)  # 100 kW, 50 kVAr
        backend.runpp()
        v = backend.get_tension(b1)
        # Validación contra cálculo analítico aproximado
        assert v.magnitud_v > 200  # Debe caer algo desde 230V
        assert v.magnitud_v < 230
        assert v.es_linea_linea is False

    def test_get_tension_bus_slack(self):
        net, b0, b1 = crear_red_referencia()
        backend = BackendPandapower(net)
        backend.runpp()
        v0 = backend.get_tension(b0)
        # Slack bus: 0.4 kV LL = 400/sqrt(3) ≈ 230.94 V LN
        assert abs(v0.magnitud_v - 400.0 / 3**0.5) < 1.0

    def test_set_generacion(self):
        net, b0, b1 = crear_red_referencia()
        backend = BackendPandapower(net)
        backend.set_generacion(b1, 150_000, 0.0)  # 150 kW inyección
        backend.runpp()
        v = backend.get_tension(b1)
        # Con generación, tensión LN debe subir sobre slack (≈230.94 V LN)
        v_slack_ln = 400.0 / 3**0.5
        assert v.magnitud_v > v_slack_ln

    def test_bus_inexistente_falla(self):
        net, b0, b1 = crear_red_referencia()
        backend = BackendPandapower(net)
        with pytest.raises(ValueError):
            backend.set_carga(999, 1000, 500)

    def test_get_tension_sin_runpp_falla(self):
        net, _, _ = crear_red_referencia()
        backend = BackendPandapower(net)
        with pytest.raises(RuntimeError):
            backend.get_tension(0)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])