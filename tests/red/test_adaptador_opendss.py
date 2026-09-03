"""
Tests para adaptador OpenDSS (Fase 2).
Se activan cuando py-dss-interface esté instalado.
"""

import pytest

try:
    import py_dss_interface
    from bloques_crysi.red.adaptador_opendss import BackendOpenDSS
    OPENDSS_DISPONIBLE = True
except ImportError:
    OPENDSS_DISPONIBLE = False
    BackendOpenDSS = None  # type: ignore

pytestmark = pytest.mark.skipif(not OPENDSS_DISPONIBLE, reason="py-dss-interface no instalado")


def crear_red_referencia_opendss():
    """Red radial simple 2 buses en OpenDSS"""
    dss = py_dss_interface.DSS()
    dss.text("Clear")
    dss.text("New Circuit.Test phases=3 basekv=0.4")
    # Vsource primero (slack bus)
    dss.text("New Vsource.Source bus1=Bus0 basekV=0.4 pu=1.0")
    # Linea
    dss.text("New Linecode.L1 phases=3 r1=0.1 x1=0.05")
    dss.text("New Line.L1 bus1=Bus0 bus2=Bus1 linecode=L1 length=0.1 units=km")
    # Carga
    dss.text("New Load.Load1 bus1=Bus1 phases=3 kV=0.4 kW=100 kvar=50 model=1")
    dss.solution.solve()
    return dss


class TestBackendOpenDSS:
    def test_inicializacion(self):
        dss = crear_red_referencia_opendss()
        backend = BackendOpenDSS(dss, v_base_kv_ll=0.4)
        assert backend.nombre_backend == "OpenDSS"

    def test_set_carga_y_runpp(self):
        dss = crear_red_referencia_opendss()
        backend = BackendOpenDSS(dss, v_base_kv_ll=0.4)
        backend.set_carga(1, 100_000, 50_000)
        backend.runpp()
        v = backend.get_tension(1)
        # Tensión fase-neutro en bus 1 debe caer desde ~230.94V
        v_slack_ln = 400.0 / 3**0.5
        assert v.magnitud_v > 200
        assert v.magnitud_v < v_slack_ln

    def test_get_tension_bus_slack(self):
        dss = crear_red_referencia_opendss()
        backend = BackendOpenDSS(dss, v_base_kv_ll=0.4)
        backend.runpp()
        v0 = backend.get_tension(0)
        assert v0.magnitud_v > 0
        # Slack ~230.94 V LN
        v_slack_ln = 400.0 / 3**0.5
        assert abs(v0.magnitud_v - v_slack_ln) < 1.0

    def test_set_generacion(self):
        dss = crear_red_referencia_opendss()
        backend = BackendOpenDSS(dss, v_base_kv_ll=0.4)
        # 500 kW generación para superar la carga existente de 100 kW
        backend.set_generacion(1, 500_000, 0.0)
        backend.runpp()
        v = backend.get_tension(1)
        # Con generación neta, tensión LN debe subir sobre slack
        v_slack_ln = 400.0 / 3**0.5
        assert v.magnitud_v > v_slack_ln


if __name__ == "__main__":
    pytest.main([__file__, "-v"])