from .unidades import (
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
from .backend_base import (
    BackendRed,
    BackendRedMock,
    TensionBus,
    PotenciaInyectada,
    ErrorFlujoPotencia,
)
def _importar_adaptador_pandapower():
    try:
        from .adaptador_pandapower import BackendPandapower
        return BackendPandapower
    except ImportError as e:
        raise ImportError(
            "BackendPandapower requiere 'pandapower'. "
            "Instala con: pip install bloques_crysi[red-pandapower]"
        ) from e
def _importar_adaptador_opendss():
    try:
        from .adaptador_opendss import BackendOpenDSS
        return BackendOpenDSS
    except ImportError as e:
        raise ImportError(
            "BackendOpenDSS requiere 'opendssdirect'. "
            "Instala con: pip install bloques_crysi[red-opendss]"
        ) from e
def _importar_red_opendss():
    try:
        from .adaptador_opendss import RedOpenDSS
        return RedOpenDSS
    except ImportError:
        return None
try:
    BackendPandapower = _importar_adaptador_pandapower()
except ImportError:
    BackendPandapower = None
try:
    BackendOpenDSS = _importar_adaptador_opendss()
except ImportError:
    BackendOpenDSS = None
try:
    RedOpenDSS = _importar_red_opendss()
except ImportError:
    RedOpenDSS = None
try:
    from .co_simulador import (
        CoSimuladorRed,
        ResultadoCoSim,
        PCCConfig,
        ResultadoPCC,
        crear_cosimulador_simple,
        comparar_backends,
    )
    from .norton import CoSimNorton, PuertoNorton
except ImportError:
    CoSimuladorRed = None
    ResultadoCoSim = None
    PCCConfig = None
    ResultadoPCC = None
    crear_cosimulador_simple = None
    comparar_backends = None
__all__ = [
    "watts_a_mw", "mw_a_watts", "watts_a_kw", "kw_a_watts",
    "vars_a_mvar", "mvar_a_vars", "vars_a_kvar", "kvar_a_vars",
    "va_a_mva", "mva_a_va", "va_a_kva", "kva_a_va",
    "volts_a_kv", "kv_a_volts",
    "volts_linea_a_fase", "fase_a_volts_linea",
    "rad_a_grados", "grados_a_rad",
    "SistemaUnidades",
    "pq_si_a_pandapower",
    "tension_pandapower_a_si",
    "tension_si_a_opendss",
    "tension_opendss_a_si",
    "pq_si_a_opendss",
    "pq_opendss_a_si",
    "BackendRed",
    "BackendRedMock",
    "TensionBus",
    "PotenciaInyectada",
    "ErrorFlujoPotencia",
    "BackendPandapower",
    "BackendOpenDSS",
    "RedOpenDSS",
    "CoSimuladorRed",
    "ResultadoCoSim",
    "PCCConfig",
    "ResultadoPCC",
    "crear_cosimulador_simple",
    "comparar_backends",
    "CoSimNorton",
    "PuertoNorton",
]
__version__ = "0.2.0"
