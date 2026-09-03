import ctypes
import os
import pathlib
import shutil
import subprocess
from . import opcodes
CORE_DIR = pathlib.Path(__file__).parent / "core"
SRC_C = CORE_DIR / "block_core.c"
DLL = CORE_DIR / ("bloques_core.dll")
GCC_PATHS = [
    os.environ.get("BLOQUES_CORE_GCC", ""),
    r"C:\msys64\ucrt64\bin\gcc.exe",
    r"C:\msys64\mingw64\bin\gcc.exe",
]
LIB = None
def _locate_gcc():
    for p in GCC_PATHS:
        if p and pathlib.Path(p).exists():
            return p
    which = shutil.which("gcc")
    if which:
        return which
    raise RuntimeError(
        "No se encontró gcc. Instala MSYS2 UCRT64 o define BLOQUES_CORE_GCC."
    )
def compilar(fuerza=False):
    if DLL.exists() and not fuerza:
        if SRC_C.stat().st_mtime <= DLL.stat().st_mtime:
            return DLL
    gcc = _locate_gcc()
    DLL.parent.mkdir(parents=True, exist_ok=True)
    import platform
    extra = ["-lwinmm"] if platform.system() == "Windows" else []
    cmd = [gcc, "-shared", "-O2", "-o", str(DLL), str(SRC_C), "-lm"] + extra
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError("Fallo al compilar la DLL:\n" + proc.stderr)
    return DLL
def _ctypes_modelo():
    return None
class BloqueC(ctypes.Structure):
    _fields_ = [
        ("op", ctypes.c_int),
        ("n_in", ctypes.c_int),
        ("in_idx", ctypes.POINTER(ctypes.c_longlong)),
        ("n_out", ctypes.c_int),
        ("out_idx", ctypes.POINTER(ctypes.c_longlong)),
        ("n_param", ctypes.c_int),
        ("param", ctypes.POINTER(ctypes.c_double)),
        ("n_state", ctypes.c_int),
        ("state", ctypes.POINTER(ctypes.c_double)),
        ("n_ws", ctypes.c_int),
        ("ws", ctypes.POINTER(ctypes.c_double)),
        ("dt", ctypes.c_double),
        ("Ts", ctypes.c_double),
        ("t_next_update", ctypes.c_double),
        ("init", ctypes.c_void_p),
        ("eval_estatico", ctypes.c_void_p),
        ("deriv", ctypes.c_void_p),
        ("out", ctypes.c_void_p),
        ("update", ctypes.c_void_p),
    ]
class ModeloC(ctypes.Structure):
    _fields_ = [
        ("n_bloques", ctypes.c_int),
        ("bloques", ctypes.POINTER(BloqueC)),
        ("n_sig", ctypes.c_int),
        ("sig", ctypes.POINTER(ctypes.c_double)),
        ("n_alg", ctypes.c_int),
        ("alg_list", ctypes.POINTER(ctypes.c_longlong)),
        ("max_iter", ctypes.c_int),
        ("tol", ctypes.c_double),
        ("w_opt", ctypes.c_double),
        ("method", ctypes.c_int),
        ("t", ctypes.c_double),
        ("t_fin", ctypes.c_double),
        ("dt", ctypes.c_double),
        ("error_flag", ctypes.c_int),
    ]
def libreria():
    global LIB
    if LIB is None:
        compilar()
        LIB = ctypes.CDLL(str(DLL))
        LIB.m_sim_run.argtypes = [
            ctypes.POINTER(ModeloC),
            ctypes.c_int,
            ctypes.c_int,
            ctypes.POINTER(ctypes.c_longlong),
            ctypes.POINTER(ctypes.c_double),
        ]
        LIB.m_sim_run.restype = None
        LIB.m_sim_iniciar.argtypes = [ctypes.POINTER(ModeloC)]
        LIB.m_sim_iniciar.restype = ctypes.c_int
        LIB.m_sim_paso.argtypes = [ctypes.POINTER(ModeloC)]
        LIB.m_sim_paso.restype = ctypes.c_int
        LIB.m_sim_guardar.argtypes = [ctypes.POINTER(ModeloC),
                                      ctypes.POINTER(ctypes.c_double)]
        LIB.m_sim_guardar.restype = ctypes.c_int
        LIB.m_sim_restaurar.argtypes = [ctypes.POINTER(ModeloC),
                                        ctypes.POINTER(ctypes.c_double)]
        LIB.m_sim_restaurar.restype = None
        LIB.m_hw_serial_cerrar.argtypes = [ctypes.POINTER(BloqueC)]
        LIB.m_hw_serial_cerrar.restype = None
        LIB.m_hil_ws_size.argtypes = []
        LIB.m_hil_ws_size.restype = ctypes.c_int
    return LIB

def hil_ws_size() -> int:
    return int(libreria().m_hil_ws_size())
