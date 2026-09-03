import pytest
from pathlib import Path
import subprocess
import ctypes
import os

from bloques_crysi import Modelo, GeneradorCodigo
from bloques_crysi.bloques import FuenteConstante, PID, Integrador, Suma, Ganancia

def test_codegen_basico(tmp_path):
    m = Modelo(dt=1e-4)
    with m:
        src = m.add(FuenteConstante('ref', 1.0))
        ctrl = m.add(PID('pid', Kp=2.0, Ki=1.0))
        planta = m.add(Integrador('planta'))
        err = m.add(Suma('err', [1.0, -1.0]))
        
        m.conectar(src.salida, err.entrada[0:1])
        m.conectar(planta.salida, err.entrada[1:2])
        m.conectar(err.salida, ctrl.entrada)
        m.conectar(ctrl.salida, planta.entrada)

    gen = GeneradorCodigo(m)
    c_code = gen.codigo("test_model")
    
    # Verificar que existen las funciones de API publicas
    assert "void test_model_init(void)" in c_code
    assert "void test_model_paso(void)" in c_code
    assert "void test_model_run(int n_steps, double *buf, int n_out, int *idx)" in c_code
    
    # Verificar que los parametros estan macros
    assert "#define _B1P0 2.0" in c_code # Kp = 2.0 del bloque 1 (PID)
    
    # Exportar y tratar de compilar para verificar sintaxis C
    c_file = tmp_path / "test_model.c"
    c_file.write_text(c_code)
    
    # Si gcc está en el path, compilar. Si no, ignorar.
    import shutil
    gcc = shutil.which("gcc")
    if gcc:
        obj_file = tmp_path / "test_model.o"
        res = subprocess.run([gcc, "-c", str(c_file), "-o", str(obj_file)], capture_output=True, text=True)
        assert res.returncode == 0, f"Error de compilación C:\n{res.stderr}"
