import ctypes
from typing import Optional, Sequence, Tuple, Union, List, Dict
import numpy as np
from . import opcodes as ops
from ._clib import BloqueC, ModeloC, libreria, hil_ws_size
from .bloques import Display
from .puertos import Puerto, Sensor
from ._contexto import _get_modelo_actual, _set_modelo_actual, _pop_modelo
class Resultado(dict):
    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self.t = None
    def get(self, key: str, default=None):
        if key not in self:
            return default
        arr = np.asarray(self[key])
        if arr.ndim == 2 and arr.shape[1] == 1:
            return arr.ravel()
        return arr
    def get_2d(self, key: str, default=None):
        if key not in self:
            return default
        arr = np.asarray(self[key])
        if arr.ndim == 1:
            return arr.reshape(-1, 1)
        return arr
    def guardar_csv(self, ruta, separador=","):
        nombres = ["t"]
        cols = [np.asarray(self.t)]
        for k, v in self.items():
            arr = np.asarray(v)
            if arr.ndim == 1:
                nombres.append(str(k))
                cols.append(arr)
            else:
                for j in range(arr.shape[1]):
                    nombres.append(f"{k}[{j}]")
                    cols.append(arr[:, j])
        data = np.column_stack(cols)
        np.savetxt(ruta, data, delimiter=separador, fmt="%.17g",
                   header=separador.join(nombres), comments="")
    def pico(self, key: str, t_max: Optional[float] = None) -> float:
        y = self.get(key)
        if y is None:
            raise KeyError(f"Señal '{key}' no encontrada en Resultado")
        if t_max is not None:
            idx = np.searchsorted(self.t, t_max, side='right')
            y = y[:idx]
        return float(np.max(np.abs(y)))
    def final(self, key: str, ventana: float = 0.0) -> float:
        y = self.get(key)
        if y is None:
            raise KeyError(f"Señal '{key}' no encontrada en Resultado")
        if ventana > 0.0 and self.t is not None:
            idx = np.searchsorted(self.t, self.t[-1] - ventana, side='left')
            y = y[idx:]
        return float(np.mean(y)) if len(y) > 0 else float(y[-1])
    def tiempo_establecimiento(self, key: str, tolerancia: float = 0.02) -> float:
        y = self.get(key)
        if y is None:
            raise KeyError(f"Señal '{key}' no encontrada en Resultado")
        y_final = self.final(key)
        umbral = tolerancia * abs(y_final)
        fuera = np.where(np.abs(y - y_final) > umbral)[0]
        if len(fuera) == 0:
            return 0.0
        idx = fuera[-1]
        return float(self.t[idx]) if self.t is not None else float(idx) * getattr(self, '_dt', 1.0)
    def thd(self, key: str, f0: float = 60.0, t_inicio: float = 0.0, n_armonicos: int = 50) -> float:
        y = self.get(key)
        if y is None:
            raise KeyError(f"Señal '{key}' no encontrada en Resultado")
        if self.t is None:
            raise ValueError("Resultado no tiene eje temporal")
        idx_inicio = np.searchsorted(self.t, t_inicio, side='left')
        y_win = y[idx_inicio:]
        t_win = self.t[idx_inicio:]
        if len(y_win) < 2:
            return 0.0
        n = len(y_win)
        dt = t_win[1] - t_win[0] if n > 1 else 1.0
        freqs = np.fft.rfftfreq(n, dt)
        fft = np.fft.rfft(y_win)
        mag = np.abs(fft)
        idx_f0 = np.argmin(np.abs(freqs - f0))
        f0_real = freqs[idx_f0]
        m0 = mag[idx_f0]
        thd_sq = 0.0
        for k in range(2, n_armonicos + 1):
            fk = k * f0_real
            if fk > freqs[-1]:
                break
            idx_k = np.argmin(np.abs(freqs - fk))
            thd_sq += mag[idx_k] ** 2
        return float(100.0 * np.sqrt(thd_sq) / m0) if m0 > 0 else 0.0
    def fft(self, key: str, t_inicio: float = 0.0) -> Tuple[np.ndarray, np.ndarray]:
        y = self.get(key)
        if y is None:
            raise KeyError(f"Señal '{key}' no encontrada en Resultado")
        if self.t is None:
            raise ValueError("Resultado no tiene eje temporal")
        idx_inicio = np.searchsorted(self.t, t_inicio, side='left')
        y_win = y[idx_inicio:]
        t_win = self.t[idx_inicio:]
        n = len(y_win)
        dt = t_win[1] - t_win[0] if n > 1 else 1.0
        freqs = np.fft.rfftfreq(n, dt)
        fft = np.fft.rfft(y_win)
        return freqs, np.abs(fft)
class Modelo:
    def __init__(self, dt: float = 1e-4, metodo: str = "euler",
                 max_iter: int = 50, tol: float = 1e-9, w_opt: float = 1.0):
        self.dt = float(dt)
        self.metodo = 1 if metodo == "rk4" else 0
        self.max_iter = int(max_iter)
        self.tol = float(tol)
        self.w_opt = float(w_opt)
        self.bloques = []
        self._conexiones = []
        self._nsig = 0
    def __enter__(self):
        _set_modelo_actual(self)
        return self
    def __exit__(self, exc_type, exc_val, exc_tb):
        _pop_modelo()
        return False
    def add(self, bloque):
        if any(b is bloque for b in self.bloques):
            raise ValueError(f"El bloque {bloque.nombre} ya esta en el modelo.")
        if any(b.nombre == bloque.nombre for b in self.bloques):
            raise ValueError(
                f"Ya existe un bloque llamado {bloque.nombre!r} en el "
                "modelo. Los nombres deben ser unicos (revisa si el bloque "
                "se esta creando/agregando dos veces por error)."
            )
        from .bloques import FuncionTransferencia
        if isinstance(bloque, FuncionTransferencia):
            bloque.param = bloque._discretizar(self.dt)
        bloque.in_idx = [-1] * bloque.n_in
        bloque.out_idx = list(range(self._nsig, self._nsig + bloque.n_out))
        self._nsig += bloque.n_out
        self.bloques.append(bloque)
        pendientes = getattr(bloque, "_pendientes_conectar", None)
        if pendientes:
            for origen, destino in pendientes:
                self.conectar(origen, destino)
        return bloque
    def scope(self, nombre: str, *senales, **kwargs):
        from .scope import Scope
        return self.add(Scope(nombre, *senales, **kwargs))
    def indice_tiempo(self, t: float) -> int:
        return int(round(t / self.dt))
    def conectar(self, origen: Puerto, destino: Puerto) -> None:
        if origen.n != destino.n:
            raise ValueError(
                f"Puertos con distinto tamano: {origen.n} != {destino.n}"
            )
        if origen.tipo not in ("sal", "ent") or destino.tipo != "ent":
            raise ValueError("Debe conectar salida -> entrada.")
        self._conexiones.append((origen, destino))
        destino_bloque = getattr(destino, "bloque", None)
        if getattr(destino_bloque, "es_scope", False):
            nombres = getattr(origen, "canales", None)
            if nombres:
                destino_bloque.canales_meta[destino.offset:destino.offset + destino.n] = list(nombres)
    def acoplar_maquinas(self, *maquinas, J_eq: float = 0.1, Bm_eq: float = 0.0,
                         TL_fuente=None):
        from .bloques import EjeMecanico, FuenteConstante
        from .puertos import Puerto
        eje = self.add(EjeMecanico("eje_mecanico_acople", n_maquinas=len(maquinas), J_eq=J_eq, Bm_eq=Bm_eq))
        for i, maq in enumerate(maquinas):
            ext = getattr(maq, "mecanica_interna", None)
            if ext is None:
                ext = maq.n_in < 4
            if ext:
                 raise ValueError(f"La maquina {maq.nombre} parece no tener mecanica externa habilitada.")
            self.conectar(maq.sensorPar().puerto, Puerto(eje, "ent", i, 1))
            self.conectar(eje.salida, maq.puerto_mecanico)
        idx_tl = len(maquinas)
        if TL_fuente is not None:
            self.conectar(TL_fuente.salida, Puerto(eje, "ent", idx_tl, 1))
        else:
            cero = self.add(FuenteConstante("TL_cero", 0.0))
            self.conectar(cero.salida, Puerto(eje, "ent", idx_tl, 1))
        return eje
    def _resolver(self):
        for origen, destino in self._conexiones:
            o = origen.indices()
            for k in range(destino.n):
                destino.bloque.in_idx[destino.offset + k] = o[k]
        for b in self.bloques:
            if hasattr(b, "T_L") and getattr(b, "mecanica_interna", True):
                tl_puerto = b.T_L
                if all(k < 0 for k in tl_puerto.indices()):
                    cero_nombre = f"_TL_cero_{b.nombre}"
                    cero = next((bl for bl in self.bloques if bl.nombre == cero_nombre), None)
                    if cero is None:
                        from .bloques import FuenteConstante
                        cero = self.add(FuenteConstante(cero_nombre, 0.0))
                    self.conectar(cero.salida, tl_puerto)
        for origen, destino in self._conexiones:
            o = origen.indices()
            for k in range(destino.n):
                destino.bloque.in_idx[destino.offset + k] = o[k]
        for b in self.bloques:
            if getattr(b, "es_scope", False):
                continue
            n_opc = getattr(b, "opcionales", 0)
            if b.n_in - n_opc and any(k < 0 for k in b.in_idx[:b.n_in - n_opc]):
                raise ValueError(
                    f"El bloque {b.nombre} tiene entradas sin conectar: {b.in_idx}"
                )
        self._detectar_feedthrough()
    def _detectar_feedthrough(self):
        def _tiene_feedthrough(b):
            if b.op in ops.ES_ESTATICO:
                return True
            if b.op == ops.OP_TF:
                return (len(b._num) == len(b._den)
                        and b.param[1] != 0.0)
            if b.op == ops.OP_PID:
                return b.param[0] != 0.0 or b.param[2] != 0.0
            return False
        es_feed = {id(b): _tiene_feedthrough(b) for b in self.bloques
                   if not getattr(b, "es_scope", False)}
        prod = {}
        for b in self.bloques:
            for k in b.out_idx:
                prod[k] = b
        grafo = {id(b): set() for b in self.bloques
                 if not getattr(b, "es_scope", False) and es_feed[id(b)]}
        for b in self.bloques:
            if id(b) not in grafo:
                continue
            for k in b.in_idx:
                if k >= 0 and k in prod and id(prod[k]) in grafo \
                        and prod[k] is not b:
                    grafo[id(b)].add(id(prod[k]))
        grado = {i: len(hijos) for i, hijos in grafo.items()}
        hijos = {i: set() for i in grafo}
        for i, padres in grafo.items():
            for p in padres:
                hijos[p].add(i)
        from collections import deque
        cola = deque(i for i in grafo if grado[i] == 0)
        while cola:
            i = cola.popleft()
            for h in hijos[i]:
                grado[h] -= 1
                if grado[h] == 0:
                    cola.append(h)
        en_ciclo = [b for b in self.bloques if not getattr(b, "es_scope", False)
                    and id(b) in grafo and grado[id(b)] > 0]
        if not any(b.op in (ops.OP_TF, ops.OP_PID) for b in en_ciclo):
            return
        nombres = ", ".join(sorted(b.nombre for b in en_ciclo))
        raise ValueError(
            "Lazo algebraico oculto a traves de bloques con feedthrough "
            f"directo: {nombres}. Los bloques dinamicos con feedthrough "
            "(PID con Kp/Kd, TF con numerador y denominador de igual grado) "
            "no pueden quedar dentro de un lazo sin estados. Agrega un "
            "bloque con estado (Integrador, maquina, filtro) para romper "
            "el lazo."
        )
    def _orden_estatico(self):
        estaticos = [i for i, b in enumerate(self.bloques) if b.op in ops.ES_ESTATICO]
        if not estaticos:
            return []
        prod = {}
        for i in estaticos:
            for k in self.bloques[i].out_idx:
                prod[k] = i
        hijos = {i: set() for i in estaticos}
        grado = {i: 0 for i in estaticos}
        for i in estaticos:
            for k in self.bloques[i].in_idx:
                if k in prod and prod[k] != i:
                    hijos[prod[k]].add(i)
        for i in estaticos:
            grado[i] = len([p for p in estaticos if i in hijos[p] and p != i])
        cola = [i for i in estaticos if grado[i] == 0]
        orden = []
        from collections import deque
        cola = deque(cola)
        while cola:
            i = cola.popleft()
            orden.append(i)
            for j in hijos[i]:
                grado[j] -= 1
                if grado[j] == 0:
                    cola.append(j)
        if len(orden) != len(estaticos):
            return estaticos
        return orden
    def _armar_modelo_c(self, t_fin: float, registrar: Sequence):
        self._resolver()
        bloques = [b for b in self.bloques if not getattr(b, "es_scope", False)]
        scopes = [b for b in self.bloques if getattr(b, "es_scope", False)]
        n = len(bloques)
        n_sig = self._nsig
        arr_in = []
        arr_out = []
        arr_param = []
        arr_state = []
        arr_ws = []
        for b in bloques:
            nin = max(b.n_in, 1)
            nout = max(b.n_out, 1)
            arr_in.append((ctypes.c_longlong * nin)(*(b.in_idx + [-1] * (nin - b.n_in))))
            arr_out.append((ctypes.c_longlong * nout)(*(b.out_idx + [0] * (nout - b.n_out))))
            npam = max(len(b.param), 1)
            arr_param.append((ctypes.c_double * npam)(*([float(x) for x in b.param] + [0.0] * (npam - len(b.param)))))
            nst = max(b.n_state, 1)
            arr_state.append((ctypes.c_double * nst)(*([float(x) for x in b.estados_iniciales] + [0.0] * (nst - b.n_state))))
            _hil_ws = hil_ws_size() if b.op == int(ops.OP_HW_SERIAL) else 5 * nst
            arr_ws.append((ctypes.c_double * _hil_ws)(*([0.0] * _hil_ws)))
        self._param_arrays = arr_param
        self._in_arrays = arr_in
        self._out_arrays = arr_out
        self._state_arrays = arr_state
        self._ws_arrays = arr_ws
        self._bloques_activos = bloques
        bloques_c = (BloqueC * n)()
        self._bloques_c = bloques_c
        self._n_bloques_c = n
        for i, b in enumerate(bloques):
            bc = bloques_c[i]
            bc.op = int(b.op)
            bc.n_in = b.n_in
            bc.in_idx = arr_in[i]
            bc.n_out = b.n_out
            bc.out_idx = arr_out[i]
            bc.n_param = len(b.param)
            bc.param = arr_param[i]
            bc.n_state = b.n_state
            bc.state = arr_state[i]
            bc.n_ws = hil_ws_size() if b.op == int(ops.OP_HW_SERIAL) else 5 * max(b.n_state, 1)
            bc.ws = arr_ws[i]
            bc.dt = self.dt
            bc.Ts = getattr(b, "Ts", 0.0)
            bc.t_next_update = 0.0
        grabar = []
        if not registrar:
            registrar = tuple(bloques)
        for item in registrar:
            if isinstance(item, tuple):
                raise TypeError("registrar acepta sensores, puertos o bloques.")
            if isinstance(item, bytes) or isinstance(item, str):
                self._grabar_por_nombre(item, grabar)
            elif getattr(item, "es_scope", False):
                idx = item._indices_conectados()
                if not idx:
                    raise ValueError(f"El Scope {item.nombre} no tiene senales conectadas.")
                grabar.append((item.nombre, item._etiquetas(len(idx)), idx))
            elif hasattr(item, "puerto"):
                p = item.puerto
                idx = p.indices()
                grabar.append((item.nombre, item.canales, [int(k) for k in idx]))
            elif hasattr(item, "tipo"):
                idx = item.indices()
                canales = [f"{item.bloque.nombre}[{item.offset + k}]" for k in range(len(idx))]
                grabar.append((item.bloque.nombre, canales, [int(k) for k in idx]))
            elif hasattr(item, "out_idx"):
                canales = getattr(item, "NOMBRES", None)
                if canales is None:
                    canales = [f"{item.nombre}[{k}]" for k in range(item.n_out)]
                grabar.append((item.nombre, list(canales), [int(k) for k in item.out_idx]))
            else:
                raise TypeError(f"No se grabar {item!r}.")
        d_rec = []
        total = 0
        for b in self.bloques:
            if isinstance(b, Display):
                d_rec.append((b.nombre, b.formato, b.out_idx[0]))
                total += 1
        if total == 0 and not grabar:
            raise ValueError("No hay senales que registrar.")
        n_steps = int(round(t_fin / self.dt)) + 1
        sig = (ctypes.c_double * n_sig)(*([0.0] * n_sig))
        self._sig_array = sig
        orden = self._orden_estatico()
        pos = {id(b): i for i, b in enumerate(bloques)}
        orden_c = [pos[id(self.bloques[k])] for k in orden]
        n_alg = len(orden_c)
        alg = (ctypes.c_longlong * max(n_alg, 1))(*([int(k) for k in orden_c] + [-1] * (max(n_alg, 1) - n_alg)))
        modelo = ModeloC()
        modelo.n_bloques = n
        modelo.bloques = bloques_c
        modelo.n_sig = n_sig
        modelo.sig = sig
        modelo.n_alg = n_alg
        modelo.alg_list = alg
        modelo.max_iter = self.max_iter
        modelo.tol = self.tol
        modelo.w_opt = self.w_opt
        modelo.method = self.metodo
        modelo.t = 0.0
        modelo.t_fin = float(t_fin)
        modelo.dt = self.dt
        rec_flat = [int(k) for g in grabar for k in g[2]]
        rec_flat += [int(idx) for _, _, idx in d_rec]
        n_rec = len(rec_flat)
        rec_idx = (ctypes.c_longlong * n_rec)(*rec_flat)
        return modelo, grabar, d_rec, rec_idx, n_steps, scopes
    def _desempaquetar(self, grabar, d_rec, datos, t):
        res = Resultado()
        res.t = np.asarray(t, dtype=float)
        off = 0
        for nombre, canales, idxs in grabar:
            if len(canales) < 2:
                res[nombre] = datos[off] * 1.0
            else:
                res[nombre] = datos[off:off + len(canales)].T * 1.0
            off += len(canales)
        for nombre, formato, idx in d_rec:
            print(f"{nombre}: {formato % datos[off, -1]}")
            off += 1
        return res
    @staticmethod
    def _incluir_scopes(registrar, bloques):
        registrar = list(registrar)
        for b in bloques:
            if getattr(b, "es_scope", False) and b not in registrar:
                registrar.append(b)
        return registrar
    def run(self, t_fin: float, registrar: Sequence = (),
            retorna_records: bool = False) -> "Resultado":
        registrar = self._incluir_scopes(registrar, self.bloques)
        rt = [s for s in registrar if getattr(s, "tiempo_real", False)]
        if rt and not retorna_records:
            chunk = max(10 * self.dt, t_fin / 50.0)
            parciales = []
            for parcial in self.iterar(t_fin, registrar=registrar,
                                       chunk=chunk):
                parciales.append(parcial)
                for sc in rt:
                    sc.actualizar(parcial.t, parcial[sc.nombre])
            res = self._concatenar_resultados(parciales)
            for sc in rt:
                sc.esperar()
            return res
        modelo, grabar, d_rec, rec_idx, n_steps, scopes = \
            self._armar_modelo_c(t_fin, registrar)
        total = len(rec_idx)
        rec_buf = (ctypes.c_double * (total * n_steps))()
        lib = libreria()
        lib.m_sim_run(ctypes.byref(modelo), n_steps, total, rec_idx, rec_buf)
        if modelo.error_flag:
            raise RuntimeError(
                "El lazo algebraico no convergio "
                f"(t = {modelo.t:.6g} s, max_iter = {self.max_iter}, "
                f"tol = {self.tol:g}). Revisa el modelo o ajusta "
                "Modelo(max_iter=..., tol=..., w_opt=...)."
            )
        datos = np.frombuffer(rec_buf, dtype=np.float64).reshape(total, n_steps)
        res = self._desempaquetar(grabar, d_rec, datos,
                                  np.arange(n_steps) * self.dt)
        if retorna_records:
            return res, modelo
        visibles = [sc for sc in scopes
                    if sc.mostrar and not getattr(sc, "tiempo_real", False)]
        for i, sc in enumerate(visibles):
            sc.mostrar_grafico(res, ahora=(len(visibles) == 1))
        if len(visibles) > 1:
            import matplotlib.pyplot as plt
            plt.show(block=any(sc.bloqueo for sc in visibles))
        return res
    @staticmethod
    def _concatenar_resultados(parciales):
        res = Resultado()
        res.t = np.concatenate([parciales[0].t] +
                                [p.t[1:] for p in parciales[1:]])
        for nombre in parciales[0]:
            if nombre == "t":
                continue
            v0 = parciales[0][nombre]
            if np.asarray(v0).ndim == 0:
                res[nombre] = v0
                continue
            res[nombre] = np.concatenate(
                [v0] + [np.asarray(p[nombre])[1:] for p in parciales[1:]])
        return res
    def iniciar(self, registrar: Sequence = ()):
        registrar = self._incluir_scopes(registrar, self.bloques)
        modelo, grabar, d_rec, rec_idx, _, _ = \
            self._armar_modelo_c(0.0, registrar)
        self._paso_ctx = (modelo, grabar, d_rec, rec_idx)
        libreria().m_sim_iniciar(ctypes.byref(modelo))
        return self._valores_actuales()
    def paso(self):
        if not hasattr(self, "_paso_ctx"):
            raise RuntimeError("Llama primero a Modelo.iniciar(registrar=...).")
        modelo, grabar, d_rec, rec_idx = self._paso_ctx
        if libreria().m_sim_paso(ctypes.byref(modelo)):
            raise RuntimeError(
                "El lazo algebraico no convergio "
                f"(t = {modelo.t:.6g} s, max_iter = {self.max_iter}, "
                f"tol = {self.tol:g}). Revisa el modelo o ajusta "
                "Modelo(max_iter=..., tol=..., w_opt=...)."
            )
        return self._valores_actuales()
    def _tam_foto(self, modelo):
        total_estados = sum(int(modelo.bloques[i].n_state)
                            for i in range(modelo.n_bloques))
        return 1 + int(modelo.n_sig) + sum(
            int(modelo.bloques[i].n_state) + int(modelo.bloques[i].n_ws) + 1
            for i in range(modelo.n_bloques))
    def guardar_estado(self, buffer=None):
        if not hasattr(self, "_paso_ctx"):
            raise RuntimeError("Llama primero a Modelo.iniciar(registrar=...).")
        modelo = self._paso_ctx[0]
        if buffer is None:
            buffer = (ctypes.c_double * max(self._tam_foto(modelo), 1))()
        libreria().m_sim_guardar(ctypes.byref(modelo), buffer)
        return buffer
    def restaurar_estado(self, buffer=None):
        if not hasattr(self, "_paso_ctx"):
            raise RuntimeError("Llama primero a Modelo.iniciar(registrar=...).")
        modelo = self._paso_ctx[0]
        if buffer is None:
            buffer = (ctypes.c_double * max(self._tam_foto(modelo), 1))()
        libreria().m_sim_restaurar(ctypes.byref(modelo), buffer)
    def _valores_actuales(self):
        modelo, grabar, d_rec, rec_idx = self._paso_ctx
        sig = np.ctypeslib.as_array(modelo.sig, shape=(modelo.n_sig,))
        out = {}
        for nombre, canales, idxs in grabar:
            v = sig[list(idxs)]
            out[nombre] = float(v[0]) if len(v) == 1 else np.array(v, dtype=float)
        for nombre, formato, idx in d_rec:
            print(f"{nombre}: {formato % sig[idx]}")
        return out
    def iterar(self, t_fin: float, registrar: Sequence = (),
               chunk: Optional[float] = None, eventos=None,
               profundidad: int = 12):
        registrar = self._incluir_scopes(registrar, self.bloques)
        if chunk is None:
            yield self.run(t_fin, registrar=registrar)
            return
        modelo, grabar, d_rec, rec_idx, n_steps, _ = \
            self._armar_modelo_c(t_fin, registrar)
        lib = libreria()
        lib.m_sim_iniciar(ctypes.byref(modelo))
        sig = np.ctypeslib.as_array(modelo.sig, shape=(modelo.n_sig,))
        idx = np.array([int(k) for k in rec_idx], dtype=np.int64)
        n_can = len(idx)
        ev_list = []
        for item, umbral in (eventos or []):
            p = item.puerto if hasattr(item, "puerto") else item
            ev_list.append((int(p.indices()[0]), float(umbral), p.bloque.nombre))
        n_chunk = max(1, int(round(chunk / self.dt)))
        n_steps = max(n_steps, 1)
        if n_steps <= 1:
            buf = np.zeros((n_can, 1))
            buf[:, 0] = sig[idx]
            yield self._desempaquetar(grabar, d_rec, buf, np.zeros(1))
            return
        total_estados = sum(int(modelo.bloques[i].n_state)
                            for i in range(modelo.n_bloques))
        tam_foto = 1 + int(modelo.n_sig) + sum(
            int(modelo.bloques[i].n_state) + int(modelo.bloques[i].n_ws) + 1
            for i in range(modelo.n_bloques))
        buf0 = (ctypes.c_double * max(tam_foto, 1))()
        k = 1
        t0 = 0.0
        while k < n_steps:
            n_win = min(n_chunk, n_steps - k + 1)
            cols = [sig[idx].copy()]
            ts = [t0]
            t_cur = t0
            for _ in range(n_win - 1):
                v_ant = {ev[0]: float(sig[ev[0]]) for ev in ev_list}
                lib.m_sim_guardar(ctypes.byref(modelo), buf0)
                if lib.m_sim_paso(ctypes.byref(modelo)):
                    raise RuntimeError(
                        "El lazo algebraico no convergio "
                        f"(t = {modelo.t:.6g} s, max_iter = {self.max_iter}, "
                        f"tol = {self.tol:g}). Revisa el modelo o ajusta "
                        "Modelo(max_iter=..., tol=..., w_opt=...)."
                    )
                t_ant, t_cur = t_cur, modelo.t
                for ev_idx, umbral, nombre in ev_list:
                    va, vb = v_ant[ev_idx], float(sig[ev_idx])
                    if (va - umbral) * (vb - umbral) < 0:
                        ref = self._bisecar_evento(
                            modelo, sig, idx, ev_idx, umbral, va,
                            t_ant, t_cur, profundidad, buf0)
                        if ref is not None:
                            cols.append(ref[1])
                            ts.append(ref[0])
                cols.append(sig[idx].copy())
                ts.append(t_cur)
            k += n_win - 1
            yield self._desempaquetar(grabar, d_rec,
                                      np.stack(cols, axis=1), ts)
            t0 += (n_win - 1) * self.dt
        return
    def _bisecar_evento(self, modelo, sig, idx, ev_idx, umbral, v_ant,
                        t_ant, t_cur, profundidad, buf0):
        lib = libreria()
        dt_orig = self.dt
        def revertir(mensaje):
            lib.m_sim_restaurar(ctypes.byref(modelo), buf0)
            modelo.t = t_ant
            modelo.dt = dt_orig
            if lib.m_sim_paso(ctypes.byref(modelo)):
                raise RuntimeError(
                    "El lazo algebraico no convergio al restaurar un "
                    f"evento (t = {modelo.t:.6g} s). {mensaje}")
            return None
        a, b = t_ant, t_cur
        fa = v_ant - umbral
        if abs(fa) < 1e-15:
            return None
        fb = float(sig[ev_idx]) - umbral
        if abs(fb) < 1e-15:
            return None
        for _ in range(profundidad):
            lib.m_sim_restaurar(ctypes.byref(modelo), buf0)
            modelo.t = a
            modelo.dt = 0.5 * (b - a)
            if lib.m_sim_paso(ctypes.byref(modelo)):
                return revertir("Prueba de biseccion sin convergencia.")
            fm = float(sig[ev_idx]) - umbral
            if fa * fm <= 0:
                b = modelo.t
            else:
                a, fa = modelo.t, fm
        tc = 0.5 * (a + b)
        dt1 = tc - t_ant
        dt2 = t_cur - tc
        if dt1 < 1e-15 or dt2 < 1e-15:
            return revertir("")
        lib.m_sim_restaurar(ctypes.byref(modelo), buf0)
        modelo.t = t_ant
        modelo.dt = dt1
        if lib.m_sim_paso(ctypes.byref(modelo)):
            return revertir("Minipaso hasta t_cruce sin convergencia.")
        vals = sig[idx].copy()
        modelo.dt = 0.0
        if lib.m_sim_paso(ctypes.byref(modelo)):
            return revertir("Conmutacion en t_cruce sin convergencia.")
        modelo.dt = dt2
        if lib.m_sim_paso(ctypes.byref(modelo)):
            return revertir("Minipaso final sin convergencia.")
        modelo.dt = dt_orig
        return (tc, vals)
    def _grabar_por_nombre(self, nombre, grabar):
        encontrado = [b for b in self.bloques if b.nombre == nombre]
        if not encontrado:
            raise ValueError(f"No existe un bloque llamado {nombre!r}.")
        b = encontrado[0]
        canales = getattr(b, "NOMBRES", None)
        if canales is None:
            canales = [f"{nombre}[{k}]" for k in range(b.n_out)]
        grabar.append((nombre, list(canales), [int(k) for k in b.out_idx]))
    def set_param(self, bloque, indice: int, valor: float) -> None:
        if not hasattr(self, "_param_arrays") or self._param_arrays is None:
            raise ValueError(
                "El modelo no ha sido armado aun. Llama run(), iniciar() o iterar() "
                "antes de usar set_param()."
            )
        if isinstance(bloque, str):
            encontrados = [b for b in self._bloques_activos if b.nombre == bloque]
            if not encontrados:
                raise ValueError(f"No existe un bloque activo llamado {bloque!r}.")
            b = encontrados[0]
        else:
            b = bloque
            if b not in self._bloques_activos:
                raise ValueError("El bloque no pertenece a este modelo o es un Scope.")
        n_params = len(b.param)
        if not (0 <= indice < n_params):
            raise ValueError(
                f"Indice de parametro {indice} fuera de rango para bloque {b.nombre} "
                f"(tiene {n_params} parametros, indices 0..{n_params-1})."
            )
        b.param[indice] = float(valor)
        idx_activo = self._bloques_activos.index(b)
        self._param_arrays[idx_activo][indice] = float(valor)
    def get_param(self, bloque, indice: int) -> float:
        if isinstance(bloque, str):
            encontrados = [b for b in self._bloques_activos if b.nombre == bloque]
            if not encontrados:
                raise ValueError(f"No existe un bloque activo llamado {bloque!r}.")
            b = encontrados[0]
        else:
            b = bloque
            if b not in self._bloques_activos:
                raise ValueError("El bloque no pertenece a este modelo o es un Scope.")
        n_params = len(b.param)
        if not (0 <= indice < n_params):
            raise ValueError(
                f"Indice de parametro {indice} fuera de rango para bloque {b.nombre}."
            )
        idx_activo = self._bloques_activos.index(b)
        return self._param_arrays[idx_activo][indice]
    def cerrar_hw(self):
        lib = libreria()
        for i in range(self._n_bloques_c):
            bc = self._bloques_c[i]
            if bc.op == int(ops.OP_HW_SERIAL):
                lib.m_hw_serial_cerrar(ctypes.byref(bc))
    def __repr__(self):
        lineas = [f"<Modelo dt={self.dt} bloques={[b.nombre for b in self.bloques]}>"]
        return "\n".join(lineas)
    def acoplar_red(self, backend, bus_pcc: Union[int, Sequence[int]],
                    elemento: Union[object, Sequence[object]],
                    v_nominal_ll: Union[float, Sequence[float]] = 400.0,
                    dt_red: float = 0.1,
                    es_generacion: Union[bool, Sequence[bool]] = False,
                    fases: Union[str, Sequence[str]] = "ABC",
                    v_nominal_ln: Optional[Union[float, Sequence[float]]] = None,
                    tol_convergencia_v: float = 1e-3,
                    max_iter_ventana: int = 20,
                    relajacion: float = 0.5,
                    reemplazar_cargas: bool = True) -> "CoSimuladorRed":
        if not isinstance(bus_pcc, (list, tuple)):
            bus_pcc = [bus_pcc]
        if not isinstance(elemento, (list, tuple)):
            elemento = [elemento]
        n = len(bus_pcc)
        if len(elemento) != n:
            raise ValueError("bus_pcc y elemento deben tener la misma longitud.")
        if any(isinstance(b, str) for b in bus_pcc):
            norm = []
            for b in bus_pcc:
                if isinstance(b, str):
                    base = str(b).split('.')[0].lower()
                    if hasattr(backend, "_bus_to_idx") and base in backend._bus_to_idx:
                        norm.append(backend._bus_to_idx[base])
                    elif hasattr(backend, "_net") and hasattr(backend._net, "bus"):
                        mask = backend._net.bus.name.astype(str).str.contains(base, na=False)
                        if mask.any():
                            norm.append(int(backend._net.bus[mask].index[0]))
                        else:
                            norm.append(b)
                    else:
                        norm.append(b)
                else:
                    norm.append(b)
            bus_pcc = norm
        if reemplazar_cargas:
            buses = bus_pcc if isinstance(bus_pcc, (list, tuple)) else [bus_pcc]
            for b in buses:
                try:
                    backend.desactivar_cargas_estaticas(b)
                except Exception:
                    pass
        if not isinstance(es_generacion, (list, tuple)):
            es_generacion = [es_generacion] * n
        if not isinstance(v_nominal_ll, (list, tuple)):
            v_nominal_ll = [v_nominal_ll] * n
        if not isinstance(fases, (list, tuple)):
            fases = [fases] * n
        if v_nominal_ln is None:
            v_nominal_ln = [v / np.sqrt(3.0) for v in v_nominal_ll]
        elif not isinstance(v_nominal_ln, (list, tuple)):
            v_nominal_ln = [v_nominal_ln] * n
        if len(es_generacion) != n or len(v_nominal_ll) != n or len(v_nominal_ln) != n or len(fases) != n:
            raise ValueError("Las listas de configuración (es_generacion, v_nominal_ll, fases) "
                             "deben tener la misma longitud que bus_pcc o ser escalares.")
        from .bloques import FuenteTrifasica, MedidorPotencia, Demultiplexor
        from .red import CoSimuladorRed, PCCConfig
        pccs = []
        for bus, elem, es_gen, v_ll, v_ln, f_str in zip(bus_pcc, elemento, es_generacion, v_nominal_ll, v_nominal_ln, fases):
            f_str = f_str.upper()
            amp = v_ll * np.sqrt(2.0 / 3.0)
            fuente = self.add(FuenteTrifasica(f"red_pcc_{bus}_{f_str}", amplitud=amp, frecuencia=50.0))
            puerto_v = getattr(elem, "terminales", getattr(elem, "entrada", None))
            if puerto_v is None:
                raise ValueError(f"No se detectó puerto de entrada de tensión para {elem.nombre}")
            if hasattr(elem, "sensor3I"):
                puerto_i = elem.sensor3I().puerto
            else:
                puerto_i = elem.salida
            if f_str == "ABC":
                medidor = self.add(MedidorPotencia(f"med_pcc_{bus}_ABC", fases=3))
                self.conectar(fuente.salida, medidor.entrada)
                self.conectar(fuente.salida, puerto_v)
                self.conectar(puerto_i, medidor.corrientes)
            else:
                idx_fase = {"A": 0, "B": 1, "C": 2}[f_str]
                dmx_v = self.add(Demultiplexor(f"dmx_v_pcc_{bus}_{f_str}", 3))
                self.conectar(fuente.salida, dmx_v.entrada)
                medidor = self.add(MedidorPotencia(f"med_pcc_{bus}_{f_str}", fases=1))
                self.conectar(dmx_v.salidas[idx_fase], medidor.entrada)
                puerto_v_1f = elem.entrada
                if puerto_v_1f is None or puerto_v_1f.n != 1:
                    raise ValueError(f"Elemento {elem.nombre} no tiene entrada monofásica válida (n={puerto_v_1f.n if puerto_v_1f else None})")
                self.conectar(dmx_v.salidas[idx_fase], puerto_v_1f)
                self.conectar(puerto_i, medidor.corrientes)
            pccs.append(PCCConfig(bus_pcc=bus, medidor=medidor, fuente_red=fuente,
                                  v_nominal_ln=v_ln, es_generacion=es_gen, fases=f_str))
        return CoSimuladorRed(
            modelo=self, backend=backend, pccs=pccs,
            dt_red=dt_red, tol_convergencia_v=tol_convergencia_v,
            max_iter_ventana=max_iter_ventana, relajacion=relajacion
        )
