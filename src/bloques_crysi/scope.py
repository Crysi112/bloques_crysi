import numpy as np
from typing import Optional, Sequence, Tuple
from .puertos import Puerto

def _decimar_minmax(t: np.ndarray, y: np.ndarray, max_puntos: int) -> Tuple[np.ndarray, np.ndarray]:
    n = len(t)
    if n <= max_puntos:
        return t, y
    n_bins = max_puntos // 2
    if n_bins < 1:
        n_bins = 1
    bin_size = int(np.ceil(n / n_bins))
    actual_bins = n // bin_size
    if actual_bins == 0:
        actual_bins = 1
        bin_size = n
    elif actual_bins > n_bins:
        actual_bins = n_bins
    n_use = actual_bins * bin_size
    t_bin = t[:n_use].reshape(actual_bins, bin_size)

    if y.ndim == 1:
        y_bin = y[:n_use].reshape(actual_bins, bin_size)
        y_min = y_bin.min(axis=1)
        y_max = y_bin.max(axis=1)
        t_dec = np.empty(actual_bins * 2)
        y_dec = np.empty(actual_bins * 2)
        t_dec[0::2] = t_bin[:, 0]
        t_dec[1::2] = t_bin[:, -1]
        y_dec[0::2] = y_min
        y_dec[1::2] = y_max
    else:
        m = y.shape[1]
        y_bin = y[:n_use, :].reshape(actual_bins, bin_size, m)
        y_min = y_bin.min(axis=1)
        y_max = y_bin.max(axis=1)
        t_dec = np.empty(actual_bins * 2)
        y_dec = np.empty((actual_bins * 2, m))
        t_dec[0::2] = t_bin[:, 0]
        t_dec[1::2] = t_bin[:, -1]
        y_dec[0::2, :] = y_min
        y_dec[1::2, :] = y_max
    return t_dec, y_dec

def _decimar_lttb(t: np.ndarray, y: np.ndarray, max_puntos: int) -> Tuple[np.ndarray, np.ndarray]:
    n = len(t)
    if n <= max_puntos:
        return t, y
    step = max(1, n // max_puntos)
    idx = np.arange(0, n, step)
    if len(idx) > max_puntos:
        idx = idx[:max_puntos]
    if idx[-1] != n - 1:
        idx = np.append(idx, n - 1)
    return t[idx], y[idx]

def decimar_datos(t: np.ndarray, y: np.ndarray, max_puntos: int = 10000,
                  metodo: str = "step") -> Tuple[np.ndarray, np.ndarray]:
    if len(t) <= max_puntos:
        return t, y
    if metodo == "step":
        step = max(1, len(t) // max_puntos)
        return t[::step], y[::step]
    if metodo == "lttb":
        return _decimar_lttb(t, y, max_puntos)
    return _decimar_minmax(t, y, max_puntos)

FORMATO_SENALES = {
    "ia": r"$i_{as}\ [\mathrm{A}]$",
    "ib": r"$i_{bs}\ [\mathrm{A}]$",
    "ic": r"$i_{cs}\ [\mathrm{A}]$",
    "ias": r"$i_{as}\ [\mathrm{A}]$",
    "ibs": r"$i_{bs}\ [\mathrm{A}]$",
    "ics": r"$i_{cs}\ [\mathrm{A}]$",
    "i'ar": r"$i'_{ar}\ [\mathrm{A}]$",
    "i'br": r"$i'_{br}\ [\mathrm{A}]$",
    "i'cr": r"$i'_{cr}\ [\mathrm{A}]$",
    "Te": r"$T_e\ [\mathrm{N\cdot m}]$",
    "wm": r"$\omega_m\ [\mathrm{rad/s}]$",
    "rpm": r"$\mathrm{Speed}\ [\mathrm{r/min}]$",
    "P_cu_s": r"$P_{\mathrm{Cu}}\ [\mathrm{W}]$",
    "P_cu": r"$P_{\mathrm{Cu}}\ [\mathrm{W}]$",
    "T": r"$T\ [^\circ\mathrm{C}]$",
    "va": r"$v_{as}\ [\mathrm{V}]$",
    "vb": r"$v_{bs}\ [\mathrm{V}]$",
    "vc": r"$v_{cs}\ [\mathrm{V}]$",
    "vdc": r"$V_{\mathrm{dc}}\ [\mathrm{V}]$",
    "Vbat": r"$V_{\mathrm{bat}}\ [\mathrm{V}]$",
}

def formalizar_etiqueta(texto: str) -> str:
    if texto in FORMATO_SENALES:
        return FORMATO_SENALES[texto]
    if "[" in texto or "$" in texto:
        return texto
    return texto

class Scope:
    es_scope = True
    op = None
    n_out = 0
    n_state = 0

    def __init__(self, nombre: str, *senales, max_canales: int = 1,
                 anchos: Optional[Sequence[int]] = None, mostrar: bool = True,
                 guiones: Optional[Sequence[str]] = None,
                 titulo: Optional[str] = None, bloqueo: bool = True,
                 xy_mode: Optional[Tuple[int, int]] = None,
                 superponer_canales: bool = True,
                 max_puntos: int = 10000,
                 cuadricula: Optional[Tuple[int, int]] = None):
        if senales:
            self._senales_a_conectar = [
                s.puerto if hasattr(s, "puerto") else s for s in senales
            ]
            anchos = [len(s) for s in self._senales_a_conectar]
        else:
            self._senales_a_conectar = []

        if anchos is None:
            anchos = [1] * max_canales

        self.nombre = nombre
        self.anchos = [int(a) for a in anchos]
        self.max_canales = len(self.anchos)
        self.n_in = sum(self.anchos)
        offset = 0
        self.canales = []
        for w in self.anchos:
            self.canales.append(Puerto(self, "ent", offset, w))
            offset += w

        self.in_idx = [-1] * self.n_in
        self.out_idx = []
        self.mostrar = bool(mostrar)
        self.bloqueo = bool(bloqueo)
        self.guiones = list(guiones) if guiones else None
        self.canales_meta = [None] * self.n_in
        self.titulo = titulo or nombre
        self.xy_mode = xy_mode
        self.superponer_canales = bool(superponer_canales)
        self.max_puntos = int(max_puntos)
        self.cuadricula = cuadricula
        self._pendientes_conectar = [
            (s, self.canales[k]) for k, s in enumerate(self._senales_a_conectar)
        ]

    def __repr__(self):
        return f"<Scope {self.nombre!r} ({self.max_canales} canales)>"

    def _indices_conectados(self):
        return [int(k) for k in self.in_idx if k >= 0]

    def _etiquetas(self, n):
        if self.guiones:
            if len(self.guiones) == n:
                return list(self.guiones)
            if len(self.guiones) == len(self.canales):
                out = [g for g, w in zip(self.guiones, self.anchos) for _ in range(w)]
                if len(out) == n:
                    return out
        meta = [c for c in self.canales_meta if c is not None]
        if len(meta) == n:
            return meta
        return [f"{self.nombre}[{k}]" for k in range(n)]

    def datos(self, res):
        t = np.asarray(res.t)
        d = np.asarray(res[self.nombre])
        if d.ndim == 1:
            d = d.reshape(-1, 1)
        return t, d

    def mostrar_grafico(self, res, ahora=True):
        if not self.mostrar:
            return
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            return

        plt.rcParams.update({
            "font.family": "serif",
            "font.serif": ["Times New Roman", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "font.size": 22,
            "axes.linewidth": 1.4,
            "lines.linewidth": 2.8,
            "xtick.major.width": 1.4,
            "ytick.major.width": 1.4,
            "xtick.labelsize": 22,
            "ytick.labelsize": 22,
            "axes.labelsize": 24,
        })

        t, d = self.datos(res)

        if self.xy_mode is not None:
            if len(t) > self.max_puntos:
                step = max(1, len(t) // self.max_puntos)
                t, d = t[::step], d[::step]
        else:
            if len(t) > self.max_puntos:
                t, d = decimar_datos(t, d, max_puntos=self.max_puntos, metodo="step")

        n = d.shape[1]
        if not n:
            return

        etiquetas = self._etiquetas(n)
        colores = ["blue", "red", "limegreen", "magenta", "black", "darkorange", "purple", "cyan", "brown"]

        if self.xy_mode is not None:
            ix, iy = self.xy_mode
            fig, ax = plt.subplots(figsize=(9, 6), layout="constrained")
            ax.plot(d[:, ix], d[:, iy], color="blue")
            ax.set_xlabel(etiquetas[ix], rotation=0, labelpad=12)
            ax.set_ylabel(etiquetas[iy], rotation=90, labelpad=12)
            ax.grid(True, linewidth=0.6, alpha=0.6)
            if self.titulo: fig.suptitle(self.titulo)
            if ahora: plt.show(block=self.bloqueo)
            return

        if self.superponer_canales:
            canales_bounds = []
            offset = 0
            for w in self.anchos:
                canales_bounds.append((offset, offset + w))
                offset += w
            n_canales = len(canales_bounds)
            fig, ejes = plt.subplots(n_canales, 1, sharex=True, figsize=(10, 2.6 * n_canales), constrained_layout=True)
            if n_canales == 1: ejes = [ejes]

            for k_canal, (ini, fin) in enumerate(canales_bounds):
                ax = ejes[k_canal]
                label_canal = self.guiones[k_canal] if self.guiones and len(self.guiones) == len(self.anchos) else (etiquetas[ini] if etiquetas[ini] else f"Canal {k_canal}")
                for j, idx_senal in enumerate(range(ini, fin)):
                    ax.plot(t, d[:, idx_senal], color=colores[(k_canal * 10 + j) % len(colores)], label=etiquetas[idx_senal] if self.guiones else None)
                ax.set_ylabel(label_canal, rotation=90, labelpad=10)
                ax.grid(True, linewidth=0.5, alpha=0.5)
                ax.set_xlim(t.min(), t.max())
                if fin - ini > 1: ax.legend(loc='upper right', fontsize=18)
                if k_canal == n_canales - 1:
                    ax.set_xlabel("Time [s]", labelpad=8)
                    ax.tick_params(axis='x', labelsize=20, labelbottom=True)
                else:
                    ax.tick_params(axis='x', labelbottom=False, bottom=False)
            if ahora: plt.show(block=self.bloqueo)
            return

        fig, ejes = plt.subplots(n, 1, sharex=True, figsize=(10, 2.6 * n), constrained_layout=True)
        if n == 1: ejes = [ejes]

        for k, ax in enumerate(ejes[:n]):
            ax.plot(t, d[:, k], color=colores[k % len(colores)])
            ax.set_ylabel(etiquetas[k], rotation=90, labelpad=10)
            ax.grid(True, linewidth=0.5, alpha=0.5)
            ax.set_xlim(t.min(), t.max())
            if k == n - 1:
                ax.set_xlabel("Time [s]", labelpad=8)
                ax.tick_params(axis='x', labelsize=20, labelbottom=True)
            else:
                ax.tick_params(axis='x', labelbottom=False, bottom=False)
        if ahora: plt.show(block=self.bloqueo)

class ScopeTiempoReal(Scope):
    tiempo_real = True
    def __init__(self, nombre: str, *senales, max_canales: int = 1,
                 anchos: Optional[Sequence[int]] = None, mostrar: bool = True,
                 guiones: Optional[Sequence[str]] = None,
                 titulo: Optional[str] = None, bloqueo: bool = True,
                 xy_mode: Optional[Tuple[int, int]] = None,
                 superponer_canales: bool = True,
                 max_puntos: int = 10000,
                 cuadricula: Optional[Tuple[int, int]] = None,
                 esperar: bool = False,
                 ventana_tiempo: Optional[float] = None):
        super().__init__(nombre, *senales, max_canales=max_canales, anchos=anchos,
                         mostrar=mostrar, guiones=guiones, titulo=titulo,
                         bloqueo=bloqueo, xy_mode=xy_mode,
                         superponer_canales=superponer_canales,
                         max_puntos=max_puntos, cuadricula=cuadricula)
        self._esperar = bool(esperar)
        self.ventana_tiempo = ventana_tiempo
        self._t_buf = []
        self._y_buf = []
    def actualizar(self, t, y):
        try:
            self._t_buf.append(np.asarray(t))
            self._y_buf.append(np.asarray(y))
        except Exception:
            pass
        if self.mostrar and not self._esperar:
            try:
                import matplotlib.pyplot as plt
                plt.pause(0.001)
            except Exception:
                pass
    def esperar(self):
        if self._esperar and self.mostrar:
            try:
                import matplotlib.pyplot as plt
                plt.show(block=True)
            except Exception:
                pass
