"""Util compartida para tests: correr un modelo por C y por el backend numpy."""

import numpy as np

from bloques_crysi import Modelo
from bloques_crysi.backend_numpy import simular


def resolver_registro(modelo, registrar):
    """Devuelve (rec_idx, grabar) igual que Modelo.run los resuelve."""
    grabar = []
    if not registrar:
        registrar = tuple(modelo.bloques)
    for item in registrar:
        if hasattr(item, "puerto"):  # Sensor
            p = item.puerto
            idx = p.indices()
            grabar.append((item.nombre, item.canales, [int(k) for k in idx]))
        elif hasattr(item, "tipo"):  # Puerto
            idx = item.indices()
            canales = [f"{item.bloque.nombre}[{item.offset + k}]" for k in range(len(idx))]
            grabar.append((item.bloque.nombre, canales, [int(k) for k in idx]))
        elif hasattr(item, "out_idx"):  # bloque
            canales = getattr(item, "NOMBRES", None)
            if canales is None:
                canales = [f"{item.nombre}[{k}]" for k in range(item.n_out)]
            grabar.append((item.nombre, list(canales), [int(k) for k in item.out_idx]))
        else:
            raise TypeError(f"No sé grabar {item!r}.")
    rec_idx = [k for _, _, idxs in grabar for k in idxs]
    return rec_idx, grabar


def correr_c_y_numpy(modelo, t_fin, registrar, retorna_c=None):
    """Ejecuta el modelo por la DLL y por numpy; devuelve (res_c, arr_np, n_steps)."""
    modelo._resolver()
    n_steps = int(round(t_fin / modelo.dt)) + 1

    res_c = modelo.run(t_fin, registrar=registrar)

    rec_idx, _ = resolver_registro(modelo, registrar)
    arr_np = simular(
        modelo.bloques, modelo.dt, t_fin, rec_idx,
        metodo=1 if modelo.metodo == 1 else 0,
        max_iter=modelo.max_iter, tol=modelo.tol, w_opt=modelo.w_opt,
        orden_estatico=modelo._orden_estatico(),
    )
    return res_c, arr_np, n_steps


def flat_c(res_c, registrar, modelo, n_steps):
    """Devuelve un array (n_rec, n_steps) con las señales de res_c en el
    mismo orden plano que usa el backend numpy."""
    _, grabar = resolver_registro(modelo, registrar)
    filas = []
    for nombre, canales, _ in grabar:
        arr = np.asarray(res_c[nombre])
        if arr.ndim == 1:
            arr = arr.reshape(1, -1)
        else:
            arr = arr.T  # (n_steps, n_canales) -> (n_canales, n_steps)
        filas.append(arr)
    return np.vstack(filas)