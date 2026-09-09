"""Curvas TCC y verificacion de coordinacion (estudios estaticos)."""

from __future__ import annotations
import csv
import math
from typing import Dict, Sequence, Tuple

import numpy as np

INF = float("inf")
K_FUSIBLE_T = 0.35 * (833.0 / 100.0) ** 2
K_FUSIBLE_K = 0.15 * (500.0 / 65.0) ** 2


def curva_51(rele, i):
    v = np.atleast_1d(np.asarray(i, dtype=float))
    t = np.array([rele.calcular_tiempo_trip(float(x)) for x in v])
    return float(t[0]) if np.asarray(i).ndim == 0 else t


def curva_fusible(i, i_rat, K=K_FUSIBLE_T, n=2.0, umbral_pu=1.0):
    v = np.asarray(i, dtype=float)
    return np.where(v > i_rat * umbral_pu, K * (i_rat / v) ** n, INF)


def curva_itm(i, ir=680.0, isd=2720.0, tsd=0.2, ii=6800.0, tr=10.0, ti=0.02):
    v = np.asarray(i, dtype=float)
    t = np.full_like(v, INF)
    t[(v >= ir) & (v < isd)] = tr * (6.0 * ir / v[(v >= ir) & (v < isd)]) ** 2
    t[(v >= isd) & (v < ii)] = tsd
    t[v >= ii] = ti
    return t


def curva_dano(i, i_base, K=1250.0, i_min_pu=2.0):
    v = np.asarray(i, dtype=float)
    return np.where(v >= i_min_pu * i_base, K / (v / i_base) ** 2, INF)


def icc_faultstudy(backend, buses: Sequence[str]):
    icc3: Dict[str, float] = {}
    icc1: Dict[str, float] = {}
    try:
        backend.dss.text("Solve mode=FaultStudy")
        ruta = backend.dss.text("Export Faultstudy").strip()
        with open(ruta, newline="") as f:
            lector = csv.DictReader(f)
            lector.fieldnames = [c.strip() for c in lector.fieldnames]
            for fila in lector:
                fila = {k.strip(): v.strip() for k, v in fila.items()}
                try:
                    icc3[fila["Bus"].upper()] = float(fila["3-Phase"])
                    icc1[fila["Bus"].upper()] = float(fila["1-Phase"])
                except (KeyError, ValueError):
                    pass
    except OSError:
        pass
    for b in buses:
        print(f"  {b}: 3F={icc3.get(b, float('nan')):.1f}  1F={icc1.get(b, float('nan')):.1f}")
    return icc3, icc1


def verificar_cti(filas: Sequence[Tuple[str, float, str, float, float]]) -> bool:
    def fmt(t):
        return f"{t:.3f} s" if np.isfinite(t) else "no arranca"

    ok_all = True
    for dw, t_dw, up, t_up, cti_min in filas:
        cti = t_up - t_dw
        ok = bool(np.isfinite(cti) and cti >= cti_min)
        ok_all = ok and ok_all
        print(f"  [{'OK' if ok else 'REVISAR'}] {dw} ({fmt(t_dw)}) < {up} ({fmt(t_up)}) "
              f"CTI={fmt(cti)} (min {cti_min:.2f} s)")
    return ok_all


def tabla_ajustes(reles) -> None:
    for r in reles:
        print(f"  {r.codigo_ansi:6s} {r.nombre:20s} param={r.param}")


def marcas_icc(icc3: Dict[str, float], colores: Dict[str, str], referral=None):
    referral = referral or {}
    marcas = []
    for b, c in colores.items():
        if b in icc3:
            v = icc3[b] / referral[b] if b in referral else icc3[b]
            suf = " ref MT" if b in referral else ""
            marcas.append((f"Icc3 @ {b} ({v:.0f} A{suf})", v, c))
    return marcas


def figura_tcc(curvas, marcas=(), inrush=None, titulo="Coordinacion TCC",
               subtitulo="", archivo=None, xr=(15.0, 20000.0), yr=(0.01, 1000.0)):
    import plotly.graph_objects as go

    fig = go.Figure()
    for nombre, x, y, color, dash, width, visible in curvas:
        y = np.asarray(y, dtype=float)
        mask = np.isfinite(y) & (y > 0.005)
        fig.add_trace(go.Scatter(
            x=np.asarray(x)[mask], y=y[mask], mode="lines", name=nombre,
            line=dict(color=color, dash=dash, width=width), visible=visible,
            hovertemplate="<b>%{name}</b><br>I: %{x:,.1f} A<br>t: %{y:.3f} s<extra></extra>"))
    for etiqueta, valor, color in marcas:
        fig.add_trace(go.Scatter(
            x=[valor, valor], y=[yr[0], yr[1] * 2.0], mode="lines", name=etiqueta,
            line=dict(color=color, dash="dot", width=1.5),
            hoverinfo="name", visible="legendonly"))
    if inrush is not None:
        nombre, xi, yi = inrush
        fig.add_trace(go.Scatter(
            x=[xi], y=[yi], mode="markers", name=nombre,
            marker=dict(color="red", size=8),
            hovertemplate="I: %{x} A<br>t: %{y} s<extra></extra>", visible="legendonly"))
    fig.update_layout(
        title=f"<b>{titulo}</b><br><sup>{subtitulo}</sup>",
        xaxis_title="<b>Corriente Primaria a 4.16 kV (A)</b>",
        yaxis_title="<b>Tiempo de Operacion (s)</b>",
        xaxis=dict(type="log", range=[math.log10(xr[0]), math.log10(xr[1])], showgrid=True,
                   minor_ticks="inside", minor_showgrid=True, gridcolor="lightgray", minor_gridcolor="whitesmoke"),
        yaxis=dict(type="log", range=[math.log10(yr[0]), math.log10(yr[1])], showgrid=True,
                   minor_ticks="inside", minor_showgrid=True, gridcolor="lightgray", minor_gridcolor="whitesmoke"),
        plot_bgcolor="white", hovermode="x unified",
        legend=dict(title="<b>Protecciones (Clic para aislar)</b>", yanchor="top", y=0.99,
                    xanchor="right", x=0.99, bgcolor="rgba(255, 255, 255, 0.8)",
                    bordercolor="black", borderwidth=1),
        height=850)
    if archivo is not None:
        fig.write_html(str(archivo), auto_open=False)
        print(f"TCC interactiva: {archivo}")
    return fig
