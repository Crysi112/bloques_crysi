# bloques_crysi — Biblioteca de Simulación Electromecánica y de Potencia de Tiempo Continuo

**Versión:** 0.3.0 · **Núcleo:** C compilado (GCC/Clang/MSVC) con interfaz Python vía *ctypes* y backend de referencia en NumPy · **Licencia:** MIT  
**Autoría:** Laboratorio CRYSI — Sistemas Ciberfísicos y Transición Energética · **Fecha:** septiembre de 2026

![Versión](https://img.shields.io/badge/versión-0.3.0-blue)
![Licencia](https://img.shields.io/badge/licencia-MIT-green)
![Python](https://img.shields.io/badge/python-%3E%3D3.9-yellow)
![Tests](https://img.shields.io/badge/tests_núcleo-259_passed-brightgreen)

> Simulación de transitorios electromecánicos con núcleo en C y API causal en Python: máquinas eléctricas, electrónica de potencia, baterías y co-simulación con redes — con paridad numérica verificada C↔NumPy.

---

## Resumen

La biblioteca `bloques_crysi` constituye un entorno integral de modelado, simulación y análisis de sistemas dinámicos híbridos, concebido para la investigación y la docencia avanzada en máquinas eléctricas, electrónica de potencia, almacenamiento electroquímico y dinámica vehicular. Su arquitectura adopta el paradigma de diagramas de bloques de flujo causal con un núcleo numérico de alto rendimiento en lenguaje C y una capa de especificación, orquestación y postprocesamiento en Python.

El sistema resuelve de forma unificada ecuaciones diferenciales ordinarias no lineales, sistemas algebraico-diferenciales con lazos algebraicos implícitos y fenómenos de conmutación y eventos de cruce por cero, garantizando paridad numérica entre la implementación nativa en C y el backend homólogo en NumPy/SciPy ($\text{rtol}=10^{-8}$, $\text{atol}=10^{-9}$). Integra, bajo una formalidad común, el modelado de máquinas de inducción, síncronas y de corriente continua, la conversión estática DC-DC y DC-AC, los circuitos equivalentes de baterías y la co-simulación EMT–fasorial con flujos de potencia.

**Palabras clave:** Simulación EMT; Control Orientado a Campo (FOC); Máquinas Eléctricas; Electrónica de Potencia; Sistemas de Almacenamiento; Integración Numérica; Análisis Nodal Modificado; Co-simulación.

---

## Instalación

**Requisitos:** Python ≥ 3.9, un compilador C (GCC vía MSYS2 UCRT64 en Windows, `gcc` o Clang en Linux/macOS), `numpy` y `scipy`. `matplotlib` es opcional (solo para gráficos); `pandapower` y `opendssdirect.py` son opcionales (solo para co-simulación con redes).

```bash
git clone https://github.com/Crysi112/bloques_crysi.git
cd bloques_crysi
pip install -e ".[test]"
```

La biblioteca dinámica del núcleo C (`core/bloques_core.dll` o `.so`) **se compila sola** en la primera simulación. Para forzar la compilación manualmente:

```bash
python -c "from bloques_crysi._clib import compilar; compilar(fuerza=True)"
```

---

## Ejemplo rápido en 30 segundos

Lazo cerrado de primer orden ($d=2$, planta integradora): la salida converge a la referencia.

```python
from bloques_crysi import Modelo, FuenteConstante, Suma, Ganancia, Integrador

m = Modelo(dt=1e-4, metodo="rk4")
ref = m.add(FuenteConstante("ref", 1.0))
s   = m.add(Suma("err", (1.0, -1.0)))
k   = m.add(Ganancia("Kp", 2.0))
p   = m.add(Integrador("planta"))

m.conectar(ref.salida, s.entrada[0:1])
m.conectar(p.salida,   s.entrada[1:2])
m.conectar(s.salida,   k.entrada)
m.conectar(k.salida,   p.entrada)

res = m.run(2.0, registrar=[p])
print(res.get("planta")[-1])   # 0.9817 ≈ 1 - e^{-4} (teórico)
```

Más casos (FOC de PMSM, micro-red IEEE-13, HIL, MNA) en [`ejemplos/`](#12-estructura-del-repositorio-y-ejemplos).

---

## Índice

- [Instalación](#instalación)
- [Ejemplo rápido en 30 segundos](#ejemplo-rápido-en-30-segundos)
1. [Fundamentación y Arquitectura](#1-fundamentación-y-arquitectura)
2. [Catálogo Exhaustivo de Bloques y Ecuaciones](#2-catálogo-exhaustivo-de-bloques-y-ecuaciones)
3. [Máquinas Eléctricas](#3-máquinas-eléctricas)
4. [Electrónica de Potencia](#4-electrónica-de-potencia)
5. [Almacenamiento, Térmica y Vehículo](#5-almacenamiento-térmica-y-vehículo)
6. [Instrumentación: Scopes](#6-instrumentación-scopes)
7. [Sistema de Puertos y Sensores](#7-sistema-de-puertos-y-sensores)
8. [Motor de Simulación: Modelo y Resultado](#8-motor-de-simulación-modelo-y-resultado)
9. [Análisis Nodal Modificado (MNA)](#9-análisis-nodal-modificado-mna)
10. [Co-simulación con Redes Eléctricas](#10-co-simulación-con-redes-eléctricas)
11. [Generación de Código Embebible](#11-generación-de-código-embebible)
12. [Estructura del Repositorio y Ejemplos](#12-estructura-del-repositorio-y-ejemplos)
13. [Validación](#13-validación)
14. [Referencias](#14-referencias)

---

## 1. Fundamentación y Arquitectura

### 1.1 Principio Causal

Todo sistema se describe como grafo dirigido de bloques. Cada bloque implementa una relación estática $\mathbf{y}=\mathbf{g}(\mathbf{u})$ o dinámica $\dot{\mathbf{x}}=\mathbf{f}(\mathbf{x},\mathbf{u}),\quad\mathbf{y}=\mathbf{h}(\mathbf{x},\mathbf{u})$. Las señales fluyen por puertos tipificados (`Puerto`), con inferencia automática de índices de señal (`in_idx`, `out_idx`) y detección de lazos algebraicos.

### 1.2 Integración Temporal

El núcleo resuelve el problema de valor inicial $\dot{\mathbf{x}}=\mathbf{f}(\mathbf{x},\mathbf{u},t)$ mediante:

*Euler explícito:* $\mathbf{x}_{k+1}=\mathbf{x}_k+h\mathbf{f}(\mathbf{x}_k,t_k)$, error $\mathcal{O}(h^2)$ local, $\mathcal{O}(h)$ global.

*Runge-Kutta 4:* $\mathbf{x}_{k+1}=\mathbf{x}_k+\tfrac{h}{6}(\mathbf{k}_1+2\mathbf{k}_2+2\mathbf{k}_3+\mathbf{k}_4)$, error $\mathcal{O}(h^5)$ local, $\mathcal{O}(h^4)$ global.

Selección vía `Modelo(dt, metodo="euler"|"rk4")`.

### 1.3 Lazos Algebraicos

Interconexiones cíclicas entre bloques estáticos generan $\mathbf{y}=\mathbf{g}(\mathbf{y})$. Resolución por Gauss-Seidel con relajación SOR:

```math
\mathbf{y}^{(k+1)}=\mathbf{y}^{(k)}+\omega(\mathbf{g}(\mathbf{y}^{(k)})-\mathbf{y}^{(k)}),\quad \|\mathbf{y}^{(k+1)}-\mathbf{y}^{(k)}\|_\infty<\varepsilon
```

con $\omega=w_{opt}$, $\varepsilon=tol$, $k_{\max}$ = `max_iter`. No convergencia → `RuntimeError` con diagnóstico.

### 1.4 Estratos

| Estrato | Módulos | Función |
|---|---|---|
| Especificación | `bloques.py`, `maquinas.py`, `potencia.py`, `pwm.py`, `baterias.py`, `solar.py`, `mna.py` | Definición declarativa de topologías |
| Resolución | `core/block_core.c/h` | Núcleo C: integración, lazo algebraico, eventos |
| Orquestación | `modelo.py`, `scope.py`, `backend_numpy.py`, `codegen.py`, `red/` | Grafo, despacho, visualización, co-simulación |

---

## 2. Catálogo Exhaustivo de Bloques y Ecuaciones

Todos los bloques derivan de `Bloque(nombre, Ts)` con atributos `op`, `n_in`, `n_out`, `n_state`, `param`, `estados_iniciales`, `entrada`/`salida` (`Puerto`). A continuación se documenta la totalidad del catálogo disponible, con su formulación matemática canónica.

### 2.1 Fuentes de Excitación

#### FuenteConstante

```math
y(t)=P_0
```

Parámetro: `valor`.

#### FuenteEscalon

```math
y(t)=\begin{cases} P_2 & t\lt P_1 \\ P_0 & t\ge P_1 \end{cases}
```

Parámetros: `valor_final`, `t_paso`, `valor_inicial`.

#### FuenteRampa

```math
y(t)=\begin{cases} P_2 & t\lt P_1 \\ P_2+P_0(t-P_1) & t\ge P_1 \end{cases}
```

Parámetros: `pendiente`, `t_inicio`, `offset`.

#### FuenteSeno

```math
y(t)=P_0\sin(2\pi P_1 t+P_2)+P_3
```

Parámetros: `amplitud`, `frecuencia`, `fase`, `offset`.

#### FuenteTrifasica

```math
\begin{aligned} y_0&=P_0\sin(2\pi P_1 t+P_2)\\ y_1&=P_0\sin(2\pi P_1 t+P_2-2\pi/3)\\ y_2&=P_0\sin(2\pi P_1 t+P_2-4\pi/3) \end{aligned}
```

$n_{out}=3$.

#### PulsoRectangular

```math
y(t)=P_4+\begin{cases} P_0 & \mathrm{mod}(t+P_3,P_1)\lt P_2P_1 \\ 0 & \text{en otro caso} \end{cases}
```

Parámetros: `amplitud`, `periodo`, `duty`, `fase`, `offset`. Validación: $periodo>0$, $0\lt duty\le1$.

#### FuenteCSV
Interpolación lineal (o retención) sobre tabla $(t_i,y_i)$ leída de archivo CSV con delimitador autodetectado `,`/`;` y codificación `utf-8-sig`:

```math
y(t)=y_k+\frac{y_{k+1}-y_k}{t_{k+1}-t_k}(t-t_k),\quad t_k\le t\lt t_{k+1}
```

Parámetro: `archivo`, `columna_t`, `columna_y`, `interpolar`.

#### FuenteTabla
Idéntica a `FuenteCSV` pero con puntos en memoria $[(t_0,y_0),\dots,(t_{n-1},y_{n-1})]$ estrictamente crecientes en $t$.

### 2.2 Operadores Estáticos Escalares y Vectoriales

#### Ganancia

```math
y = P_0\cdot u
```

$n_{in}=1$.

#### Suma

```math
y = \sum_{i=0}^{n-1} P_i\,u_i
```

$n_{in}=|signos|$, $P_i$ son los signos/pesos. Puertos: `entrada` vectorial.

#### Saturar

```math
y = \min(\max(u,P_0),P_1)
```

Parámetros: `u_min`, `u_max`.

#### Multiplicador

```math
y = u_0\cdot u_1
```

$n_{in}=2$, $n_{out}=1$, canales `["y"]`.

#### SaturarVectorial
Dado $\mathbf{v}^\ast=(v_d^\ast,v_q^\ast)$ y $V_{\max}=P_0$:

```math
\mathbf{v}=\begin{cases} \mathbf{v}^* & \|\mathbf{v}^*\|\le V_{\max} \\ V_{\max}\mathbf{v}^*/\|\mathbf{v}^*\| & \text{en otro caso} \end{cases}
```

Preserva fase, escala magnitud. Fundamental en FOC para respetar hexágono SVPWM ($V_{\max}=V_{dc}/\sqrt3$).

#### Transformaciones de Referencia

**Clarke** ($3\to2$, amplitud invariante):

```math
\begin{bmatrix}\alpha\\\beta\end{bmatrix}=\frac23\begin{bmatrix}1&-1/2&-1/2\\0&\sqrt3/2&-\sqrt3/2\end{bmatrix}\begin{bmatrix}a\\b\\c\end{bmatrix}
```

**InvClarke:**

```math
\begin{aligned} a&=\alpha\\ b&=-0.5\alpha+0.866\beta\\ c&=-0.5\alpha-0.866\beta \end{aligned}
```

**Park** ($\alpha\beta\to dq$, ángulo $\theta$):

```math
\begin{bmatrix}d\\q\end{bmatrix}=\begin{bmatrix}\cos\theta&\sin\theta\\-\sin\theta&\cos\theta\end{bmatrix}\begin{bmatrix}\alpha\\\beta\end{bmatrix}
```

**InvPark:**

```math
\begin{bmatrix}\alpha\\\beta\end{bmatrix}=\begin{bmatrix}\cos\theta&-\sin\theta\\\sin\theta&\cos\theta\end{bmatrix}\begin{bmatrix}d\\q\end{bmatrix}
```

**TransformadaQD:** Bloque compuesto $n_{in}=7$ ($v_{abc},i_{abc},\theta$), $n_{out}=4$ ($v_{qs},v_{ds},i_{qs},i_{ds}$) que encapsula Clarke+Park para medición FOC.

#### Multiplexor / Demultiplexor
Concatenación y separación de buses vectoriales sin dinámica:

```math
\mathbf{y}=[u_0,\dots,u_{n-1}]^T,\quad y_k=u_k
```

#### Tabla 1D / 2D / 3D (LUT)
Interpolación multilineal sobre mallas estrictamente crecientes.

1D: Dado $x$, con $x_k\le x\lt x_{k+1}$:

```math
y=y_k+\frac{y_{k+1}-y_k}{x_{k+1}-x_k}(x-x_k)
```

2D: Bilineal sobre $(x_i,y_j,z_{ij})$.

3D: Trilineal sobre $(x_i,y_j,z_k,w_{ijk})$.

Parámetros almacenan $n_x$, $x_i$, $y_i$, $z_{ij\dots}$.

#### Logico y Relacional

**Logico:** $y=\text{OP}(u_0,\dots,u_{n-1}>P_1)$, OP $\in\lbrace\text{AND},\text{OR},\text{NAND},\text{NOR},\text{XOR},\text{XNOR},\text{NOT}\rbrace$, $n_{in}$ configurable, umbral $P_1$.

**Relacional:** $y=1$ si $u_0\text{ OP }u_1 \pm tol$, OP $\in\lbrace=,\neq,<,\le,>,\ge\rbrace$.

#### LimitadorRapidez (Rate Limiter)
Estado $x$, $h=dt$:

```math
x_{k+1}=x_k+\mathrm{clip}(u_k-x_k,\,-P_1h,\,P_0h)
```

Parámetros: `subida`, `bajada` [unidades/s].

#### RetenedorDisparado (Triggered Hold)
Estados $x_0$ (valor retenido), $x_1$ (memoria de disparo). Si $u_{trig}>P_0$ y flanco:

```math
x_0\gets u_{sig}
```

#### MaquinaEstados
Autómata finito determinista con $n_{estados}$, $n_{entradas}$ y tabla de transiciones $(desde,hacia,idx_{señal},cond,umbral)$, cond $\in\lbrace<,\le,>,\ge,=,\neq\rbrace$.

```math
s_{k+1}=\begin{cases} hacia & s_k=desde \land (u_{idx}\,cond\,umbral)\\ s_k & \text{en otro caso} \end{cases}
```

#### Filtros como FuncionTransferencia

**FuncionTransferencia:** Discretización Tustin (bilineal) $s\to\frac{2}{h}\frac{z-1}{z+1}$ vía `scipy.signal.bilinear`, orden $n=\max(|num|,|den|)-1$, $n_{state}=2n$.

**FiltroPasoBajo:** $H(s)=\frac{\omega_c}{s+\omega_c}$ (orden1) o $\frac{\omega_c^2}{s^2+2\zeta\omega_cs+\omega_c^2}$ (orden2), $\omega_c=2\pi f_c$.

**FiltroPasoAlto:** Orden1 $\frac{s}{s+\omega_c}$, orden2 $\frac{s^2}{s^2+2\zeta\omega_cs+\omega_c^2}$.

**FiltroNotch:** $H(s)=\frac{s^2+\omega_0^2}{s^2+2\zeta\omega_0s+\omega_0^2}$, $\omega_0=2\pi f_n$.

### 2.3 Bloques Dinámicos de Control y Sincronización

#### Integrador

```math
\dot x = u,\quad y=x,\quad x(0)=x_0
```

#### PID con Anti-Windup y Derivada Filtrada
Estados: $x_0=\int e$, $x_1=e_{prev}$, $x_2=u_{d,filt}$.

```math
\begin{aligned} u_d &= \frac{K_d(e_k-x_1)+T_f x_2}{T_f+h}\\ u_{raw}&=K_pe_k+K_ix_0+u_d\\ u&=\mathrm{clip}(u_{raw},u_{\min},u_{\max})\\ \dot x_0&= \begin{cases}0 & (u_{raw}>u_{\max}\land e_k>0)\lor(u_{raw}\lt u_{\min}\land e_k\lt 0)\\ e_k & \text{en otro caso}\end{cases} \end{aligned}
```

Parámetros: $K_p,K_i,K_d,T_f,u_{\min},u_{\max}$.

#### PLLTrifasico
Estados $\theta$, $w_{int}$. Entradas $v_a,v_b,v_c$, salidas $\omega,\theta$:

```math
\begin{aligned} (\alpha,\beta)&=\text{Clarke}(v_{abc})\\ (v_d,v_q)&=\text{Park}(\alpha,\beta,\theta)\\ e&=v_q\\ w_{int}&\gets w_{int}+hK_i e\\ \omega&=2\pi f_{ff}+K_pe+w_{int}\\ \theta&\gets\theta+h\omega \end{aligned}
```

#### EjeMecanico
Múltiples máquinas acopladas rígidamente. Estados $\omega_m,\theta_m$. Entradas $T_{e,i}$ ($n_{maq}$) y $T_L$:

```math
J_{eq}\dot\omega_m=\sum_i T_{e,i}-T_L-B_{m,eq}\omega_m,\quad \dot\theta_m=\omega_m
```

#### MasaTermica

```math
C_{th}\dot T = \sum_k P_k - (T-T_{amb})/R_{amb},\quad T(0)=T_{inicial}
```

Parámetros: $C_{th},T_{amb},R_{amb}$, $n_{in}$ configurable.

#### ResistenciaTermica

```math
q = (T_1-T_2)/R
```

#### Engranaje

```math
\omega_2=\omega_1/a,\quad T_2 = a\,T_1
```

#### EjeFlexible
Estados $\theta_1,\theta_2$:

```math
\begin{aligned} T_{eje}&=K(\theta_1-\theta_2)+B(\omega_1-\omega_2)\\ \dot\theta_i&=\omega_i \end{aligned}
```

#### Embrague

```math
T_{out}=\begin{cases} \min(T_{in},T_{\max}) & u_{ctrl}>umbral\\ 0 & \text{en otro caso} \end{cases}
```

### 2.4 Bloques de Fallas

**FalloProgramado:**

```math
y=\begin{cases}u & t\lt t_f\\ valor & t\ge t_f\land modo=0\\ u+valor & modo=1\end{cases}
```

**FalloEvento:** Idem pero disparado por $u_{trig}>umbral$.

---

## 3. Máquinas Eléctricas

Todas derivan de `Maquina` con $n_{out}=10$ (o específico), puertos `terminales` (3), `T_L`, y sensores: `sensor3V/3I`, `sensorVelocidad` ($w_m$), `sensorPosicion` ($\theta_{rm}$), `sensorPosicionElectrica` ($\theta_e$), `sensorPar` ($T_e$), `sensorCorrienteD/Q` ($i_{ds},i_{qs}$), `sensorPerdidasEstator` ($P_{cu,s}$), `sensorCorrienteRotor` ($i'\sb{ar},i'\sb{br},i'\sb{cr}$), `resumen()`.

### 3.1 MaquinaImanesPermanentes (PMSM/PMAC) — `OP_MAQ_PMAC`

Estados $[i_{qs},i_{ds},\omega_m,\theta_e]$, entradas $[v_a,v_b,v_c,T_L]$ (o 5 con puerto mecánico externo si `mecanica_interna=False`).

Ecuaciones $dq$ síncronas:

```math
\begin{aligned} v_d &= R_s i_d + L_d\dot i_d - \omega_e L_q i_q \\ v_q &= R_s i_q + L_q\dot i_q + \omega_e(L_d i_d+\lambda_m) \\ T_e &= \tfrac32\tfrac{P}{2}[\lambda_m i_q+(L_d-L_q)i_d i_q] \\ J\dot\omega_m &= T_e-T_L-B_m\omega_m,\quad \dot\theta_e=\tfrac{P}{2}\omega_m \end{aligned}
```

Parámetros: $R_s,L_d,L_q,\lambda_m,P,J,B_m$, $\theta_0$, `saturacion` opcional como LUT $i_d\to L_d(i_d)$.

### 3.2 MaquinaInduccion — `OP_MAQ_INDUCCION`

Estados 6: $[\lambda\sb{qs},\lambda\sb{ds},\lambda'\sb{qr},\lambda'\sb{dr},\omega\sb{m},\theta\sb{m}]$, salidas 13 incluyendo $i'\sb{ar},i'\sb{br},i'\sb{cr}$.

Modelo $\alpha\beta$ estacionario (Krause):

```math
\begin{aligned} \dot{\boldsymbol{\lambda}}_s &= \mathbf{v}_s-R_s\mathbf{i}_s \\ \dot{\boldsymbol{\lambda}}'_r &= -R'_r\mathbf{i}'_r+j\omega_r\boldsymbol{\lambda}'_r \\ \mathbf{i}_s &= L_{i00}\boldsymbol{\lambda}_s+L_{i01}\boldsymbol{\lambda}'_r \\ \mathbf{i}'_r &= L_{i01}\boldsymbol{\lambda}_s+L_{i11}\boldsymbol{\lambda}'_r \end{aligned}
```

con $L_{i}=L^{-1}$, $L_s=L_{ls}+L_m$, $L_r=L_{lr}+L_m$. Par $T_e=\tfrac32\tfrac{P}{2}\Im\lbrace\boldsymbol{\lambda}_s^*\mathbf{i}_s\rbrace$. `velocidad_sincronica`= $4\pi f/P$ con $f=60$ Hz por defecto.

### 3.3 MaquinaSincrona — `OP_MAQ_SINCRONA`

Estados 8: flujos $q/d$ (estator, campo, amortiguadores) + $\omega,\theta$. Entradas 5: $v_{abc},v_{fd},T_L$. Parámetros $R_s,R_{fd},R_{kq1},R_{kq2},R_{kd},L_{ls},L_{mq},L_{lkq1},L_{lkq2},L_{md},L_{lf},L_{lkd},P,J,B_m$. Matrices $3\times3$ $L_q,L_d$ invertidas a $L_{iq},L_{id}$ para derivadas. Sensor `sensorVoltajeCampo`.

### 3.4 MaquinaCorrienteContinua — `OP_MAQ_CC`

Estados $[i_a,i_f,\omega_m,\theta_m]$, salidas 8: $i_a,i_f,\omega_m,\theta_m,T_e,E_a,V_t,P_{cu}$.

```math
\begin{aligned} L_a\dot i_a &= V_t - R_a i_a - L_{AF}i_f\omega_m \\ L_f\dot i_f &= V_f - R_f i_f \\ T_e &= L_{AF}i_f i_a \\ J\dot\omega_m &= T_e-T_L-B_m\omega_m \end{aligned}
```

### 3.5 MaquinaDCImanesPermanentes — `OP_MAQ_DC_PM`

Estados $[i_a,\omega_m,\theta_m]$, salidas 7.

```math
\begin{aligned} L_a\dot i_a &= V_a - R_a i_a - K_t\omega_m \\ T_e &= K_t i_a \\ J\dot\omega_m &= T_e-T_L-B_m\omega_m \end{aligned}
```

Sensores `sensorCorriente`, `sensorVelocidad`, `sensorPar`, `sensorEa` ($K_t\omega_m$), `sensorPerdidas`.

---

## 4. Electrónica de Potencia

### 4.1 Convertidores DC-DC (Monolíticos)

Evitan el lazo algebraico $V\to I$ de interconectar inductor e interruptor causalmente. Estados $[i_L,v_C]$, entradas $[V_{in},d]$, salidas $[v_{out},i_L]$.

**Buck:**

```math
L\dot i_L = dV_{in}-v_C,\quad C\dot v_C = i_L-v_C/R
```

**Boost:**

```math
L\dot i_L = V_{in}-(1-d)v_C,\quad C\dot v_C = (1-d)i_L-v_C/R
```

**BuckBoost:**

```math
L\dot i_L = dV_{in}-(1-d)v_C,\quad C\dot v_C = (1-d)i_L-v_C/R
```

Parámetros $L,C,R$.

### 4.2 RectificadorTrifasico

Estado $v_{dc}$, entradas $v_{abc}$, salidas $[v_{dc},i_{dc}]$.

```math
v_{rec}= \max(v_{abc})-\min(v_{abc}),\quad i_{ch}=(v_{rec}-v_{dc})/R_{int}\,( \ge0),\quad C\dot v_{dc}=i_{ch}-v_{dc}/R
```

### 4.3 InversorTrifasico / Monofasico

Modelos promediado vs. conmutado seleccionables por `conmutada`. Estados 6 (3F) /2 (1F): tensiones y corrientes de filtro $LC$.

Entradas $m$ (índice), salidas $[v_{Ca},v_{Cb},v_{Cc},i_{La},i_{Lb},i_{Lc}]$. Parámetros $f_{out},f_{sw},m_{start},m_{end},t_{ramp},L_f,C_f,R$.

Sensores `sensorVoltajesSalida`, `sensorCorrientesFase`.

### 4.4 Cargas y Elementos Pasivos

**CargaRLTrifasica:** Estados $[i_a,i_b]$, $i_c=-(i_a+i_b)$.

```math
\dot i_k = (v_k-v_n - R i_k)/L,\quad v_n=\tfrac13\sum v_k
```

Constructor alternativo `desde_pq(p_w,q_var,v_ll,f)` calcula $R=3V_{ln}^2P/(P^2+Q^2)$, $L=3V_{ln}^2Q/(P^2+Q^2)/(2\pi f)$.

**CargaPQTrifasica/Monofasica:** Modelo algebraico fasorial $i = (P-jQ)/V^*$ proyectado a $\alpha\beta$.

**Transformador:** Ideal $v_2=v_1/a$, $i_1=i_2/a$ (monofásico) o $6\times6$ para trifásico.

**Resistencia:** $i=v/R$.

**Inductor (MNA):** $L\dot i = v$, estado $i$.

**Capacitor (MNA):** $C\dot v = i$, estado $v$.

---

## 5. Almacenamiento, Térmica y Vehículo

### 5.1 BateriaECM — `OP_BATERIA_ECM`

Estados $[SOC, v_1, v_2, T_{cell}]$, entradas $[i_{load}]$, salidas 6: $[V_{term},SOC,T_{cell},P_{loss},I_{chg,lim},I_{dch,lim}]$.

Ecuaciones:

```math
\begin{aligned} \dot{SOC} &= -i_{cell}/(Q_{nom}3600)\\ \dot v_1 &= (i_{cell}-v_1/R_1)/C_1\\ \dot v_2 &= (i_{cell}-v_2/R_2)/C_2\\ C_{th}\dot T &= P_{loss}-(T-T_{amb})/R_{th}\\ V_{oc}&= f_{LUT}(SOC)\\ V_{cell}&=V_{oc}-i_{cell}R_0-v_1-v_2\\ P_{loss,cell}&=i_{cell}^2R_0+v_1^2/R_1+v_2^2/R_2 \end{aligned}
```

Parámetros: $Q_{nom},V_{nom},R_0,R_1,C_1,R_2,C_2,N_s,N_p$, tabla OCV-SOC (14 puntos NMC por defecto), límites $I_{chg},I_{dch},T_{min/max}$, $R_{th},C_{th},T_{amb}$, degradación. Escalado pack $V_{term}=N_sV_{cell}$, $i_{cell}=i_{load}/N_p$.

### 5.2 Panel Solar, Vehiculo, CalculoIdc

*Solar* y *Vehículo* siguen modelos de primer principio (irradiancia–temperatura y dinámica longitudinal $m_{eq}\dot V = F_t - \tfrac12\rho C_dAV|V|-C_{rr}mg - mg\sin\alpha$). `CalculoIdc` estima corriente DC del bus: $I_{dc}= (P_{ac}/\eta)/V_{dc}$.

---

## 6. Instrumentación: Scopes

### 6.1 Scope

Bloque sumidero $n_{in}=\sum anchos$, $n_{out}=0$, `es_scope=True`. No participa en el lazo algebraico. Parámetros: `anchos`, `max_canales`, `mostrar`, `guiones`, `titulo`, `bloqueo`, `xy_mode`, `superponer_canales`, `max_puntos`, `cuadricula`.

Métodos:

*`datos(res)`*: retorna $(t,d)$ con $d\in\mathbb{R}^{N\times n}$.

*`mostrar_grafico(res)`*: Renderiza con `matplotlib` (estilo Times New Roman/Stix, $\ge 22$ pt), con diezmado adaptativo si $N>$ `max_puntos`:

- `decimar_datos(t,y,metodo="step"|"lttb"|"minmax")`: *step* $t[::k]$, *LTTB* preserva triángulos, *minmax* conserva envolvente por bins ($y_{\min},y_{\max}$ por bin). Complejidad $\mathcal{O}(N)$.

Modos: XY (`xy_mode=(ix,iy)`), superposición por canal (`superponer_canales=True` genera $n_{canales}$ subplots compartidos en X), o panel por señal.

Atajos: `m.scope(nombre, *señales, anchos=[...])` crea, conecta y registra en una línea; `formalizar_etiqueta` mapea $i_a\to i_{as}$ [A], $T_e$ [N·m], etc.

### 6.2 ScopeTiempoReal

Subclase con `tiempo_real=True`, `esperar`, `ventana_tiempo`. Métodos `actualizar(t,y)` (acumula en `_t_buf/_y_buf` y `plt.pause(0.001)` si `mostrar` y no `esperar`) y `esperar()` (`plt.show(block=True)`). Integración en `Modelo.run` vía `iterar` por *chunks* ($chunk=\max(10h, t_{fin}/50)$) y `concatenar_resultados`.

---

## 7. Sistema de Puertos y Sensores

**Puerto** (`puertos.py`): descriptor tipificado `(bloque, tipo, offset, n, canales)`. Tipos `ent`/`sal`. Atributos `n`, `puerto`. Conexión `m.conectar(src, dst)` resuelve `in_idx/out_idx` y valida rangos. Soporta *slicing* `par(bloque,"ent",off,n)` y buses `Multiplexor`.

**Sensor** (subclase de `Puerto` virtual): vista sin costo sobre salidas de `Maquina` u otros bloques. Ejemplos: `maq.sensorVelocidad()` → `Sensor("wm", maq, "sal",5,1)`, `maq.sensor3I()` → 3 canales. No añade ecuaciones, solo alias de índices.

---

## 8. Motor de Simulación: Modelo y Resultado

### 8.1 Modelo

```python
m = Modelo(dt=1e-4, metodo="rk4", max_iter=50, tol=1e-9, w_opt=1.0)
b = m.add(Bloque(...))
m.conectar(src.salida, dst.entrada)
m.scope("nombre", sig1, sig2, anchos=[1,2])
res = m.run(t_fin, registrar=[sensor, bloque])
```

Métodos:

*`add`*: registra bloque, asigna `in_idx/out_idx` provisionales.
*`conectar`*: cableado causal, resuelve `_resolver()` (orden topológico, detección de feedthrough, asignación de `n_sig`, `n_alg`, `alg_list`).
*`_resolver`*: análisis de grafo, es usado por `codegen` y `backend_numpy`.
*`run`*, `iterar`, `paso`, `iniciar`*: despacho a núcleo C (`_clib.BloqueC/ModeloC`) o fallback NumPy. Manejo de eventos de cruce por cero con bisección y `Modelo.set_param(bloque,param,valor)` seguro.
*`acoplar_red`, `acoplar_pcc`*: atajos a `red.CoSimuladorRed`.
*`indice_tiempo(t)`*, `scope()`*, `resumen()`*.

### 8.2 Resultado

Subclase de `dict` con `t` (vector tiempo) y matrices por clave. Helpers:

*`get(key)` / `get_2d(key)`*: acceso con *squeeze* automático.
*`guardar_csv(ruta, separador)`*: exportación con cabecera `t,señal`.
*`pico(key, t_max)`*: $\max|y|$.
*`final(key, ventana)`*: media en ventana final.
*`tiempo_establecimiento(key, tolerancia)`*: último instante fuera de banda $\pm tol\cdot|y_{final}|$.
*`thd(key,f0,t_inicio,n_armonicos)`*, `fft(key)`*: análisis espectral para THD y espectro unilateral.

---

## 9. Análisis Nodal Modificado (MNA)

Bloque `SubredMNA` (`mna.py`) encapsula subredes eléctricas no causales. Primitivas: `Nodo`, `Resistor(R)`, `Capacitor(C)`, `Inductor(L)`, `FuenteTension`, `FuenteCorriente`, `Interruptor(R_{on},R_{off},control)`, `DiodoIdeal(R_{on},R_{off},V_f)`, `VCVS`, `VCCS`, `MutualInductor`.

Estampación Dommel: $G\mathbf{x}+C\dot{\mathbf{x}}=\mathbf{b}$, con $C_{eq}=C/h$ (Euler hacia atrás, método fijado por el núcleo). Sistema $A\mathbf{x}\sb{k+1}=\mathbf{b}\sb{eq}$ resuelto por LU densa $128\times128$ con pivoteo parcial. Diodos: iteración de complementariedad $v_d\ge V_f \perp i_d\ge0$ con $G_{on}=1/R_{on}$. Exposición al exterior vía `n_out` tensiones diferenciales e corrientes de rama.

---

## 10. Co-simulación con Redes Eléctricas

Módulo `red/` desacoplado:

*`BackendRed`* (ABC): `set_carga(bus,P,Q)`, `set_generacion`, `runpp()`, `get_tension(bus)->TensionBus(magnitud_v,angulo_rad)`, `get_corriente_linea`. Unidades SI en frontera, conversiones en `unidades.py` (`SistemaUnidades`).

*`CoSimuladorRed`*: orquestación por ventanas $H=dt_{red}$ con promediado $P,Q$, `runpp`, sub-relajación $V^\alpha=\alpha V_{new}+(1-\alpha)V_{old}$, convergencia $\lVert V^\alpha-V_{old}\rVert<\varepsilon$.

Atajos `Modelo.acoplar_red(backend, bus_pcc, elemento, v_nominal_ll, dt_red)` y `acoplar_pcc(maquina,bus_idx)`.

Adaptadores opcionales `adaptador_pandapower`, `adaptador_opendss` (extras `red-pandapower`, `red-opendss`).

---

## 11. Generación de Código Embebible

`codegen.GeneradorCodigo(Modelo).codigo(nombre)` / `.generar(ruta)` emite C autónomo con cabecera, `#include <math.h>`, `static double _sig[]`, `static double _bX_x[]`, `lu_solve_dense`, funciones `_actualizar_fuentes`, `_lazo_algebraico`, `_actualizar_dinamicos` y API `nombre_init/paso/run/t/sig`, compilable con `gcc -O2 -o sim modelo.c -lm` y apto para firmware/HIL.

---

## 12. Estructura del Repositorio y Ejemplos

```
bloques_crysi/
├── src/bloques_crysi/
│   ├── core/            # Núcleo C (block_core.c/h, serial_win.h)
│   ├── red/             # Co-simulación (BackendRed, adaptadores, unidades)
│   ├── bloques.py       # Bloques elementales, control y fallas
│   ├── maquinas.py      # Máquinas eléctricas
│   ├── potencia.py      # Convertidores, cargas y pasivos
│   ├── mna.py           # Análisis Nodal Modificado
│   ├── modelo.py        # Orquestador (Modelo, Resultado)
│   ├── scope.py         # Instrumentación
│   ├── codegen.py       # Generación de C embebible
│   ├── backend_numpy.py # Backend de referencia NumPy
│   └── pwm, solar, baterias, hardware, subsistemas, puertos, opcodes, ...
├── docs/{teoria_modelado.md,diseno_red.md,img/}
├── ejemplos/01_maquinas_electricas ... 07_mna/ + data/  # .py ejecutables
├── tests/               # 314 pruebas (259 núcleo + 55 red)
└── pyproject.toml, LICENSE, build_core.py
```

Ejemplos cubren: FOC PMSM, micro-red IEEE-13, HIL, almacenamiento, protecciones y MNA (buck con diodo).

---

## 13. Validación

Paridad C–NumPy $\text{rtol}=10^{-8},\text{atol}=10^{-9}$, soluciones analíticas (RLC, lazo $y=G r/(1+G)$), benchmarks Krause/Mohan/WLTP. Ejecución del núcleo (sin dependencias de red): `pytest tests/ -k "not red" -v` → `259 passed`. Suite completa (`pip install -e ".[red-pandapower,red-opendss]"`): `pytest tests/ -v` → 314 pruebas.

---

## 14. Referencias

Krause *et al.* (2013) *Analysis of Electric Machinery*; Mohan *et al.* (2003) *Power Electronics*; Park (1929) *AIEE Trans.*; Clarke (1943); Dommel (1969) *IEEE Trans. PAS*.

---

*Fin del catálogo. Para el desarrollo teórico completo véase `docs/teoria_modelado.md`; para la especificación de co-simulación, `docs/diseno_red.md`.*

