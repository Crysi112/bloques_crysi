# Fundamentos Teóricos y Modelado Matemático en bloques_crysi

## Tratado Teórico–Pedagógico del Núcleo de Simulación

**Programa:** Sistemas Electromecánicos y Electrónica de Potencia — Laboratorio CRYSI  
**Versión del documento:** 2.0 · **Fecha de emisión:** septiembre de 2026  
**Naturaleza:** Documento de referencia académica para docencia de posgrado e investigación aplicada  
**Prerrequisitos:** Cálculo diferencial e integral, álgebra lineal, teoría de circuitos, fundamentos de máquinas eléctricas

---

### Resumen

El presente tratado expone, con rigor formal y vocación pedagógica, los fundamentos físicos y numéricos que sustentan la biblioteca de simulación `bloques_crysi`. Se articulan en un cuerpo coherente la teoría de integración de ecuaciones diferenciales ordinarias, la resolución de sistemas algebraico-diferenciales, las transformaciones de referencia para sistemas trifásicos, el modelado dinámico de máquinas eléctricas rotativas y el análisis de convertidores estáticos de energía. Cada modelo se deduce a partir de primeros principios, se expresa en forma canónica de espacio de estados y se vincula explícitamente con su implementación algorítmica en el núcleo compilado en C. El documento aspira a constituir un puente entre la formulación analítica clásica y su realización computacional determinista.

**Palabras clave:** Modelado dinámico; Integración numérica; Transformada de Park; Máquina síncrona de imanes permanentes; Máquina de inducción; Convertidores conmutados y promediados; Co-simulación EMT–fasorial.

---

### Índice

1. [Análisis Numérico y Resolución Temporal](#1-análisis-numérico-y-resolución-temporal)
2. [Teoría de Transformaciones en Sistemas Trifásicos](#2-teoría-de-transformaciones-en-sistemas-trifásicos)
3. [Modelado Dinámico de Máquinas Eléctricas](#3-modelado-dinámico-de-máquinas-eléctricas)
4. [Electrónica de Potencia y Convertidores Estáticos](#4-electrónica-de-potencia-y-convertidores-estáticos)
5. [Ingeniería de Vehículos Eléctricos y Almacenamiento](#5-ingeniería-de-vehículos-eléctricos-y-almacenamiento)
6. [Análisis Nodal Modificado y Subredes Eléctricas](#6-análisis-nodal-modificado-y-subredes-eléctricas)
7. [Interfaces de Instrumentación y Análisis de Resultados](#7-interfaces-de-instrumentación-y-análisis-de-resultados)
8. [Validación, Verificación y Calidad Numérica](#8-validación-verificación-y-calidad-numérica)
9. [Referencias Bibliográficas](#9-referencias-bibliográficas)

---

## 1. Análisis Numérico y Resolución Temporal

### 1.1 Naturaleza del Problema de Valor Inicial

Los sistemas físicos de interés se describen, tras la aplicación de las leyes de Kirchhoff, Newton y Faraday, como un sistema de ecuaciones diferenciales ordinarias (EDO) de primer orden, eventualmente acoplado a restricciones algebraicas:

$$
\dot{\mathbf{x}}(t) = \mathbf{f}\bigl(\mathbf{x}(t), \mathbf{u}(t), t\bigr), \qquad \mathbf{x}(t_0) = \mathbf{x}_0, \qquad \mathbf{y}(t) = \mathbf{h}(\mathbf{x}(t), \mathbf{u}(t))
$$

donde $\mathbf{x} \in \mathbb{R}^{n}$ denota el vector de estados (corrientes de inductores, tensiones de condensadores, velocidades y posiciones mecánicas), $\mathbf{u} \in \mathbb{R}^{m}$ el vector de entradas exógenas y $\mathbf{y} \in \mathbb{R}^{p}$ el vector de salidas observables. La simulación digital impone la discretización del eje temporal en instantes $t_k = k \cdot h$, con paso de integración $h = dt$.

### 1.2 Métodos de Integración Implementados

#### 1.2.1 Método de Euler Explícito (Primer Orden)

El esquema de Euler hacia adelante aproxima la evolución del estado mediante la expansión de Taylor truncada en el término lineal:

$$
\mathbf{x}_{k+1} = \mathbf{x}_k + h \cdot \mathbf{f}(\mathbf{x}_k, \mathbf{u}_k, t_k) + \mathcal{O}(h^2)
$$

**Propiedades formales:** Error local de truncamiento $\mathcal{O}(h^2)$, error global $\mathcal{O}(h)$, región de estabilidad absoluta limitada al disco $|1 + h\lambda| < 1$, donde $\lambda$ recorre el espectro del Jacobiano $\partial\mathbf{f}/\partial\mathbf{x}$. Su principal virtud reside en la mínima carga computacional por paso (una única evaluación de $\mathbf{f}$), lo que lo hace idóneo para prototipado y para sistemas no rígidos con $h$ suficientemente pequeño.

**Limitación:** En sistemas rígidos —caracterizados por constantes de tiempo dispares (p. ej., dinámicas eléctricas rápidas frente a mecánicas lentas)— la condición de estabilidad obliga a pasos extremadamente reducidos, incrementando el costo global.

#### 1.2.2 Método de Runge-Kutta de Cuarto Orden (RK4)

El método clásico de Runge-Kutta evalúa el campo vectorial en cuatro puntos del intervalo $[t_k, t_{k+1}]$:

$$
\begin{aligned}
\mathbf{k}_1 &= \mathbf{f}(\mathbf{x}_k, t_k) \\
\mathbf{k}_2 &= \mathbf{f}(\mathbf{x}_k + \tfrac{h}{2}\mathbf{k}_1, t_k + \tfrac{h}{2}) \\
\mathbf{k}_3 &= \mathbf{f}(\mathbf{x}_k + \tfrac{h}{2}\mathbf{k}_2, t_k + \tfrac{h}{2}) \\
\mathbf{k}_4 &= \mathbf{f}(\mathbf{x}_k + h\mathbf{k}_3, t_k + h) \\
\mathbf{x}_{k+1} &= \mathbf{x}_k + \tfrac{h}{6}(\mathbf{k}_1 + 2\mathbf{k}_2 + 2\mathbf{k}_3 + \mathbf{k}_4)
\end{aligned}
$$

**Propiedades formales:** Error local $\mathcal{O}(h^5)$, error global $\mathcal{O}(h^4)$, región de estabilidad significativamente más amplia. A igualdad de precisión, permite pasos entre 5 y 10 veces superiores a Euler, compensando el costo de cuatro evaluaciones por paso. Constituye el método por defecto para estudios de precisión y para la validación cruzada C–NumPy.

La selección del método se realiza mediante el parámetro `Modelo(metodo="rk4"|"euler")` y queda registrada en la cabecera del código C generado.

### 1.3 Causalidad Computacional y Lazos Algebraicos

#### 1.3.1 Definición Formal

En un formalismo de diagramas de bloques causal, cada bloque consume señales de entrada y produce señales de salida mediante una función estática $\mathbf{y} = \mathbf{g}(\mathbf{u})$ o dinámica $\dot{\mathbf{x}} = \mathbf{f}(\mathbf{x},\mathbf{u})$. Una interconexión cíclica entre bloques exclusivamente estáticos induce una ecuación implícita:

$$
\mathbf{y} = \mathbf{g}(\mathbf{y}, \mathbf{u}_{\text{ext}})
$$

Esta estructura, denominada *lazo algebraico*, carece de estados que rompan la dependencia temporal y, por tanto, no admite evaluación secuencial directa.

#### 1.3.2 Algoritmo de Resolución: Gauss-Seidel con Relajación Sucesiva

El núcleo implementa un método iterativo de punto fijo con relajación:

$$
\mathbf{y}^{(k+1)} = \mathbf{y}^{(k)} + \omega \bigl(\mathbf{g}(\mathbf{y}^{(k)}) - \mathbf{y}^{(k)}\bigr), \qquad 0 < \omega \leq 1
$$

donde $\omega = w_{\text{opt}}$ es el factor de relajación (SOR). La iteración se repite hasta satisfacer $\|\mathbf{y}^{(k+1)}-\mathbf{y}^{(k)}\|_{\infty} < \varepsilon$, con $\varepsilon = \text{tol}$ (típicamente $10^{-9}$) y cota superior $k_{\max} = \text{max\_iter}$ (50 por defecto). La no convergencia se notifica como excepción `RuntimeError` con diagnóstico del residuo y de la ventana temporal, evitando la propagación silenciosa de resultados espurios.

La detección estructural del lazo se basa en el análisis del grafo de dependencias y en la ordenación topológica de Kahn. La existencia de un ciclo entre bloques con *feedthrough* directo (p. ej., `Ganancia`, `Suma`, `PID` con $K_p \neq 0$) se valida en `Modelo._resolver()`, emitiendo un error constructivo si el lazo no es resoluble.

---

## 2. Teoría de Transformaciones en Sistemas Trifásicos

### 2.1 Motivación: Limitaciones del Control en el Dominio Temporal

Las magnitudes trifásicas equilibradas evolucionan como funciones sinusoidales $x_a(t) = X_m \cos(\omega t + \phi)$, $x_b(t) = X_m \cos(\omega t + \phi - 120^\circ)$, $x_c(t) = X_m \cos(\omega t + \phi + 120^\circ)$. Un regulador proporcional-integral (PI) clásico, cuya función de transferencia es $C(s)=K_p + K_i/s$, presenta ganancia infinita únicamente en $s=0$ (corriente continua). En consecuencia, resulta estructuralmente incapaz de seguir referencias sinusoidales con error nulo. La solución, debida a Park (1929) y Clarke (1943), consiste en cambiar el sistema de coordenadas.

### 2.2 Transformada de Clarke (Componentes $\alpha\beta0$)

Bajo la hipótesis de sistema equilibrado ($\sum x_{abc}=0$), la transformada de Clarke proyecta el espacio tridimensional sobre un plano ortogonal estacionario:

$$
\begin{bmatrix} x_{\alpha} \\ x_{\beta} \end{bmatrix} = \frac{2}{3}\begin{bmatrix} 1 & -\tfrac12 & -\tfrac12 \\ 0 & \tfrac{\sqrt{3}}{2} & -\tfrac{\sqrt{3}}{2} \end{bmatrix} \begin{bmatrix} x_a \\ x_b \\ x_c \end{bmatrix}, \qquad x_0 = \tfrac13(x_a+x_b+x_c)
$$

La componente homopolar $x_0$ se anula en condiciones equilibradas y se conserva únicamente para el análisis de fallas. La transformación es conservativa en amplitud (invariante en amplitud) y lineal, con inversa inmediata.

### 2.3 Transformada de Park (Componentes $dq0$)

La transformada de Park introduce un referencial rotatorio solidario al flujo rotórico, con ángulo eléctrico $\theta_e(t) = \int \omega_e(\tau)\,d\tau$:

$$
\begin{bmatrix} x_d \\ x_q \end{bmatrix} = \begin{bmatrix} \cos\theta_e & \sin\theta_e \\ -\sin\theta_e & \cos\theta_e \end{bmatrix} \begin{bmatrix} x_{\alpha} \\ x_{\beta} \end{bmatrix}
$$

**Consecuencia fundamental:** Una terna sinusoidal equilibrada de frecuencia $\omega_e$ se percibe, en el referencial $dq$, como un vector constante $(\bar{x}_d, \bar{x}_q)$. La regulación se reduce, por tanto, a un problema de control de variables continuas, donde el PI alcanza error nulo. La asignación canónica es: eje $d$ (directo) — control de flujo; eje $q$ (cuadratura) — control de par.

La implementación en `bloques_crysi` respeta la convención de Krause (amplitud invariante) y provee los bloques `Clarke`, `InvClarke`, `Park`, `InvPark` y `TransformadaQD` con validación de paridad C–NumPy.

---

## 3. Modelado Dinámico de Máquinas Eléctricas

### 3.1 Máquina Síncrona de Imanes Permanentes (PMSM/PMAC)

#### 3.1.1 Ecuaciones en el Referencial Síncrono $dq$

Adoptando la convención motora y referencial rotórico, las ecuaciones de tensión son:

$$
\begin{aligned}
v_d &= R_s i_d + L_d \frac{di_d}{dt} - \omega_e L_q i_q \\
v_q &= R_s i_q + L_q \frac{di_q}{dt} + \omega_e \bigl(L_d i_d + \lambda_m\bigr)
\end{aligned}
$$

donde $R_s$ es la resistencia estatórica, $L_d, L_q$ las inductancias síncronas, $\lambda_m$ el flujo concatenado de los imanes y $\omega_e = \tfrac{P}{2}\omega_m$ la pulsación eléctrica ($P$: número de polos).

**Acoplamiento cruzado:** Los términos $\omega_e L_q i_q$ y $\omega_e L_d i_d$ evidencian la interacción entre ejes. Un control FOC de alto desempeño requiere su compensación mediante acción *feedforward*, implementable con el bloque `Multiplicador`.

#### 3.1.2 Par Electromagnético

$$
T_e = \frac{3}{2}\,\frac{P}{2}\,\Bigl[\lambda_m i_q + (L_d - L_q)i_d i_q\Bigr]
$$

El primer sumando corresponde al par de alineación; el segundo, al par de reluctancia, explotable en máquinas de imanes interiores (IPMSM) mediante la estrategia de máximo par por amperio (MTPA). La dinámica mecánica se rige por:

$$
J \frac{d\omega_m}{dt} = T_e - T_L - B_m \omega_m, \qquad \frac{d\theta_m}{dt} = \omega_m
$$

#### 3.1.3 Instrumentación Asociada

*   `sensorPerdidasEstator()`: $P_{cu}=R_s(i_{as}^2+i_{bs}^2+i_{cs}^2)$ — pérdidas Joule instantáneas.
*   `sensorCorrienteRotor()`: Corrientes rotóricas referidas al estator $(i'_{ar},i'_{br},i'_{cr})$ mediante transformación inversa.
*   `velocidad_sincronica`: $\omega_s = 4\pi f / P$.

### 3.2 Máquina de Inducción Trifásica

El modelo en referencial estacionario $\alpha\beta$ (Krause, cap. 6) se expresa en términos de flujos concatenados:

$$
\begin{aligned}
\frac{d\boldsymbol{\lambda}_s}{dt} &= \mathbf{v}_s - R_s \mathbf{i}_s \\
\frac{d\boldsymbol{\lambda}_r}{dt} &= -R_r \mathbf{i}_r + j\omega_r \boldsymbol{\lambda}_r \\
T_e &= \frac{3}{2}\frac{P}{2}\,\Im\{\boldsymbol{\lambda}_s^{*} \mathbf{i}_s\}
\end{aligned}
$$

con relaciones constitutivas $\boldsymbol{\lambda}_s = L_s \mathbf{i}_s + L_m \mathbf{i}_r$, $\boldsymbol{\lambda}_r = L_r \mathbf{i}_r + L_m \mathbf{i}_s$. El modelo captura los fenómenos de deslizamiento, saturación incipiente y pérdidas en el cobre con fidelidad transitoria completa.

### 3.3 Máquina Síncrona de Polos Salientes y Máquina de Corriente Continua

Se incluyen formulaciones análogas para la máquina síncrona con devanado de excitación y para la máquina de corriente continua con excitación independiente, con parametrización directa $R_a, L_a, K_t, J, B_m$ y salidas de par, velocidad y posición.

---

## 4. Electrónica de Potencia y Convertidores Estáticos

### 4.1 Dicotomía: Modelos Conmutados y Promediados

| Atributo | Modelo Conmutado | Modelo Promediado |
|:---|:---|:---|
| Representación | Tren de pulsos discontinuo generado por comparación moduladora–portadora | Variable continua $d \in [0,1]$ (ciclo de trabajo) |
| Fidelidad | Reproduce rizado, armónicos de conmutación y pérdidas de conmutación | Filtra el rizado; conserva la dinámica de baja frecuencia |
| Requisito de paso | $h \ll (2f_{sw})^{-1}$ (típ. $10^{-6}$ s) | $h$ determinado por la dinámica del filtro $LC$ (típ. $10^{-5}$–$10^{-4}$ s) |
| Aplicación | Validación de EMI, diseño de filtros, HIL | Sintonía de lazos de control, estudios de sistema |

La biblioteca provee ambas variantes bajo un mismo bloque monolítico, conmutables mediante el parámetro `promediado`.

### 4.2 Convertidores DC–DC

Cada topología se modela como sistema de segundo orden con estados $i_L$ y $v_C$. A título ilustrativo, el convertidor Buck promediado:

$$
\begin{aligned}
L \frac{di_L}{dt} &= d\,V_{in} - v_C \\
C \frac{dv_C}{dt} &= i_L - \frac{v_C}{R}
\end{aligned}
$$

Las topologías Boost y Buck-Boost se deducen por dualidad, con la particularidad de la no linealidad bilineal $ (1-d)v_C$. La encapsulación monolítica evita la formación de lazos algebraicos irresolubles que surgirían al interconectar un `Inductor` ($V \to I$) en serie con un `InterruptorIdeal` mediante puertos causales.

### 4.3 Inversores DC–AC

#### 4.3.1 Puente Monofásico

Con modulación bipolar, la tensión diferencial promediada es $v_{out}(t) = m_a V_{dc} \sin(\omega t)$, con índice de modulación $m_a \in [-1,1]$.

#### 4.3.2 Puente Trifásico de Dos Niveles

La tensión de fase promediada, con inyección de secuencia cero, se expresa como:

$$
v_{kN} = \frac{V_{dc}}{2} m_k, \qquad k \in \{a,b,c\}, \qquad v_{k0} = v_{kN} - \tfrac13\sum_j v_{jN}
$$

La variante conmutada genera los estados discretos $s_k \in \{0,1\}$ y $v_{kN} = V_{dc}(2s_k - s_j - s_l)/3$. La biblioteca implementa las estrategias SPWM y SVPWM, con límite de saturación vectorial $V_{\max}=V_{dc}/\sqrt{3}$.

---

## 5. Ingeniería de Vehículos Eléctricos y Almacenamiento

### 5.1 Dinámica Longitudinal

Conforme a la segunda ley de Newton, la aceleración del centro de masa satisface:

$$
m_{eq} \frac{dV}{dt} = F_t - \bigl(F_{aero} + F_{rod} + F_{pend}\bigr)
$$

donde $F_t = T_m i_g \eta_g / r_w$ (fuerza tractiva), $F_{aero}=\tfrac12 \rho C_d A V|V|$ (arrastre aerodinámico cuadrático), $F_{rod}=C_{rr} m g \cos\alpha$ (rodadura) y $F_{pend}=m g \sin\alpha$ (pendiente). La masa equivalente $m_{eq}$ incorpora la inercia rotacional del tren motriz.

### 5.2 Modelo Electroquímico Equivalente de Batería (ECM R0–2RC)

El modelo adopta una arquitectura de circuito equivalente de segundo orden con acoplamiento térmico:

1.  **Tensión de circuito abierto (OCV):** Función no lineal $V_{oc}=f(\text{SOC})$ tabulada mediante *lookup table* unidimensional con interpolación lineal y extrapolación constante. La tabla por defecto reproduce la característica de una celda NMC de 4,2 V.
2.  **Polarización y difusión:** Dos redes $R_1$–$C_1$ y $R_2$–$C_2$ en serie con $R_0$ modelan las constantes de tiempo electroquímicas (segundos y minutos, respectivamente): $\dot{v}_{1}= (i - v_1/R_1)/C_1$.
3.  **Balance térmico:** $C_{th}\dot{T}_{cell}= P_{loss} - (T_{cell}-T_{amb})/R_{th}$, con $P_{loss}=i^2R_0 + v_1^2/R_1 + v_2^2/R_2$.
4.  **Gestión operativa:** Límites de corriente de carga/descarga dependientes de SOC y temperatura, con histéresis y degradación cíclica $\Delta Q = k_{deg} \cdot |i|dt$.

El bloque `BateriaECM` expone seis salidas: tensión terminal, SOC, temperatura, potencia disipada y límites de corriente, habilitando el diseño de estrategias de gestión energética con restricciones reales.

### 5.3 Saturación Vectorial y SVPWM

Cuando la referencia de tensión $(v_d^*, v_q^*)$ excede el hexágono de factibilidad del inversor, el bloque `SaturarVectorial` aplica una homotecia:

$$
(v_d, v_q) = \frac{V_{\max}}{\sqrt{{v_d^*}^2+{v_q^*}^2}}(v_d^*, v_q^*), \quad \text{si } \| \mathbf{v}^*\| > V_{\max}
$$

La preservación de la fase resulta crítica: una alteración angular implicaría una orientación subóptima del vector de corriente, con la consiguiente reducción del par y generación de armónicos de orden bajo.

---

## 6. Análisis Nodal Modificado y Subredes Eléctricas

### 6.1 Motivación

Los formalismos causales resultan insuficientes para topologías eléctricas arbitrarias con mallas, interruptores y diodos ideales. El Análisis Nodal Modificado (MNA) proporciona un marco no causal que asamblea sistemáticamente las ecuaciones de Kirchhoff (KCL, KVL) y las relaciones constitutivas en un sistema algebraico-diferencial de la forma $\mathbf{G}\mathbf{x} + \mathbf{C}\dot{\mathbf{x}} = \mathbf{b}(t)$.

### 6.2 Formulación e Implementación

El bloque `SubredMNA` encapsula una subred eléctrica descrita mediante primitivas `Nodo`, `Resistor`, `Capacitor`, `Inductor`, `Fuente de Tensión/Corriente`, `Interruptor` y `Diodo Ideal`. La estampación (*stamping*) sigue la metodología de Dommel para elementos dinámicos (integración trapezoidal o Euler hacia atrás) y el tratamiento de complementariedad para diodos ($v_d \geq V_f \perp i_d \geq 0$). La resolución por paso combina factorización LU dispersa con iteración de punto fijo para la conmutación de semiconductores, con detección de cruces por cero y adaptación del paso.

La subred expone al grafo causal exterior únicamente las tensiones nodales y corrientes de rama seleccionadas, preservando la compatibilidad con el motor de flujo de señales.

---

## 7. Interfaces de Instrumentación y Análisis de Resultados

### 7.1 Adquisición y Diezmado

La clase `Scope` implementa una estrategia de adquisición con ventana deslizante y diezmado adaptativo (*step*, *LTTB*, *min-max*) que garantiza la preservación de extremos locales con complejidad $\mathcal{O}(n)$. El modo XY habilita la representación de retratos de fase y características par–velocidad.

### 7.2 Métricas Automáticas

El objeto `Resultado` provee operadores de análisis pos-simulación: valor de pico, valor medio en ventana, tiempo de establecimiento al 2 %, tasa de distorsión armónica (THD) mediante FFT con ventana de Hann, y transformada rápida de Fourier unilateral. Dichas métricas se emplean en la validación de criterios de desempeño y en la sintonización de controladores.

---

## 8. Validación, Verificación y Calidad Numérica

La estrategia de aseguramiento de la calidad se articula en tres niveles:

1.  **Paridad C–NumPy:** Cada bloque dinámico posee una implementación de referencia en `backend_numpy.py`. La suite de regresión (314 pruebas) verifica la equivalencia con tolerancias $\text{rtol}=10^{-8}$, $\text{atol}=10^{-9}$ en trayectorias completas.
2.  **Soluciones analíticas:** Casos con solución cerrada (circuito RLC serie, lazo estático $y = G\cdot r/(1+G)$) se emplean como oráculos.
3.  **Benchmarks bibliográficos:** Reproducción de los casos de Krause (máquina síncrona), Mohan (convertidor Buck) y perfiles de conducción WLTP (vehículo).

La ejecución se realiza mediante `pytest`, con cobertura de ramas y reporte de no convergencia de lazos algebraicos como fallo explícito.

---

## 9. Referencias Bibliográficas

*   Clarke, E. (1943). *Circuit Analysis of A-C Power Systems: Symmetrical and Related Components*. Wiley.
*   Dommel, H. W. (1969). Digital computer solution of electromagnetic transients in single- and multiphase networks. *IEEE Transactions on Power Apparatus and Systems*, 88(4), 388–399.
*   Krause, P. C., Wasynczuk, O., Sudhoff, S. D., & Pekarek, S. (2013). *Analysis of Electric Machinery and Drive Systems* (3.ª ed.). Wiley-IEEE Press.
*   Mohan, N., Undeland, T. M., & Robbins, W. P. (2003). *Power Electronics: Converters, Applications, and Design* (3.ª ed.). Wiley.
*   Park, R. H. (1929). Two-reaction theory of synchronous machines: Generalized method of analysis—Part I. *Transactions of the AIEE*, 48(3), 716–727.
*   Rashid, M. H. (2017). *Power Electronics Handbook* (4.ª ed.). Butterworth-Heinemann.
*   Ehsani, M., Gao, Y., Longo, S., & Ebrahimi, K. (2018). *Modern Electric, Hybrid Electric, and Fuel Cell Vehicles* (3.ª ed.). CRC Press.

---

*Fin del tratado. Para la especificación de la interfaz de co-simulación con flujos de potencia, véase `diseno_red.md`.*

