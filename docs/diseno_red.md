# Especificación Formal de Diseño: Integración entre el Simulador Electromecánico bloques_crysi y Backends de Flujo de Potencia

## Co-simulación EMT–Fasorial con Interfaces Desacopladas

**Documento:** Especificación de Requisitos y Diseño de Arquitectura (ERDA) — Módulo `red`  
**Versión:** 1.0 · **Estado:** Aprobado para implementación · **Clasificación:** Especificación técnica formal  
**Autoría:** Laboratorio CRYSI — Área de Sistemas Eléctricos y Co-simulación  
**Fecha de emisión:** septiembre de 2026 · **Revisión:** Anual o ante cambio mayor de contrato

---

### Resumen Ejecutivo

La presente especificación define, con carácter normativo, la arquitectura de integración entre el simulador de transitorios electromecánicos `bloques_crysi` —de naturaleza EMT, con resolución micro–milisegunda— y los motores de flujo de potencia cuasi-estático (pandapower, OpenDSS) que operan en régimen fasorial. El diseño persigue tres objetivos fundamentales: (i) garantizar el desacoplamiento estricto entre el núcleo EMT y los backends de red, (ii) establecer un sistema único y canónico de unidades físicas en el Sistema Internacional, y (iii) proveer un orquestador determinista que asegure la convergencia de la co-simulación mediante ventanas de comunicación con relajación. El documento constituye el contrato inmutable (Fase 0) sobre el cual se desarrollan las fases subsiguientes de implementación y validación.

**Palabras clave:** Co-simulación; Flujo de potencia; Interfaz backend-agnóstica; Sistema Internacional de Unidades; Orquestación por ventanas; Convergencia numérica.

---

### Índice

1. [Introducción y Alcance](#1-introducción-y-alcance)
2. [Principios Rectores de Diseño](#2-principios-rectores-de-diseño)
3. [Sistema Canónico de Unidades](#3-sistema-canónico-de-unidades)
4. [Contrato Formal de la Interfaz BackendRed](#4-contrato-formal-de-la-interfaz-backendred)
5. [Estructura Modular y Organización de Paquetes](#5-estructura-modular-y-organización-de-paquetes)
6. [Orquestador de Co-simulación](#6-orquestador-de-co-simulación)
7. [Metodología de Validación y Criterios de Aceptación](#7-metodología-de-validación-y-criterios-de-aceptación)
8. [Gestión de Dependencias y Versionado](#8-gestión-de-dependencias-y-versionado)
9. [Análisis de Riesgos y Estrategias de Mitigación](#9-análisis-de-riesgos-y-estrategias-de-mitigación)
10. [Hoja de Ruta y Plan de Trabajo](#10-hoja-de-ruta-y-plan-de-trabajo)
11. [Trazabilidad con el Núcleo Existente](#11-trazabilidad-con-el-núcleo-existente)
12. [Referencias Normativas y Bibliográficas](#12-referencias-normativas-y-bibliográficas)

---

## 1. Introducción y Alcance

### 1.1 Contexto y Motivación

La simulación integral de sistemas eléctricos modernos —que comprenden generación distribuida, accionamientos de velocidad variable y recursos de almacenamiento— exige la conciliación de dos escalas temporales disímiles. Por un lado, los fenómenos electromagnéticos transitorios (corrientes de magnetización, conmutaciones de semiconductores, dinámicas de control) requieren pasos de integración del orden de microsegundos. Por otro, el análisis de redes de distribución y transporte se efectúa tradicionalmente mediante flujos de potencia fasoriales, con resolución de segundos a minutos y bajo hipótesis de régimen sinusoidal permanente.

La co-simulación EMT–fasorial emerge como la solución metodológicamente rigurosa para capturar la interacción bidireccional entre ambas escalas sin incurrir en la prohibitiva carga computacional de una simulación EMT monolítica de la red completa.

### 1.2 Alcance del Presente Documento

Esta especificación cubre exclusivamente la **Fase 0: Fundaciones**, que comprende:

*   La definición del contrato abstracto `BackendRed` y de los tipos de datos inmutables asociados.
*   La formalización del sistema de unidades y de las convenciones de signo.
*   La estructura de paquetes y la política de dependencias opcionales.
*   La especificación funcional del orquestador `CoSimuladorRed` (Fase 3), a modo de vista previa normativa.

Quedan fuera del alcance la implementación concreta de los adaptadores para pandapower (Fase 1) y OpenDSS (Fase 2), así como las extensiones trifásicas desequilibradas y multi-PCC (Fase 4).

### 1.3 Audiencia y Uso Normativo

El documento se dirige al equipo de desarrollo, a revisores externos y a futuros mantenedores. Su contenido posee carácter prescriptivo: toda implementación que aspire a la conformidad debe satisfacer íntegramente los requisitos aquí consignados. Cualquier modificación del contrato `BackendRed` exige el incremento de la versión mayor del paquete y un ciclo completo de revisión formal.

---

## 2. Principios Rectores de Diseño

La arquitectura se fundamenta en los siguientes principios, cuya observancia es de cumplimiento obligatorio:

| Principio | Formulación Normativa | Justificación |
|:---|:---|:---|
| **Separación de responsabilidades** | El simulador `bloques_crysi` asume exclusivamente la simulación EMT con paso $h \in [10^{-6},10^{-3}]$ s. El backend de red resuelve el flujo de potencia cuasi-estático mediante instantáneas desacopladas con periodo de comunicación $H \in [10^{-1},1]$ s, con $H \gg h$. | Evita la contaminación de modelos y permite la sustitución independiente de cada subsistema. |
| **Abstracción agnóstica al backend** | Un único orquestador, `CoSimuladorRed`, interactúa con los backends a través de la clase abstracta `BackendRed`. Las implementaciones concretas (`BackendPandapower`, `BackendOpenDSS`) constituyen meros adaptadores (*plugins*) intercambiables. | Garantiza la extensibilidad y la testabilidad mediante dobles de prueba (*mocks*). |
| **Unidades SI como única interfaz pública** | La totalidad de la API pública de `BackendRed` se expresa exclusivamente en unidades del Sistema Internacional (V, A, W, VAr, rad). Las conversiones a unidades nativas (p.u., MW, kV) se encapsulan internamente. | Elimina ambigüedades dimensionales y centraliza la trazabilidad metrológica. |
| **Dependencias opcionales y núcleo autónomo** | Las bibliotecas `pandapower` y `opendssdirect` se declaran como dependencias opcionales (`extras_require`). El núcleo `bloques_crysi` no las importa en tiempo de carga. | Preserva la instalabilidad mínima y evita el acoplamiento obligatorio a ecosistemas externos. |
| **Inmutabilidad del contrato** | La interfaz `BackendRed`, una vez aprobada en Fase 0, permanece congelada. Toda evolución compatible se realiza mediante adición, nunca mediante modificación, y toda ruptura exige versión mayor. | Asegura la estabilidad de las fases 1–3 y la confianza de los consumidores de la API. |
| **Coherencia lingüística y estilística** | La nomenclatura sigue las convenciones del proyecto: `PascalCase` para clases, `snake_case` para módulos y funciones, identificadores en español para conceptos de dominio. | Mantiene la homogeneidad del código base y su alineación con la documentación académica en español. |

---

## 3. Sistema Canónico de Unidades

### 3.1 Fundamento Metrológico

La heterogeneidad de convenciones entre herramientas de flujo de potencia constituye una fuente histórica de errores sistemáticos. Pandapower adopta el sistema por unidad (p.u.) con base $S_{base}$ y $V_{base}$; OpenDSS emplea kV, kW y kVAr con tensiones línea–neutro por defecto en ciertos contextos. La presente especificación impone el Sistema Internacional como única fuente de verdad en la frontera entre `bloques_crysi` y los backends.

### 3.2 Tabla de Correspondencias

| Magnitud física | Representación canónica (bloques_crysi) | Representación nativa pandapower | Representación nativa OpenDSS | Observaciones |
|:---|:---|:---|:---|:---|
| Tensión | V, fase–neutro, [V] | p.u. sobre $V_{base}$ [kV LL] | kV, línea–neutro por defecto | Conversión LN–LL mediante $\sqrt{3}$ |
| Potencia activa | W, [W] | MW, [MW] | kW, [kW] | Factor $10^{6}$ / $10^{3}$ |
| Potencia reactiva | VAr, [VAr] | MVAr, [MVAr] | kVAr, [kVAr] | Idem |
| Potencia aparente | VA, [VA] | MVA, [MVA] | kVA, [kVA] | Idem |
| Ángulo de fase | rad, [rad] | grados, [°] | grados, [°] | Factor $\pi/180$ |
| Frecuencia | Hz, [Hz] | Hz, [Hz] | Hz, [Hz] | Sin conversión |
| Impedancia | $\Omega$, [$\Omega$] | p.u. | $\Omega$ | Base $Z_{base}=V_{base}^2/S_{base}$ |

### 3.3 Reglas Normativas de Conversión

1.  **Centralización absoluta:** La totalidad de las transformaciones dimensionales reside en el módulo `red/unidades.py`. Queda proscrita cualquier operación aritmética de conversión dispersa en adaptadores ($*10^{6}$, $/10^{3}$, etc.).
2.  **Encapsulación de base:** La clase `SistemaUnidades(v_base_kv: float, s_base_mva: float)` encapsula la base del sistema y provee métodos `a_pu()`, `a_si()`, `v_ll_a_ln()`, con validación de rangos y propagación de incertidumbres.
3.  **Convención de signo de potencia:** Se adopta la convención de carga como referencia canónica:
    *   **Carga:** $P > 0$, $Q > 0$ denota consumo de potencia activa y reactiva inductiva.
    *   **Generación:** En la interfaz `BackendRed`, `set_generacion(bus, P, Q)` emplea $P > 0$, $Q > 0$ para inyección. El adaptador traduce al signo nativo del backend (p. ej., inyección negativa en pandapower si este adopta convención de carga).
    *   La distinción entre `set_carga` y `set_generacion` es semántica y no meramente aritmética, habilitando validaciones de coherencia.
4.  **Tensión de retorno:** `BackendRed.get_tension(bus)` retorna una instancia inmutable `TensionBus` con magnitud fase–neutro en voltios y ángulo en radianes. El atributo `es_linea_linea` indica si la magnitud original del backend era línea–línea, a efectos de trazabilidad.

---

## 4. Contrato Formal de la Interfaz BackendRed

### 4.1 Definición Abstracta

```python
from abc import ABC, abstractmethod
from dataclasses import dataclass

class BackendRed(ABC):
    """Contrato abstracto para backends de flujo de potencia.

    Todas las magnitudes en SI. Implementaciones deben ser reentrantes
    y no retener estado mutable entre invocaciones de runpp() salvo el
    modelo de red subyacente.
    """

    @abstractmethod
    def set_carga(self, bus_idx: int, p_w: float, q_var: float) -> None:
        """Establece la demanda en el bus indicado (convención carga)."""
        ...

    @abstractmethod
    def set_generacion(self, bus_idx: int, p_w: float, q_var: float) -> None:
        """Establece la inyección en el bus indicado (convención generación)."""
        ...

    @abstractmethod
    def runpp(self) -> None:
        """Ejecuta el flujo de potencia. Lanza ErrorFlujoPotencia si no converge."""
        ...

    @abstractmethod
    def get_tension(self, bus_idx: int) -> "TensionBus":
        """Retorna la tensión nodal tras runpp(). Requiere runpp() previo."""
        ...

    @abstractmethod
    def get_corriente_linea(self, linea_idx: int) -> float:
        """Retorna la magnitud de corriente en la línea indicada [A]."""
        ...

    @property
    @abstractmethod
    def nombre_backend(self) -> str:
        """Identificador canónico del backend (p. ej., 'pandapower', 'opendss')."""
        ...
```

### 4.2 Tipos de Datos Inmutables

#### 4.2.1 TensionBus

```python
from dataclasses import dataclass

@dataclass(frozen=True, slots=True)
class TensionBus:
    magnitud_v: float       # Tensión fase–neutro [V], salvo es_linea_linea
    angulo_rad: float       # Ángulo de fase [rad], referencia slack = 0
    es_linea_linea: bool = False  # Verdadero si la magnitud original era LL
```

La inmutabilidad (`frozen=True`) garantiza la ausencia de efectos laterales y la aptitud para su uso como clave en estructuras de datos funcionales.

#### 4.2.2 PotenciaInyectada

```python
@dataclass(frozen=True, slots=True)
class PotenciaInyectada:
    p_w: float      # Potencia activa [W], >0 = generación/inyección
    q_var: float    # Potencia reactiva [VAr], >0 = inyección capacitiva
```

### 4.3 Taxonomía de Excepciones

| Excepción | Condición de Lanzamiento | Atributos | Tratamiento Recomendado |
|:---|:---|:---|:---|
| `ErrorFlujoPotencia` | No convergencia de `runpp()` tras iteraciones del backend | `backend: str`, `detalles: str`, `iteraciones: int` | El orquestador decide: reducir $H$, reintentar con relajación, abortar con diagnóstico |
| `ValueError` | Índice de bus/línea fuera de rango, magnitud no física | — | Error de programación; no recuperable |
| `RuntimeError` | Invocación de `get_tension` sin `runpp()` previo | — | Violación de protocolo; no recuperable |

---

## 5. Estructura Modular y Organización de Paquetes

```
bloques_crysi/
├── src/bloques_crysi/red/
│   ├── __init__.py                 # Re-exportación pública y registro de backends
│   ├── unidades.py                 # Sistema SI ↔ p.u./kV/MW (fuente única de conversión)
│   ├── backend_base.py             # Definición de BackendRed, TensionBus, BackendRedMock
│   ├── adaptador_pandapower.py     # Adaptador Fase 1 (dependencia opcional)
│   ├── adaptador_opendss.py        # Adaptador Fase 2 (dependencia opcional)
│   └── co_simulador.py             # Orquestador Fase 3 (CoSimuladorRed)
└── tests/red/
    ├── test_unidades.py
    ├── test_backend_base.py
    ├── test_adaptador_pandapower.py
    ├── test_adaptador_opendss.py
    └── test_co_simulador.py
```

**Regla de importación:** Ningún módulo del núcleo (`modelo.py`, `bloques.py`) importa `red.*` en tiempo de carga. La co-simulación se activa exclusivamente mediante importación explícita del consumidor, preservando la autonomía del núcleo EMT.

---

## 6. Orquestador de Co-simulación

### 6.1 Especificación Funcional

```python
class CoSimuladorRed:
    """Orquestador de co-simulación EMT–fasorial por ventanas de comunicación."""

    def __init__(
        self,
        modelo: Modelo,
        backend: BackendRed,
        bus_pcc: int,
        medidor: MedidorPotencia,
        fuente_red: FuenteTrifasica,
        dt_red: float = 0.1,                # Periodo de comunicación [s]
        tol_convergencia_v: float = 1e-3,   # Tolerancia relativa en tensión [p.u.]
        max_iter_ventana: int = 10,         # Iteraciones máximas por ventana
        relajacion: float = 0.5,            # Factor de sub-relajación α ∈ (0,1]
    ) -> None: ...

    def run(self, t_fin: float) -> ResultadoCoSim: ...
```

Los parámetros `dt_red`, `tol_convergencia_v` y `relajacion` se exponen como parte del contrato público y quedan sujetos a validación de rangos ($dt_{red} > 10 \cdot dt_{emt}$, $0 < \alpha \leq 1$).

### 6.2 Algoritmo Formal por Ventana de Comunicación

Para cada ventana $W_k = [t_k, t_k + H]$, con $H = dt_{red}$, se ejecuta el siguiente procedimiento iterativo de punto fijo con sub-relajación:

```
Entrada: Estado EMT en t_k, Backend en estado convergido previo
Salida: Estado EMT en t_{k+1}, Tensión PCC actualizada

1. Simular EMT en W_k:
   (P_k, Q_k) ← promedio temporal de MedidorPotencia en W_k
                obtenido mediante Modelo.run_generador(H, chunk=H)

2. Transferir a red:
   backend.set_carga(bus_pcc, P_k, Q_k)
   backend.runpp()  // puede lanzar ErrorFlujoPotencia

3. Adquirir tensión:
   V_k ← backend.get_tension(bus_pcc)  // TensionBus en SI

4. Sub-relajación:
   V_k^{(\alpha)} ← α · V_k + (1-α) · V_{k-1}^{(\alpha)}

5. Retroalimentar a EMT:
   Modelo.set_param(fuente_red, "amplitud", V_k^{(\alpha)}.magnitud_ll_v)
   Modelo.set_param(fuente_red, "fase", V_k^{(\alpha)}.angulo_rad)

6. Criterio de convergencia:
   Si ‖V_k^{(\alpha)} - V_{k-1}^{(\alpha)}‖ / V_{nom} < ε_v → avanzar a W_{k+1}
   En caso contrario, repetir desde (1) con V_{k-1} ← V_k^{(\alpha)}
   hasta agotar max_iter_ventana → emitir advertencia y avanzar con último valor
```

La sub-relajación con $\alpha = 0{,}5$ (valor por defecto) garantiza la estabilidad del acoplamiento débil, mitigando oscilaciones numéricas cuando $H$ es comparable a las constantes de tiempo de la red.

### 6.3 Inicialización (*Bootstrapping*)

Previo a la primera invocación de `runpp()`, la fuente trifásica `FuenteTrifasica` opera con los valores paramétricos iniciales especificados en su construcción. Se prescribe inicializarla con la tensión nominal línea–línea del punto de acoplamiento común (PCC):

$$
V_{init} = V_{nom,LL} = \sqrt{3} \cdot V_{nom,LN}
$$

Esta elección minimiza el transitorio numérico de la primera ventana.

### 6.4 Requisito de Mutabilidad Controlada

La actualización de parámetros en tiempo de ejecución se realiza exclusivamente mediante `Modelo.set_param(bloque, nombre_parametro, valor)`, que valida índices, rangos y tipos antes de escribir en los arreglos `_param_arrays` del núcleo C. Queda proscrito el acceso directo a dichos arreglos, por el riesgo de corrupción de memoria y violación de la encapsulación.

---

## 7. Metodología de Validación y Criterios de Aceptación

### 7.1 Circuito de Referencia Canónico

A fin de garantizar la paridad entre backends, se define un circuito de referencia analíticamente tratable:

```
Bus 0 (Slack) ── Línea (R = 0,10 Ω, X = 0,05 Ω) ── Bus 1 (Carga: P = 100 kW, Q = 50 kVAr)
V0 = 230 V LN (400 V LL), f = 50 Hz, Sbase = 1 MVA, Vbase = 0,4 kV LL
```

**Resultados analíticos esperados (flujo monofásico equivalente):**

*   $V_1 \approx 224{,}3\ \text{V LN}$, $\angle V_1 \approx -1{,}2^{\circ}$
*   $I_{línea} \approx 289\ \text{A}$ (magnitud)

### 7.2 Criterios de Aceptación

| Prueba | Backend | Métrica | Tolerancia |
|:---|:---|:---|:---|
| `test_adaptador_pandapower` | pandapower | $\|V_{sim} - V_{ref}\| / V_{ref}$ | $< 1\%$ |
| `test_adaptador_opendss` | OpenDSS | Idem | $< 1\%$ |
| `test_co_simulador` (ventana única) | Mock | Convergencia en $< 5$ iteraciones | Determinista |
| `test_unidades` | — | Round-trip SI→p.u.→SI | Error $< 10^{-12}$ |

La no superación de cualquiera de estos criterios se considera fallo bloqueante para la promoción de la fase correspondiente.

---

## 8. Gestión de Dependencias y Versionado

### 8.1 Declaración en pyproject.toml

```toml
[project.optional-dependencies]
red-pandapower = ["pandapower>=2.14,<3.0"]
red-opendss    = ["opendssdirect.py>=0.9,<1.0"]
red            = ["bloques_crysi[red-pandapower]", "bloques_crysi[red-opendss]"]

[project.urls]
"Documentación" = "https://crysi.github.io/bloques_crysi/red/"
```

La cota superior `<3.0` / `<1.0` previene la incorporación automática de versiones mayores con ruptura potencial de API, sin perjuicio de su evaluación explícita en ciclos de mantenimiento.

### 8.2 Política de Versionado Semántico

*   **Versión mayor (X.0.0):** Modificación del contrato `BackendRed` o de las convenciones de unidades.
*   **Versión menor (0.X.0):** Adición de backends o funcionalidades compatibles (p. ej., soporte trifásico desequilibrado).
*   **Versión de parche (0.0.X):** Correcciones de adaptadores sin cambio de interfaz.

---

## 9. Análisis de Riesgos y Estrategias de Mitigación

| Riesgo Identificado | Probabilidad | Impacto | Estrategia de Mitigación | Responsable |
|:---|:---|:---|:---|:---|
| Exposición directa de `_param_arrays` → corrupción de memoria C | Media | Crítico | Encapsulación mediante `Modelo.set_param()` con validación de límites y tipos; pruebas de inyección de fallos | Núcleo |
| Periodo de comunicación $H$ excesivo → oscilación del punto fijo | Media | Alto | Sub-relajación $\alpha$ configurable (0,5 por defecto); validación $H \leq 10 \cdot \tau_{red}$; log de residuo por ventana | Orquestador |
| No convergencia del flujo en una ventana | Baja | Alto | Propagación de `ErrorFlujoPotencia`; política de reintento con $H$ reducido o aborto con diagnóstico | Orquestador |
| Ruptura de compatibilidad por actualización de pandapower/OpenDSS | Media | Medio | Dependencias con cota mínima probada; suite de paridad ejecutada en CI contra versiones *latest* | DevOps |
| Limitación a red equilibrada en v1.0 | Alta | Medio | Alcance explícitamente documentado; extensión trifásica desequilibrada planificada como Fase 4 con modelo de secuencia | Diseño |

---

## 10. Hoja de Ruta y Plan de Trabajo

| Fase | Entregables | Criterio de Entrada | Criterio de Salida |
|:---|:---|:---|:---|
| **Fase 0** | Aprobación de la presente ERDA; stubs `unidades.py`, `backend_base.py` con `BackendRedMock` | Revisión por pares del documento | Firma de aprobación y congelación del contrato |
| **Fase 1** | `adaptador_pandapower.py` + `test_adaptador_pandapower.py` (paridad <1 %) | Fase 0 aprobada | Pruebas de paridad superadas en CI |
| **Fase 2** | `adaptador_opendss.py` + `test_adaptador_opendss.py` (paralela a Fase 1) | Fase 0 aprobada | Idem |
| **Fase 3** | `co_simulador.py` + `Modelo.set_param()` + ejemplo end-to-end (motor + red IEEE-13) | Fases 1–2 completadas | Ejemplo reproducible con convergencia documentada |
| **Fase 4** | Soporte trifásico desequilibrado, multi-PCC, perfiles de carga anual | Fase 3 estable | Validación contra OpenDSS en modo armónico |

---

## 11. Trazabilidad con el Núcleo Existente

| Funcionalidad Requerida | Localización en Código Base | Relevancia para Co-simulación |
|:---|:---|:---|:---|
| `Modelo.run_generador(t, chunk)` | `src/bloques_crysi/modelo.py:380` | Ejecución por ventanas con guardado/restauración de estado EMT |
| `Modelo._param_arrays` | `src/bloques_crysi/modelo.py:120` | Base para `set_param()` (mutación segura) |
| `MedidorPotencia` | `src/bloques_crysi/medidores.py` | Cálculo de $P$ y $Q$ promedio en PCC |
| `FuenteTrifasica` | `src/bloques_crysi/fuentes.py` | Inyección de tensión PCC en EMT |
| `TransformadaQD` | `src/bloques_crysi/transformadas.py` | Medición fasorial para control FOC en el lado EMT |
| `crear_control_foc` | `src/bloques_crysi/subsistemas.py` | Caso de carga EMT realista para validación integrada |

---

## 12. Referencias Normativas y Bibliográficas

*   IEEE Std 1459-2010. *Definitions for the Measurement of Electric Power Quantities Under Sinusoidal, Nonsinusoidal, Balanced, or Unbalanced Conditions*.
*   Milano, F. (2010). *Power System Modelling and Scripting*. Springer.
*   Monti, A., et al. (2008). Co-simulation of heterogeneous power systems. *IEEE Transactions on Power Systems*, 23(4), 1747–1755.
*   Palensky, P., et al. (2017). Co-simulation of intelligent power systems. *IEEE Industrial Electronics Magazine*, 11(1), 34–50.

---

**Aprobación:**

| Rol | Nombre | Firma | Fecha |
|:---|:---|:---|:---|
| Arquitecto de Sistemas | Equipo bloques_crysi | — | 2026-09-02 |
| Revisor Técnico | — | — | — |
| Responsable de Calidad | — | — | — |

*Fin de la Especificación Formal de Diseño — Fase 0.*

