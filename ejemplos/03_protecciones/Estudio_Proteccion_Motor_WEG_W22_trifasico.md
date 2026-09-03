# Estudio de Protección — Motor WEG W22 7.5 HP (Trifásico)
### Fuente 13.2 kV → Transformador 500 kVA → Cable 8 AWG → Fusible + Relé 50/51 → Motor (00718ET3E213T-W22)


## Diagrama Unifilar

```mermaid
graph LR
    A["Fuente<br/>13.2 kV"] -->|Transformador<br/>500 kVA, 13,200/220-127 V<br/>Z = 4.5%, X/R = 5| B["Barra 220 V"]
    B -->|Cable 8 AWG<br/>≈19.6 m<br/>R=2.56 Ω/km, X=0.171 Ω/km| C["Fusible<br/>FRN-R-35<br/>(dual-element)"]
    C --> D["Relé<br/>50 / 51"]
    D --> E(("Motor WEG W22<br/>7.5 HP, 220V, 3Ø, 4 polos<br/>FLA=19.46A, SF=1.25"))
```

---

## Datos del Sistema

| Parámetro | Valor |
|---|---|
| Transformador | 500 kVA |
| Tensión secundaria | **220/127 V** (13,200–220/127 V, dato de placa del transformador) |
| Impedancia del transformador | 4.5% (dentro del rango típico 4–5.5%, NMX-J-205) |
| X/R asumido (transformador) | 5 |
| **Conductor** | **8 AWG** — R = 2.56 Ω/km, X = 0.171 Ω/km (NOM-001-SEDE-2012, Tabla 9, conducto PVC, cobre sin recubrimiento, 75 °C) — **verificado contra Tabla 9, renglón 8 AWG, columna PVC: coincide exactamente** |
| **Motor** | **WEG W22, modelo 00718ET3E213T-W22, 7.5 HP (5.5 kW), NEMA, TRIFÁSICO** |
| Tensión / polos | 208-230/460 V, 3Ø, 4 polos, 1765 RPM, 60 Hz, marco 213/5T |
| Corriente nominal de placa | 208 V: 20.5 A · 230 V: 18.6 A · 460 V: 9.29 A *(ver nota de verificación arriba)* |
| **Corriente nominal (interpolada a 220 V)** | **19.46 A** |
| Corriente de arranque de placa (LRC 7.1×, Código H) | 208 V: 145.55 A (calc.) · 230 V: 132.06 A (calc.) |
| **Corriente de arranque (interpolada a 220 V)** | **138.19 A** |
| Eficiencia / FP (dato de placa, 100% carga) | 91.7% / 0.81 |
| Factor de servicio | **1.25** |
| Protección / aislamiento | IP55, Clase F |

**Metodología de interpolación a 220 V** (lineal entre 208 V y 230 V, fracción = (220−208)/(230−208) = 0.5455):
$$FLA_{220} = 20.5 + (18.6-20.5)\times 0.5455 = 19.46\ A$$

Para el LRA, WEG no publica el valor discriminado a 208 V y 230 V para esta variante tri-tensión; se aplicó el múltiplo de código de rotor bloqueado **LRC = 7.1× (Código H)**, confirmado en la ficha técnica de la versión gemela 230/460 V (130 A / 18.4 A = 7.07× ≈ 7.1×), a cada nivel de tensión antes de interpolar:
$$LRA_{208}=20.5\times7.1=145.55\ A,\quad LRA_{230}=18.6\times7.1=132.06\ A$$
$$LRA_{220} = 145.55 + (132.06-145.55)\times 0.5455 = 138.19\ A$$

Verificación de consistencia: $138.19/19.46 = 7.10\times$

---

## Ampacidad del Conductor

NEC/NOM **430.22** exige que la ampacidad del conductor del circuito derivado a un solo motor sea **≥ 125% del FLA**:

$$I_{ampacidad,min} = FLA \times 1.25 = 19.46\ A \times 1.25 = \mathbf{24.33\ A}$$

Ampacidad de referencia — NOM-001-SEDE-2012, Tabla 310-15(b)(16), THHN/THWN-2 a 75 °C, 3 conductores portadores de corriente en conducto, 30 °C ambiente. El conductor 8 AWG tiene ampacidad de **50 A**, muy por encima del mínimo requerido (50/24.33 = **205.5%** del mínimo) — el conductor existente queda con amplio margen frente al motor trifásico real (menor FLA que el valor monofásico usado por error anteriormente).

---

## Caída de Tensión (NOM-001-SEDE-2012)

Para el circuito derivado a un motor, el límite es **3%**; considerando además el alimentador, el límite total del sistema es **5%** (NOM-001-SEDE-2012, notas de la Tabla 9). **Este límite del 3%/5% está definido para condición nominal de operación de plena carga (régimen permanente), no para la corriente de arranque a rotor bloqueado** 
**Fórmula (circuito trifásico, NOM-001-SEDE-2012):**

$$\%V = \frac{I \cdot L \cdot Z_{ef} \cdot \sqrt{3}}{V_{ff}} \times 100$$

donde $Z_{ef} = R\cos\varphi + X\sin\varphi$, con FP del motor = 0.81 (dato de placa):

$$\sin\varphi = \sqrt{1-0.81^2} = 0.5864 $$

$$Z_{ef} = (2.56\ \Omega/km\times0.81) + (0.171\ \Omega/km\times0.5864) = 2.0736\ \Omega/km + 0.1003\ \Omega/km = 2.1739\ \Omega/km$$

**Cálculo a corriente nominal (FLA = 19.46 A):**

$$\%V = \frac{19.46\ A \times 0.0196\ km \times 2.1739\ \Omega/km \times \sqrt3}{220\ V} \times 100 = \mathbf{0.65\%}$$

**Caída durante el arranque (LRA = 138.19 A):**

$$\%V_{arranque} = \frac{138.19\ A \times 0.0196\ km \times 2.1739\ \Omega/km \times \sqrt3}{220\ V} \times 100 = \mathbf{4.64\%}$$

---

## Cálculo de Impedancias

### Impedancia del Transformador (referida a 220 V)

$$Z_{base} = \frac{V_{base}^2}{S_{base}} = \frac{(220\ V)^2}{500{,}000\ VA} = 0.09680\ \Omega$$

$$|Z_{trafo}| = 0.045 \times 0.09680\ \Omega = 0.004356\ \Omega$$

Con X/R = 5 (adimensional):

$$X_s = \frac{|Z_{trafo}|}{\sqrt{1+\frac{1}{(X/R)^2}}} = \frac{0.004356\ \Omega}{\sqrt{1.04}} = 0.004271\ \Omega$$

$$R_s = \frac{X_s}{X/R} = \frac{0.004271\ \Omega}{5} = 0.000854\ \Omega$$

**Z_trafo = 0.000854 + j0.004271 Ω**

### Impedancia del Cable (8 AWG, ≈19.6 m)

$$R_{cable} = 2.56\ \Omega/km \times 0.0196\ km = 0.05018\ \Omega$$
$$X_{cable} = 0.171\ \Omega/km \times 0.0196\ km = 0.00335\ \Omega$$

**Z_cable = 0.05018 + j0.00335 Ω**

### Impedancia Total hasta el Motor

$$R_{1,T} = R_s + R_{cable} = 0.000854\ \Omega + 0.05018\ \Omega = 0.05103\ \Omega$$
$$X_{1,T} = X_s + X_{cable} = 0.004271\ \Omega + 0.00335\ \Omega = 0.00762\ \Omega$$

$$|Z_{1,T}| = \sqrt{(0.05103\ \Omega)^2 + (0.00762\ \Omega)^2} = 0.05160\ \Omega$$

**Asunción de puesta a tierra:** transformador sólidamente aterrizado → **Z₀ ≈ Z₁** (aproximación de cálculo manual; el Z₀ real depende de la construcción del núcleo del transformador y puede diferir de Z₁).

---

## Corriente de Cortocircuito en Terminales del Motor

$$V_{LN} = \frac{220\ V}{\sqrt{3}} = 127.017\ V$$

### Falla Trifásica

$$I_{cc,3\phi} = \frac{V_{LN}}{|Z_{1,T}|} = \frac{127.017\ V}{0.05160\ \Omega} = \mathbf{2461.5697\ A}$$

### Falla Monofásica a Tierra (SLG)

Como Z₁ = Z₂ = Z₀ (puesta a tierra sólida, aproximación):

$$I_{cc,SLG} = \frac{3 \cdot V_{LN}}{2|Z_1|+|Z_0|} = \frac{V_{LN}}{|Z_1|} = \frac{127.017\ V}{0.05160\ \Omega} = I_{cc,3\phi} = \mathbf{2461.5697\ A}$$

---

## Corriente Nominal del Motor (FLA) y Ajustes del Relé 50/51

$$FLA = \mathbf{19.46\ A} $$
$$LRA = \mathbf{138.19\ A} \quad \$$

$$\frac{LRA}{FLA} = \frac{138.19\ A}{19.46\ A} = 7.10\times\ $$

Con **FS = 1.25** de placa (≥1.15), NEC 430.32(A)(1) exige usar **125%**:

| Ajuste | Fórmula | Resultado (primario) |
|---|---|---|
| Pickup 51 (sobrecarga) | FLA × 1.25 | **24.33 A** |
| LRA (arranque, dato de placa) | — | **138.19 A** |
| Pickup 50 (instantáneo) | LRA × 1.7 | **234.92 A** |

$$Pickup_{51} = 19.46\ A \times 1.25 = 24.33\ A$$
$$Pickup_{50} = 138.19\ A \times 1.7 = 234.92\ A$$

*(El factor 1.7 en el pickup instantáneo es un margen de coordinación de protecciones de práctica común — 1.6 a 2.0× LRA — para evitar disparo por la asimetría del arranque; no proviene de una tabla NEC/NOM específica, sino de guías de aplicación de relevadores tipo IEEE C37.96.)*

---

## Selección de Fusible de Respaldo

$$I_{min,fusible} = FLA \times 1.75 = 19.46\ A \times 1.75 = 34.06\ A$$

**Fusible seleccionado: Fusetron FRN-R-35** 

$$\frac{35\ A}{19.46\ A}=179.9\%\ FLA$$

Se aplica la **Excepción 1** de NEC 430.52(C)(1): el valor calculado (34.06 A) no coincide con un tamaño estándar de NEC 240.6(A) (...30, 35, 40...), por lo que se permite subir al siguiente tamaño estándar — 35 A. *(Verificado: NEC Tabla 430.52, fusible dual-element = 175% FLA como techo de cálculo; la Excepción 1 permite explícitamente exceder ese porcentaje al redondear al tamaño estándar superior — 179.9% .)*

**Calibración de la curva del fusible** (K es propiedad de la familia Fusetron, normalizada en múltiplos de la corriente nominal del fusible; el mismo K aplica a cualquier tamaño de la familia — modelo I²t simplificado; la curva TCC real del fabricante no es perfectamente I²t en todo el rango, pero es una aproximación razonable para verificación de coordinación):

$$X_{ref}=5.0\times I_{fusible}=5.0\times 35\ A = 175\ A,\quad t_{ref}=10\ s$$

$$K=t_{ref}\times\left(\frac{X_{ref}}{I_{fusible}}\right)^2 = 10\ s\times(5.0)^2\ \text{(adim.)} = 250\ s$$

$$t_{fusible}(I) = K \times \left(\frac{I_{fusible}}{I}\right)^2 = 250\ s \times \left(\frac{35\ A}{I}\right)^2$$

**Verificación de arranque:**
$$t_{fusible}(138.19\ A) = 250\ s\times\left(\frac{35\ A}{138.19\ A}\right)^2 = 16.04\ s$$

Margen amplio frente a un arranque típico de unos pocos segundos.

---

## Selección y Verificación del CT

$$I_{CT,min\ por\ carga} = FLA \times 1.25 = 19.46\ A\times1.25=\mathbf{24.33\ A}$$

Este valor por sí solo sugeriría un CT pequeño (p. ej. 25/5 A o 30/5 A). **Sin embargo**, el CT debe además evitar saturarse durante la falla real del sistema (2461.7 A). Para una clase de precisión 5P40 (factor límite de precisión = 40×), el CT debe cumplir:

$$I_{CT,primario} \geq \frac{I_{cc,3\phi}}{40} = \frac{2461.7\ A}{40} = 61.5\ A$$

El siguiente tamaño estándar por encima de 61.5 A es **75 A** — por eso el CT seleccionado es **75/5 A** (RTC=15), clase **5P40**.

**Nota de diseño:** esto es un compromiso típico entre precisión de medición en carga normal y precisión de protección durante fallas. Con el motor trifásico real (FLA=19.46 A), la corriente en el secundario a plena carga es de solo 19.46/15 = **1.30 A** (26% del secundario nominal de 5 A) — más baja resolución para medición de carga normal, pero necesaria para que el CT no sature durante una falla real. Si la prioridad fuera solo medición de carga, un CT 25/5 o 30/5 sería más adecuado, pero no protegería correctamente la exactitud del relevador durante la falla.

Como I_cc,3φ = 2461.7 A < 3000 A (I_ALF = 40 × 75 A = 3000 A), el CT se mantiene dentro de su clase de precisión (±5%) hasta la corriente de falla real.

**Pickup 51 en secundario:** 24.33 A / 15 (RTC) = **1.622 A secundario**
**Pickup 50 en secundario:** 234.92 A / 15 (RTC) = **15.66 A secundario**

---

## Curva de Daño del Transformador (ANSI/IEEE C57.109)

### Categoría del Transformador

| Categoría | Monofásico | Trifásico |
|---|---|---|
| I | 5–500 kVA | 15–500 kVA |
| II | 501–1667 kVA | 501–5000 kVA |
| III | 1668–10,000 kVA | 5001–30,000 kVA |
| IV | > 10,000 kVA | > 30,000 kVA |

Este transformador (**500 kVA, trifásico**) cae justo en el límite superior de **Categoría I**. *(Verificado contra IEEE C57.12.00/C57.109: los rangos de la tabla coinciden.)*

### Corriente Base del Transformador

$$I_{FLA,trafo} = \frac{S_{trafo}}{\sqrt{3}\cdot V_{LL}} = \frac{500{,}000\ VA}{\sqrt{3}\times 220\ V} = \mathbf{1312.16\ A}$$

### Corriente Máxima de Falla Posible (bornes BT del transformador, sin cable)

$$I_{max} = \frac{I_{FLA,trafo}}{Z_{pu}} = \frac{1312.16\ A}{0.045\ \text{(adim.)}} = \mathbf{29{,}159\ A} = 22.22\times I_{FLA,trafo}$$

Como esta corriente máxima (22.2×) es menor que 40×, la curva estándar de Categoría I cubre el rango completo de fallas posibles en este transformador.

### Ecuación de la Curva (Categoría I)

$$I_{pu}^2 \cdot t = K = (25\ \text{adim.})^2 \times 2\ s = 1250\ s \quad\Rightarrow\quad t(I_{pu}) = \frac{1250\ s}{I_{pu}^2}$$

Válida entre **2× y 40×** de I_FLA,trafo. *(Verificado contra literatura IEEE/SEL sobre C57.109: K=1250 y el rango de validez 2×–40× para Categoría I son correctos.)*

### Tabla de la Curva

| Múltiplo de I_FLA,trafo | Corriente (A) | Tiempo de daño (s) |
|---:|---:|---:|
| 2× | 2,624.3 | 312.50 |
| 3× | 3,936.5 | 138.89 |
| 4× | 5,248.6 | 78.12 |
| 5× | 6,560.8 | 50.00 |
| 6× | 7,873.0 | 34.72 |
| 8× | 10,497.3 | 19.53 |
| 10× | 13,121.6 | 12.50 |
| 15× | 19,682.4 | 5.56 |
| 20× | 26,243.2 | 3.12 |
| 22.22× (I_max, Z=4.5%) | 29,159.1 | 2.53 |
| 25× | 32,804.0 | 2.00 |
| 30× | 39,364.8 | 1.39 |
| 40× | 52,486.4 | 0.78 |

---

## Resumen de Ajustes Finales

| Elemento | Ajuste |
|---|---|
| Sistema | 220Y/127 V |
| Motor | WEG W22 7.5HP, **00718ET3E213T-W22 (trifásico)**, FLA=19.46A (a 220V)  |
| **Conductor** | **8 AWG**, R=2.56 Ω/km, X=0.171 Ω/km, ≈19.6 m — ampacidad 50A @75°C ✓ (205.5% ≥ requerido) |
| **Caída de tensión a FLA** | **0.65%** — cumple 3% (derivado) |
| **Caída de tensión en arranque** | **4.64%** |
| **CT** | **75/5 A** (RTC=15), clase 5P40, I_ALF=3000A — dimensionado por la falla, no por la carga |
| Fusible de respaldo | **FRN-R-35** (dual-element, Excepción 1 del 430.52) |
| Pickup 51 | 24.33 A primario → **1.622 A secundario** (RTC=15) |
| Pickup 50 | 234.92 A primario → 15.66 A secundario |
| **I_cc,3φ en el motor** | **2461.7 A**  |
| **I_cc,SLG en el motor** | **2461.7 A**  |
| Categoría transformador (ANSI/IEEE C57.109) | Categoría I (trifásico, límite 500 kVA) |
| I_FLA,trafo (base 220V) | 1312.16 A |
| I_cc,3φ en pu de I_FLA,trafo | **1.876×** — fuera del rango válido de la curva (2×–40×); el transformador no está en riesgo de daño por falla en este punto |
| t_daño extrapolado (no normativo) | 355.1 s |