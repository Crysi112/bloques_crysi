# Estudio de Protecciones Complementarias (Voltaje, Frecuencia, Distancia y Térmicas)
### Alimentador IEEE 13 Nodos — Funciones Avanzadas de Subestación y Cargas Críticas

---

## 6. Protección de Línea de Subtransmisión (Acometida 115 kV)

### 6.1 Protección de Distancia (ANSI 21/21N) — (Ref. IEEE C37.113)
Aplicada en el relé de entrada de 115 kV para proteger la línea contra fallas entre fases y a tierra, midiendo la impedancia $Z_{medida} = V/I$.

**Parámetros asumidos para la línea de 115 kV:**
* Impedancia de secuencia positiva de la línea ($Z_{1L}$): $0.5\ \Omega/km$
* Longitud hipotética de la línea fuente: $20\ km \rightarrow Z_{linea} = \mathbf{10\ \Omega\ primarios}$
* RTC = 400/5 (80), PTC = 115kV/115V (1000)
* Razón de impedancias $Z_{ratio} = PTC / RTC = 1000 / 80 = \mathbf{12.5}$

**Ajustes de Zonas Mho (Valores Secundarios):**
* **Zona 1 (Instantánea - 50 a 80% de la línea):**
  $$Z_{1} = 0.80 \times \left(\frac{10\ \Omega}{12.5}\right) = \mathbf{0.64\ \Omega\ secundarios}$$
  $$t_{Z1} = \mathbf{0.0\ s}\ (\text{Instantáneo})$$

* **Zona 2 (Respaldo Trafo y barra local - 120% de la línea):**
  $$Z_{2} = 1.20 \times \left(\frac{10\ \Omega}{12.5}\right) = \mathbf{0.96\ \Omega\ secundarios}$$
  $$t_{Z2} = \mathbf{0.3\ s\ a\ 0.5\ s}$$

> **Notas del Alumno sobre IEEE C37.113 (Distancia):**
> *[Investigar por qué la Zona 1 nunca se ajusta al 100% de la línea y el efecto del "Infeed" en la medición de impedancia...]*

---

## 7. Protección de Estabilidad de Barra (Subestación 650 y Nodo 671)

### 7.1 Bajo y Sobre Voltaje (ANSI 27/59) — (Ref. ANSI C84.1 / IEEE 1547)
Ajustado para proteger equipos electrónicos, evitar sobrecalentamiento en motores por caída de tensión y detectar islas no intencionales.

**Ajustes (Base $V_{nom} = 4.16\ kV$ o $1.0\ pu$):**
* **59 - Sobrevoltaje (Nivel 1 - Alarma):**
  $$V_{59-1} = 1.10 \times V_{nom} = \mathbf{4.57\ kV} \quad (t = 10.0\ s)$$
* **59 - Sobrevoltaje (Nivel 2 - Disparo):**
  $$V_{59-2} = 1.20 \times V_{nom} = \mathbf{4.99\ kV} \quad (t = 0.16\ s)$$

* **27 - Bajo Voltaje (Nivel 1 - Motores Nodo 671):**
  $$V_{27-1} = 0.90 \times V_{nom} = \mathbf{3.74\ kV} \quad (t = 5.0\ s)$$
* **27 - Bajo Voltaje (Nivel 2 - Aislamiento Falla):**
  $$V_{27-2} = 0.80 \times V_{nom} = \mathbf{3.32\ kV} \quad (t = 1.0\ s)$$

> **Notas del Alumno sobre ANSI C84.1 (Límites de Voltaje):**
> *[Justificar los retardos de tiempo en las funciones 27/59 respecto a la duración típica de un cortocircuito (ride-through)...]*

---

## 8. Protección de Frecuencia y Deslastre de Carga

### 8.1 Baja y Sobre Frecuencia (ANSI 81U/81O) — (Ref. NERC PRC-006 / IEEE C37.117)
Utilizado a nivel subestación para UFLS (Underfrequency Load Shedding) y proteger a las turbinas/generadores del sistema interconectado.

**Ajustes (Frecuencia Nominal $f_n = 60.0\ Hz$):**
* **81O - Sobre Frecuencia (Pérdida súbita de carga grande):**
  $$f_{81O} = \mathbf{60.5\ Hz} \quad (t = 2.0\ s)$$

* **81U - Baja Frecuencia (Esquema de Deslastre - Nivel 1):**
  $$f_{81U-1} = \mathbf{59.5\ Hz} \quad (t = 10.0\ s\ \text{- Alarma})$$
* **81U - Baja Frecuencia (Esquema de Deslastre - Nivel 2, Desconexión Industrial 671):**
  $$f_{81U-2} = \mathbf{59.0\ Hz} \quad (t = 0.2\ s\ \text{- Disparo})$$
* **81U - Baja Frecuencia (Esquema de Deslastre - Nivel 3, Desconexión Total):**
  $$f_{81U-3} = \mathbf{58.5\ Hz} \quad (t = 0.1\ s\ \text{- Disparo})$$

> **Notas del Alumno sobre NERC PRC-006 / IEEE C37.117:**
> *[Explicar qué es un esquema UFLS y por qué la industria en el nodo 671 se debe desconectar antes que los circuitos residenciales...]*

---

## 9. Protecciones Físicas del Transformador de Subestación (5 MVA)

### 9.1 Protección Térmica de Devanados y Aceite (ANSI 49 / 26) — (Ref. IEEE C57.91)
Monitorea la temperatura real del transformador basada en sondas PT100 y modelos de imagen térmica de la corriente pasante.

**Ajustes (Aislamiento Clase A, $65^\circ C$ de elevación):**
* Temperatura ambiente base: $30^\circ C$
* Temperatura máxima Top-Oil (Aceite Superior):
  $$T_{oil, alarma} = \mathbf{90^\circ C} \quad | \quad T_{oil, disparo} = \mathbf{105^\circ C}$$
* Temperatura de Punto Caliente del Devanado (Winding Hot-Spot):
  $$T_{hs, alarma} = \mathbf{110^\circ C} \quad | \quad T_{hs, disparo} = \mathbf{120^\circ C}$$

### 9.2 Protección Mecánica de Gases / Presión (ANSI 63 / 71) — (Ref. IEEE C57.104)
Operación autónoma sin retardos, conectada directamente a la bobina de disparo del interruptor (86 Lockout).

* **Relevador Buchholz (Trafo con tanque conservador):**
  * Nivel 1: Acumulación lenta de gases combustibles $\rightarrow$ **Alarma**
  * Nivel 2: Flujo violento de aceite (arco interno) $\rightarrow$ **Disparo Instantáneo**
* **Válvula de Sobrepresión Súbita (Sudden Pressure Relay - 63):**
  * Tasa de incremento de presión: $\Delta P / \Delta t > \text{Límite de calibración mecánica} \rightarrow$ **Disparo Instantáneo**

> **Notas del Alumno sobre IEEE C57.91 y ANSI 63:**
> *[Investigar por qué una protección diferencial (87T) podría no detectar fallas incipientes entre espiras, y por qué el relé Buchholz o 63 es vital en este caso...]*
