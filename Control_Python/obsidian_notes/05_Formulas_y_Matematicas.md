# Fórmulas y Matemáticas del Control

Esta nota detalla las ecuaciones matemáticas clave utilizadas en la simulación, desde el diseño de los controladores (Asignación de polos) hasta los límites de saturación física.

## 1. Frecuencias de Muestreo (Separación de Lazos en Cascada)
Para que los lazos de control externos (ej. Posición) funcionen correctamente sin interferir con los lazos internos (ej. Corriente), es una regla de diseño separarlos por un factor de 10 en su frecuencia de operación.
- **Frecuencia PWM ($f_{carrier}$)**: 20 kHz
- **Lazo de Corriente ($f_{current}$)**: $f_{carrier} / 2 = 10$ kHz ($T_s = 100 \mu s$)
- **Lazo de Velocidad ($f_{wm}$)**: $f_{current} / 10 = 1$ kHz ($T_s = 1$ ms)
- **Lazo de Posición ($f_{pos}$)**: $f_{wm} / 10 = 100$ Hz ($T_s = 10$ ms)

## 2. Saturación de Voltaje (Inversor)
El inversor bifásico que alimenta al motor stepper está limitado físicamente por el voltaje del bus de corriente continua ($V_{dc}$). 
En el simulador (`inverter.py`), los voltajes de referencia de fase A y B son saturados antes de modular los ciclos de trabajo (Duty Cycles).
- **Voltaje de Bus:** $V_{dc} = 24.0$ V
- **Saturación:**
  $$v_{a, sat} = \max(-V_{dc}, \min(V_{dc}, v_{a, ref}))$$
  $$v_{b, sat} = \max(-V_{dc}, \min(V_{dc}, v_{b, ref}))$$
- **Cálculo del Ciclo de Trabajo Bipolar ($d \in [0, 1]$):**
  $$d_a = \frac{v_{a, sat}}{2 V_{dc}} + 0.5$$

## 3. Saturación de Corriente
El motor tiene restricciones térmicas y magnéticas. 
- **Corriente Nominal ($I_{nom}$):** 1.7 A
- **Límite Absoluto ($I_{max}$):** Se define como $I_{nom} \sqrt{2} \approx 2.4$ A.
En el controlador de velocidad (`simulation.py`), la salida de torque deseado, que se traduce a una referencia de corriente $i_q^*$, está limitada a $\pm I_{max}$ para no sobrecalentar las bobinas.

## 4. Cálculo de Ganancias PI (Asignación de Polos)
Las ganancias de los controladores PI se sintonizan en `parameters.py` imponiendo un ancho de banda ($BW$) y un coeficiente de amortiguamiento ($\xi$).

### A. Lazo de Corriente (Ejes $d$ y $q$)
Dado el modelo eléctrico de la bobina (RL), imponemos $\omega_d = \omega_q = 2\pi(500 \text{ Hz})$ y $\xi = 0.707$:
- **Eje Directo ($d$):**
  $$K_{p,d} = 2 \xi_d \omega_d L_d - R_s$$
  $$K_{i,d} = \omega_d^2 L_d$$
- **Eje Cuadratura ($q$):**
  $$K_{p,q} = 2 \xi_q \omega_q L_q - R_s$$
  $$K_{i,q} = \omega_q^2 L_q$$

### B. Lazo de Velocidad
Utilizamos la ecuación mecánica del motor con la inercia ($J$) y fricción viscosa ($B$). Imponemos $\omega_{n,w} = 2\pi(50 \text{ Hz})$ y amortiguamiento crítico ($\xi_w = 1.0$).
Primero, necesitamos la constante electromagnética $K_e$:
$$K_e = \frac{T_{hold}}{\sqrt{2} I_{nom}}$$
Las ganancias PI resultantes son:
$$K_{p,w} = \frac{2 \xi_w \omega_{n,w} J - B}{K_e}$$
$$K_{i,w} = \frac{\omega_{n,w}^2 J}{K_e}$$

### C. Lazo de Posición
El controlador de posición es un PD (Proporcional - Derivativo). Al ser el lazo más externo, es el más lento.
- $K_{p, pos} = 1.0$ (Rigidez de la posición)
- $K_{d, pos} = 0.05$ (Amortiguamiento para suavizar el comportamiento cerca de la meta y contrarrestar el par de retención o cogging).

## 5. Filtro de Observación de Perturbaciones (DOB)
Para estimar la perturbación de carga (el torque variable del clinostato), el EKF está asistido por un DOB de tercer orden sintonizado a $f_{bw, DO} = 50$ Hz ($\omega_{n,do}$).
Teniendo $\hat{BJ} = B / J$:
$$B_0 = 3 \omega_{n,do} - \hat{BJ}$$
$$B_1 = 3 \omega_{n,do}^2 - B_0 \hat{BJ}$$
$$B_2 = \omega_{n,do}^3 \cdot \frac{J}{P}$$
Estas ganancias ($\beta_0, \beta_1, \beta_2$) alimentan la dinámica interna del observador.

## 6. Inyección de Alta Frecuencia (HFI)
A bajas velocidades, inyectamos una pequeña señal de voltaje sinusoidal para rastrear la saliencia magnética.
- **Amplitud Inyectada:** $V_{hfi} = 2.0$ V
- **Frecuencia Inyectada:** $f_h = 1000$ Hz
Esta señal genera un rizado de alta frecuencia en la corriente, el cual es demodulado a través de un Filtro Pasa-Banda (BPF) Butterworth centrado en la misma frecuencia para extraer el ángulo eléctrico estimado $\hat{\theta}_e$.
