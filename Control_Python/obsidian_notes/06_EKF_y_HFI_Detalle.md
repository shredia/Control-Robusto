# Crítica y Detalle Técnico: EKF y HFI en el Motor Stepper

Como revisor técnico del proyecto, una de las áreas críticas que requiere profunda comprensión es cómo el **Filtro de Kalman Extendido (EKF)** se fusiona con la **Inyección de Alta Frecuencia (HFI)** para garantizar el control sin sensores a cualquier velocidad. 

Esta nota desglosa matemáticamente la implementación contenida en el archivo `observers.py`.

## 1. El Vector de Estados del EKF
El filtro estima 5 estados internos del motor en cada paso de tiempo. El vector de estados $\hat{x}$ está definido como:
$$ \hat{x} = \begin{bmatrix} i_d \\ i_q \\ \omega_m \\ \theta_m \\ T_x \end{bmatrix} $$
Donde:
- $i_d, i_q$: Corrientes en el marco rotatorio (A).
- $\omega_m$: Velocidad mecánica del rotor (rad/s).
- $\theta_m$: Posición mecánica del rotor (rad).
- $T_x$: Torque de perturbación de carga extendido (N*m).

## 2. Etapa de Predicción (EKF)
Usando los voltajes aplicados $v_d, v_q$ y el modelo no lineal del motor stepper, predecimos cómo evolucionarán los estados. En `observers.py`, esto se hace mediante integración de Euler explícita:

$$ \frac{di_d}{dt} = \frac{1}{L_d} \left( v_d - R_s i_d + N_r \omega_m L_q i_q \right) $$
$$ \frac{di_q}{dt} = \frac{1}{L_q} \left( v_q - R_s i_q - N_r \omega_m L_d i_d - K_e \omega_m \right) $$
$$ \frac{d\omega_m}{dt} = \frac{1}{J} \left( K_t i_q - N_r(L_d - L_q)i_d i_q - T_x - B \omega_m \right) $$
$$ \frac{d\theta_m}{dt} = \omega_m $$
$$ \frac{dT_x}{dt} = 0 \quad \text{(Se asume constante entre muestreos)} $$

**Matriz Jacobiana ($A$ o $F_k$):**
El EKF requiere linealizar este modelo en cada instante. La matriz Jacobiana evaluada en el estado actual actualiza la matriz de covarianza de error $P$:
$$ P_{pred} = \Phi P_{k-1} \Phi^T + Q_k $$
*(Donde $\Phi \approx I + A T_s$ y $Q_k$ es el ruido de proceso).*

## 3. Implementación de Inyección de Alta Frecuencia (HFI)
El EKF por sí solo pierde observabilidad a $\omega_m \approx 0$ porque la BEMF (Fuerza Contraelectromotriz) se vuelve cero. HFI resuelve esto explotando la **saliencia magnética** ($L_d \neq L_q$).

### A. Filtrado y Extracción
Se inyecta un voltaje de alta frecuencia $V_{hfi}$ a $f_h = 1000$ Hz. En el software, la corriente resultante $I_{\alpha, \beta}$ se pasa por un **Filtro Pasa-Banda (BPF) Butterworth** centrado en $f_h$ para extraer únicamente el rizado inducido ($I_{\alpha, hfi}, I_{\beta, hfi}$).

### B. Demodulación Ortogonal (I/Q)
El rizado de corriente se rota de vuelta a $d-q$ estimado. La componente cuadratura ($I_{q, hfi}$) contiene el error de alineación. Se multiplica por el seno de la frecuencia de inyección y se filtra (Filtro Pasa-Bajos) para extraer el error base:
$$ \varepsilon_f = \text{LPF}( I_{q, hfi} \cdot \sin(\omega_h t) ) $$

### C. Normalización y Cálculo de Ángulo
Para que el error no dependa de la amplitud absoluta, se normaliza por el factor de saliencia:
$$ \text{Saliencia} = \frac{L_q - L_d}{0.5(L_q + L_d)} $$
$$ \varepsilon_{norm} = \frac{\varepsilon_f}{|\text{Saliencia}|} $$
Este error normalizado se utiliza para corregir el ángulo eléctrico estimado por el EKF ($\hat{\theta}_{ekf}$):
$$ \theta_{hfi} = \hat{\theta}_{ekf} + K_{pll} \cdot \varepsilon_{norm} $$

## 4. Etapa de Corrección (EKF)
El EKF de este proyecto usa una **corrección en dos etapas** para asimilar toda la información disponible.

### Corrección 1: Mediante Corrientes Fundamentales
Se usan las corrientes medidas en las fases (limpias de ruido HFI): $z_1 = [i_{\alpha}, i_{\beta}]^T$.
La matriz Jacobiana de observación $C_k$ se calcula rotando las derivadas parciales de $i_d, i_q$ mediante la matriz de Park.
$$ \hat{x}_{corr1} = \hat{x}_{pred} + K_{curr} (z_1 - h_1(\hat{x}_{pred})) $$

### Corrección 2: Mediante Ángulo HFI (Adaptativo)
Si el HFI está activo (la amplitud del ruido es mayor al umbral de fondo), se usa $\theta_{hfi}$ como un sensor virtual directo del ángulo eléctrico.
La matriz de observación es simplemente $H_{hfi} = [0, 0, 0, N_r, 0]$.
$$ \hat{x}_{corr2} = \hat{x}_{corr1} + K_{hfi} (\theta_{hfi} - N_r \hat{\theta}_{m, corr1}) $$

Esta arquitectura dual garantiza que el EKF tenga alta precisión dinámica a altas velocidades (Corrección 1), pero nunca pierda la sincronía a cero velocidad o bajo alto torque (Corrección 2 con HFI).
