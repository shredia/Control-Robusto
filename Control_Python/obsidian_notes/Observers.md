# Estimación de Perturbaciones y Observadores

Este documento detalla los algoritmos de estimación de torque perturbador externo.

> [!NOTE]
> **Estado del SMO:**
> El Observador de Modos Deslizantes (SMO) con PLL para estimación sensorless ha sido **desactivado temporalmente** a petición del usuario. La simulación utiliza el sensor físico directo del motor para la realimentación de posición y velocidad.

---

## ⚖️ Observador de Torque de Carga (LOB)

El LOB implementado en [[observers.py|observers.py]] estima el torque de carga externo $T_L$ basándose en la velocidad mecánica medida $\omega_m$ y en el torque electromagnético calculado $T_e$:

### 1. Formulación Dinámica
El error de velocidad es:

$$e_\omega = \omega_m - \hat{\omega}_m$$

La dinámica de estimación se rige por:

$$\frac{d\hat{\omega}_m}{dt} = \frac{1}{J_{\text{var}}} \left( T_e - \hat{T}_L - B_{\text{var}} \hat{\omega}_m \right) + L_1 e_\omega$$

$$\frac{d\hat{T}_L}{dt} = -L_2 e_\omega$$

Donde $J_{\text{var}}$ y $B_{\text{var}}$ son los parámetros nominales del modelo del controlador.

### 2. Discretización en Código (Euler hacia adelante)
```python
err_omega = omega_m - self.omega_m_est

d_omega = (torque_e - self.T_L_est - self.B * self.omega_m_est) / self.J + self.L1 * err_omega
d_TL = -self.L2 * err_omega

self.omega_m_est += self.dt * d_omega
self.T_L_est += self.dt * d_TL
```

### 3. Ganancias por Polo Doble Real
Las ganancias $L_1$ y $L_2$ se calculan para imponer un polo doble en $-p_{\text{obs}}$ (con $p_{\text{obs}} = 50.0\text{ rad/s}$):

$$L_1 = 2 p_{\text{obs}} - \frac{B_{\text{var}}}{J_{\text{var}}}$$

$$L_2 = J_{\text{var}} \cdot p_{\text{obs}}^2$$

---

## 🌪️ Observador de Perturbaciones (DOB - Matlab/Simulink)

El script de inicialización `InitFcn.m` revela el diseño de un **Observador de Perturbaciones de 3er Orden** para su uso en la simulación de Simulink, basado en el polinomio de Hurwitz:

### 1. Ecuaciones de Diseño en MATLAB
Se sintoniza con una frecuencia natural deseada $f_{\text{bw,DO}} = 25\text{ Hz}$ ($\omega_{\text{do}} = 2\pi f_{\text{bw,DO}} \approx 157.08\text{ rad/s}$):

$$\hat{a} = \frac{B_{\text{var}}}{J_{\text{var}}}$$

Las ganancias de estimación del DOB se calculan como:

$$\beta_0 = 3\omega_{\text{do}} - \hat{a}$$

$$\beta_1 = 3\omega_{\text{do}}^2 - \beta_0 \hat{a}$$

$$\beta_2 = \omega_{\text{do}}^3 \frac{J_{\text{var}}}{P}$$

**Valores Numéricos (Calculados en `parameters.py`):**
Sustituyendo $\omega_{\text{do}} \approx 157.08\text{ rad/s}$, $J_{\text{var}} = 54\times 10^{-7}$ y $B_{\text{var}} = 1\times 10^{-4}$:
* $\hat{a} \approx 18.52$
* $\beta_0 \approx 452.7$
* $\beta_1 \approx 65646$
* $\beta_2 \approx 418$

Este estimador de 3er orden está directamente unificado con los estados mecánicos del EKF, estimando el rizado y permitiendo un Feedforward de perturbaciones robusto.
