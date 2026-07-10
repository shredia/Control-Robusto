# Control Vectorial (FOC) en Cascada de 3 Lazos

Este documento describe la estructura y ecuaciones de control vectorial en cascada de tres niveles implementada en [[controllers.py|controllers.py]] para el motor stepper bifásico.

---

## 🎛️ Estructura del Control en Cascada Multitasa

El control de movimiento se organiza en tres lazos concéntricos discretos con diferentes períodos de muestreo:

```
          dt_pos = 10 ms                  dt_speed = 1 ms                 dt_current = 0.1 ms
      +--------------------+          +--------------------+          +--------------------+
θ_ref | Lazo de Posición   |  ω_ref   | Lazo de Velocidad  |  i_q_ref | Lazo de Corriente  | v_dq_cmd
----->| PD                 |--------->| PI                 |--------->| PI + Desacople     |--------> Inverter
      | (Kp=1.0, Kd=0.05)  |          | (Kp_w, Ki_w)       |          | (Saliencia)        |
      +--------------------+          +--------------------+          +--------------------+
               ^                               ^                               ^
               | (θ_m)                         | (ω_m)                         | (i_d, i_q)
               +-------------------------------+-------------------------------+ (Feedback Real)
```

---

## ⚡ Formulación del Controlador PID Discreto (`DiscretePID`)

Se utiliza una estructura unificada de regulador lineal discreto utilizando el método de integración de **Euler hacia atrás** y una **derivada filtrada** para atenuar el ruido de alta frecuencia:

### 1. Ecuaciones del Algoritmo
* **Término Proporcional:**
  $$u_p(k) = K_p \cdot e(k)$$

* **Término Integral (Euler hacia atrás):**
  $$u_i(k) = u_i(k-1) + K_i \cdot dt \cdot e(k)$$

* **Término Derivativo Filtrado:**
  Para un filtro de derivada con coeficiente $N_{\text{filter}}$ (típicamente $10$), la constante de tiempo del filtro es $T_f = K_d / (K_p \cdot N_{\text{filter}})$. La derivada se discretiza como:
  $$u_d(k) = \frac{T_f}{T_f + dt} u_d(k-1) + \frac{K_d}{T_f + dt} \left( e(k) - e(k-1) \right)$$

* **Salida Bruta y Clamping (Anti-Windup):**
  $$u_{\text{raw}}(k) = u_p(k) + u_i(k) + u_d(k)$$
  $$u_{\text{out}}(k) = \text{clip}\left( u_{\text{raw}}(k), u_{\min}, u_{\max} \right)$$
  $$\text{Si } (u_{\text{raw}} \text{ saturado AND } \text{signo}(e(k)) == \text{signo}(u_{\text{raw}})): \quad u_i(k) = u_i(k-1)$$

---

## 📐 Diseño y Sintonía de los 3 Lazos (Asignación de Polos)

Las ganancias del control se calculan dinámicamente en [[parameters.py|parameters.py]] empleando el diseño por asignación de polos.

### 1. Lazo de Corriente Interno ($dt_{\text{current}} = 100\ \mu\text{s}$)
Se sintonizan lazos independientes para los ejes $d$ y $q$ debido a la saliencia del motor ($L_d \neq L_q$), logrando un ancho de banda de $f_{bw} = 500\text{ Hz}$ ($\omega_d = 2\pi f_{bw} \approx 3141.6\text{ rad/s}$) y amortiguamiento $\zeta = 0.707$:

* **Eje $d$ (Control de Flujo, $i_{d,\text{ref}} = 0$):**
  $$K_{p,d} = 2 \cdot \zeta_d \cdot \omega_d \cdot L_d - R_s = 2 \cdot 0.707 \cdot 3141.6 \cdot 3.66\times 10^{-3} - 2.5 \approx 13.76$$
  $$K_{i,d} = \omega_d^2 \cdot L_d = (3141.6)^2 \cdot 3.66\times 10^{-3} \approx 36123$$

* **Eje $q$ (Control de Torque/Aceleración):**
  $$K_{p,q} = 2 \cdot \zeta_q \cdot \omega_d \cdot L_q - R_s = 2 \cdot 0.707 \cdot 3141.6 \cdot 9.61\times 10^{-3} - 2.5 \approx 40.19$$
  $$K_{i,q} = \omega_d^2 \cdot L_q = (3141.6)^2 \cdot 9.61\times 10^{-3} \approx 94848$$

* **Desacoplamiento Cruzado y Compensación:**
  $$v_{d,\text{cmd}} = v_{d,PI} - \omega_e L_q i_q$$
  $$v_{q,\text{cmd}} = v_{q,PI} + \omega_e \left( L_d i_d + \psi_m \right)$$

### 2. Lazo de Velocidad Intermedio ($dt_{\text{speed}} = 1\text{ ms}$)
Sintonizado para un ancho de banda de $f_{bw} = 50\text{ Hz}$ ($\omega_{n,w} = 2\pi \cdot 50 \approx 314.16\text{ rad/s}$) y amortiguamiento crítico $\zeta_w = 1.0$, utilizando los parámetros nominales del modelo ($J_{\text{var}}, B_{\text{var}}$):

$$K_{p,\omega} = \frac{2 \cdot \zeta_w \cdot \omega_{n,w} \cdot J_{\text{var}} - B_{\text{var}}}{K_e} \approx 0.0202$$

$$K_{i,\omega} = \frac{\omega_{n,w}^2 \cdot J_{\text{var}}}{K_e} \approx 3.268$$

### 3. Lazo de Posición Externo ($dt_{\text{pos}} = 10\text{ ms}$)
Regulado por un controlador PD clásico diseñado para atenuar perturbaciones espaciales como el torque de detent/cogging:

$$K_{p,\text{pos}} = 1.0$$
$$K_{i,\text{pos}} = 0.0$$
$$K_{d,\text{pos}} = 0.05$$
