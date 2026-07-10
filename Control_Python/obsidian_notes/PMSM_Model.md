# Modelo Matemático del Motor PMSM

Este documento detalla el modelo matemático continuo del Motor Síncrono de Imanes Permanentes (PMSM) implementado en [[pmsm_model.py|pmsm_model.py]]. El modelo es válido tanto para motores de imanes superficiales (SPMSM) como de imanes interiores con saliencia (IPMSM). Cabe considerar, que se controlará un motor bifásico stepper, por lo que una transformada de park e inversa de park es lo único necesario

---

## ⚡ Ecuaciones Eléctricas (Marco Giratorio $d$-$q$)

El modelo eléctrico del motor se expresa en el marco de referencia síncrono $d$-$q$ orientado al flujo del rotor:

$$v_d = R_s i_d + L_d \frac{di_d}{dt} - \omega_e L_q i_q$$

$$v_q = R_s i_q + L_q \frac{di_q}{dt} + \omega_e L_d i_d + \omega_e \psi_m$$

Despejando las derivadas de las corrientes para simulación en variables de estado:

$$\frac{di_d}{dt} = \frac{1}{L_d} \left( v_d - R_s i_d + \omega_e L_q i_q \right)$$

$$\frac{di_q}{dt} = \frac{1}{L_q} \left( v_q - R_s i_q - \omega_e L_d i_d - \omega_e \psi_m \right)$$

Donde:
* $v_d, v_q$: Tensiones de estator en el eje $d$ y $q$ $[V]$.
* $i_d, i_q$: Corrientes de estator en el eje $d$ y $q$ $[A]$.
* $R_s$: Resistencia del devanado estatórico $[\Omega]$.
* $L_d, L_q$: Inductancias de los ejes $d$ y $q$ $[H]$. Mi motor es un stepper con saliencia, por lo que $L_d \neq L_q$.
* $\psi_m$: Flujo de los imanes permanentes del rotor $[Wb]$.
* $\omega_e$: Velocidad angular eléctrica del rotor $[rad/s]$.

---

## ⚙️ Ecuaciones Mecánicas

La dinámica mecánica del rotor se modela mediante:

$$\frac{d\omega_m}{dt} = \frac{1}{J} \left( T_e - T_L - B \omega_m \right)$$

$$\frac{d\theta_m}{dt} = \omega_m$$

Donde:
* $\omega_m$: Velocidad mecánica angular del rotor $[rad/s]$.
* $\theta_m$: Posición mecánica angular del rotor $[rad]$.
* $J$: Inercia equivalente del rotor y la carga $[kg\cdot m^2]$.
* $B$: Coeficiente de fricción viscosa $[N\cdot m\cdot s/rad]$.
* $T_L$: Torque de carga (perturbación externa) $[N\cdot m]$.
* $T_e$: Torque electromagnético generado por el motor $[N\cdot m]$.

---

## ✅ Consistencia en el Torque Electromagnético ($T_e$)

Durante el proceso de migración y refactorización, se corrigieron discrepancias históricas en la definición del torque. Ahora, tanto la definición teórica (para telemetría y observadores) como la utilizada en la dinámica de integración física (`derivatives`), son matemáticamente consistentes:

$$T_e = P \cdot \psi_m i_q - P \cdot (L_d - L_q) i_d i_q$$

> [!SUCCESS]
> **Corrección Aplicada:**
> El número de pares de polos $P$ (en este caso 50) ha sido correctamente factorizado en los términos de imanes permanentes y de reluctancia cruzada en `pmsm_model.py`. La simulación de la planta reacciona ahora de forma fiel al torque modelado teóricamente, eliminando la necesidad de sobrecompensación en las ganancias del controlador por discrepancias de modelo.

---

## 🕒 Integración Numérica (RK4)

Para simular la física continua de la planta en un sistema discreto, se utiliza el método de **Runge-Kutta de 4° orden (RK4)** implementado en `step_rk4`. Los estados integrados son:

$$x = \begin{bmatrix} i_d & i_q & \omega_m & \theta_m \end{bmatrix}^T$$

La actualización en cada paso de tiempo continuo $dt_{\text{sim}}$ se realiza como:

$$k_1 = f(x(t), u(t))$$
$$k_2 = f\left(x(t) + \frac{dt_{\text{sim}}}{2}k_1, u(t)\right)$$
$$k_3 = f\left(x(t) + \frac{dt_{\text{sim}}}{2}k_2, u(t)\right)$$
$$k_4 = f(x(t) + dt_{\text{sim}}k_3, u(t))$$
$$x(t + dt_{\text{sim}}) = x(t) + \frac{dt_{\text{sim}}}{6}(k_1 + 2k_2 + 2k_3 + k_4)$$
