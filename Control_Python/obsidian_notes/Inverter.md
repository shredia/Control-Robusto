# Modelo del Inversor Bifásico (Puentes H)

Este documento detalla la implementación del inversor bifásico en [[inverter.py|inverter.py]], el cual sustituye a la topología trifásica con SVPWM.

---

## 🛠️ Configuración de Puentes H Completos

El motor stepper bifásico se alimenta de dos puentes H independientes, uno por fase (A y B). Cada puente H es alimentado por un voltaje de bus DC de $V_{dc} = 24.0\text{ V}$ y puede aplicar un rango dinámico completo de tensión de $[-V_{dc}, V_{dc}]$.

---

## 🎛️ Modos de Operación

### 1. Modelo Promedio (`average`)
En este modo (utilizado por defecto para simulación de alta velocidad de cálculo), las tensiones de referencia calculadas por el controlador $v_{a,\text{cmd}}, v_{b,\text{cmd}}$ se limitan directamente al rango dinámico disponible:

$$v_{a} = \text{clip}\left(v_{a,\text{cmd}}, -V_{dc}, V_{dc}\right)$$

$$v_{b} = \text{clip}\left(v_{b,\text{cmd}}, -V_{dc}, V_{dc}\right)$$

Este modelo elimina el rizado de alta frecuencia y es ideal para validar las dinámicas de control macroscópicas.

### 2. Modelo de Conmutación Bipolar (`switching`)
Este modo simula la conmutación a alta frecuencia ($f_{sw} = 10\text{ kHz}$) de los transistores del puente H mediante comparación con una portadora triangular.

* **Ciclos de Trabajo (Duty Cycles):**
  Las consignas se mapean linealmente del rango $[-V_{dc}, V_{dc}]$ al intervalo adimensional $[0, 1]$:

  $$d_a = \frac{v_{a,\text{cmd,sat}}}{2 V_{dc}} + 0.5$$

  $$d_b = \frac{v_{b,\text{cmd,sat}}}{2 V_{dc}} + 0.5$$

* **Señales de Conmutación Bipolar ($S_a, S_b \in \{-1, 1\}$):**
  Se define una señal triangular portadora simétrica $C(t) \in [0, 1]$ de frecuencia $f_{sw}$. La conmutación de los brazos del puente H se establece como:

  $$S_i = \begin{cases} 1 & \text{si } d_i > C(t) \\ -1 & \text{si } d_i \leq C(t) \end{cases} \quad \text{para } i \in \{a, b\}$$

* **Voltaje de Fase Salida:**
  Las tensiones aplicadas directamente sobre los devanados A y B son:

  $$v_a = S_a \cdot V_{dc}$$

  $$v_b = S_b \cdot V_{dc}$$

Este modo introduce el chattering de conmutación de alta frecuencia, ideal para verificar la robustez de los filtros y de estimadores del sistema.
