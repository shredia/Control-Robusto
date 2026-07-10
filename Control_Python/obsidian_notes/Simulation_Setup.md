# Configuración de Simulación y Resultados

Este documento describe la arquitectura de simulación multitasa en tiempo mixto (continuo-discreto), el desajuste de parámetros simulado, y la visualización interactiva de resultados en la carpeta de salidas.

---

## 🕒 Esquema Multitasa y Tiempos de Muestreo

La simulación coordina cuatro escalas de tiempo diferentes en [[simulation.py|simulation.py]], utilizando retenedores de orden cero (ZOH) para las consignas discretas intermedias:

1. **Planta Continua ($dt_{\text{sim}} = 10\ \mu\text{s}$):**
   * Integración numérica RK4 de las ecuaciones dinámicas eléctricas y mecánicas del motor y del puente H.
2. **Lazo de Corriente ($dt_{\text{current}} = 100\ \mu\text{s} \implies 10\text{ kHz}$):**
   * Lectura de corrientes y actualización de voltajes $v_d, v_q$ con desacoplamiento.
3. **Lazo de Velocidad ($dt_{\text{speed}} = 1\text{ ms} \implies 1\text{ kHz}$):**
   * Lectura de velocidad y actualización de consigna de corriente $i_{q,\text{ref}}$.
4. **Lazo de Posición ($dt_{\text{pos}} = 10\text{ ms} \implies 100\text{ Hz}$):**
   * Comparación de posición y actualización de consigna de velocidad $\omega_{\text{ref}}$.

---

## ⚖️ Desajuste de Parámetros (Plant vs Model Mismatch)

Para evaluar la robustez del control, se ha introducido un desajuste del **100% en la inercia y fricción** del sistema:

* **Planta Física Real (Simulada):**
   * $J_{\text{real}} = 1.08 \times 10^{-5}\text{ kg}\cdot\text{m}^2$ (Motor + Carga real)
   * $B_{\text{real}} = 2.0 \times 10^{-4}\text{ N}\cdot\text{m}\cdot\text{s/rad}$
* **Controlador y Observadores (Modelo Estimado):**
   * $J_{\text{var}} = 5.4 \times 10^{-6}\text{ kg}\cdot\text{m}^2$
   * $B_{\text{var}} = 1.0 \times 10^{-4}\text{ N}\cdot\text{m}\cdot\text{s/rad}$

Este desajuste obliga a los lazos cerrados a compensar los errores de modelado, verificando la estabilidad del control robusto en cascada.

---

## 📂 Directorio de Resultados (`/outputs`)

Todas las salidas generadas se almacenan automáticamente en la carpeta:
👉 `d:\Utal\Memoria_Clinostato\Memoria-clinostato\Control Robusto\Control_Python\outputs`

### 1. Gráfica Estática (`sim_stepper_sensado.png`)
Muestra la respuesta temporal graficada en Matplotlib con subplots de posición, velocidad, corrientes dq y corrientes de fase estatóricas A y B.

### 2. Visor Interactivo (Scope) (`scope_stepper.html`)
Un archivo HTML autocontenido generado con **Plotly** que se comporta como un **Scope de Simulink**:
* **Zoom Dinámico:** Permite seleccionar cualquier rango de tiempo con el cursor para inspeccionar transitorios y rizados finos.
* **Leyenda Interactiva:** Haciendo clic en cualquier elemento de la leyenda se puede ocultar/mostrar curvas de referencia o medición en tiempo real.
* **Hover de Datos:** Muestra las lecturas de todas las señales simultáneamente para el instante en que se posiciona el mouse.

---

## 📈 Trayectoria de Prueba

* **Referencia de Posición:** Rampa lineal $\theta_{\text{ref}}(t) = 10 \cdot t$.
* **Velocidad de Régimen:** $10\text{ rad/s}$.
* **Fuerza Externa:** $T_L = 0.0\text{ N}\cdot\text{m}$ (motor en vacío).
* El motor es capaz de seguir la rampa con un error de velocidad de prácticamente cero y un desfase permanente de posición (velocity lag) de $\approx 12.8\text{ rad}$, estabilizando la corriente de par $i_q$ en un rizado centrado en la corriente media requerida para vencer el cogging.
