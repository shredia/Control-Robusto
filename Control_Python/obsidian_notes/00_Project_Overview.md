# Estructura del Proyecto: Control Cascaded de Stepper Bifásico

Este directorio contiene las notas estructuradas en formato Markdown compatibles con **Obsidian** para documentar el modelo físico, control y simulaciones de un Motor Paso a Paso Bifásico (Stepper Motor) con saliencia y torque de detent.

---

## 🗺️ Mapa de Contenidos (MOC)

### 1. Modelo de Planta
* **[[PMSM_Model]]**: Ecuaciones eléctricas $d$-$q$ y mecánicas del motor stepper bifásico. Modelo físico del torque de reluctancia con signo de saliencia y dinámica de cogging/detent.
* **[[Inverter]]**: Modelo del inversor bifásico compuesto por dos puentes H independientes, soportando modulación promedio y conmutación bipolar.

### 2. Algoritmos de Control
* **[[FOC_Control]]**: Teoría del control vectorial cascada de tres niveles (**Posición $\rightarrow$ Velocidad $\rightarrow$ Corrientes $I_d, I_q$**). Sintonización de reguladores PID/PD discretos y compensación de acoplamiento.

### 3. Estimadores y Observadores
* **[[Observers]]**: Observador de Torque de Carga (LOB) para estimar perturbaciones mecánicas de carga en el eje del motor. (Nota: Observador SMO desactivado temporalmente).

### 4. Simulación e Integración
* **[[Simulation_Setup]]**: Estructura de simulación multitasa con desajuste de parámetros (Plant vs Model). Carpeta de resultados y visor dinámico interactivo HTML (tipo Simulink Scope).

---

## 📁 Archivos de Código Fuente en el Workspace

* [[parameters.py|parameters.py]]: Archivo central con los parámetros nominales/reales del motor, condiciones iniciales, tiempos de muestreo multitasa y sintonía automática de ganancias.
* [[controllers.py|controllers.py]]: Clase `CascadedStepperController` (control en cascada de 3 niveles) y `DiscretePID` (Euler hacia atrás con filtro de derivada y anti-windup clamping).
* [[inverter.py|inverter.py]]: Clase `TwoPhaseInverter` que modela la entrega de voltaje a las fases estatóricas A y B.
* [[pmsm_model.py|pmsm_model.py]]: Clase `StepperMotor2Phase` que integra la física continua del motor con torque de saliencia, fricción y cogging mediante RK4.
* [[simulation.py|simulation.py]]: Clase `PMSMSimulation` que coordina los lazos continuos e independientes lazos discretos (multitasa ZOH).
* [[main.py|main.py]]: Punto de entrada que define la rampa de referencia para $10\text{ rad/s}$ y ejecuta la simulación guardando gráficos en `/outputs`.
