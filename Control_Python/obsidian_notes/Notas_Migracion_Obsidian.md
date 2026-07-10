# Migración de Control Robusto (Simulink/MATLAB a Python)

Este documento resume de forma detallada toda la arquitectura, estrategias de control y módulos implementados durante la migración del esquema de control de motor Stepper Bifásico desde MATLAB/Simulink hacia Python.

## 1. Arquitectura de Simulación Multitasa (Tiempo Mixto)
El núcleo de la simulación en `simulation.py` (`PMSMSimulation`) replica el comportamiento de pasos discretos (Zero-Order Hold) de los microcontroladores y la dinámica continua de la planta:

- **Planta Continua (RK4)**: `dt_sim = 10 μs`.
- **Lazo de Corriente e Inversor**: `dt_current = 100 μs` (10 kHz).
- **Observador (EKF + DOB)**: `dt_observer = 100 μs` (10 kHz).
- **Lazo de Velocidad**: `dt_speed = 1000 μs` (1 kHz).
- **Lazo de Posición**: `dt_pos = 10000 μs` (100 Hz).

La simulación coordina estos lazos ejecutándolos secuencialmente para garantizar que no existan retrasos artificiales entre la captura de estado y la acción de control.

---

## 2. Bloques de Control en Cascada
La estructura de control se modularizó en clases de Python (`control_blocks/`) manteniendo fidelidad 1 a 1 con el código de las S-Functions originales de MATLAB:

### A. Lazo de Posición (`PID_posicion.py`)
- Recibe la referencia de posición objetivo.
- Genera una trayectoria suavizada usando un generador de rampa limitado por `Wtraj_max`. Este parámetro puede actualizarse dinámicamente.
- Incorpora protección `wrapPi` para evitar saltos en la aritmética angular.

### B. Lazo de Velocidad (`PID_velocidad.py`)
- Control PI de velocidad clásico con integración hacia la referencia de corriente cuadrática (`Iq_ref`).
- **Lógica MTPA**: Calcula los límites dinámicos de corriente (`Id_mtpa_max`, etc.) utilizando la resolución de ecuaciones de segundo grado para maximizar el torque por amperio.
- **Anti-Windup Condicional en 2 Pasadas (`pass=1:2`)**: 
  - *Pasada 1*: Calcula la salida tentativa. Si hay saturación por los límites de voltaje/corriente del inversor y el integrador empeora dicha saturación, se "congela" el integrador (Back-calculation estricto).
  - *Pasada 2*: Recalcula la salida final con el estado del integrador estabilizado.

### C. Lazo de Corriente (`PID_corrientes.py`)
- Control vectorial PI sobre los ejes sincrónicos d-q.
- **Límite de Voltaje Dinámico**: Si la Inyección de Alta Frecuencia (HFI) está activa, reduce el margen de tensión fundamental disponible: `Vfund_max = max(0, Vmax - abs(Amplitud_HFI))`.
- **Anti-Windup por Back-Calculation**: Calculado independientemente para los ejes d y q: `Kaw_d = Ki_d / abs(Kp_d)`.

---

## 3. Observadores y Sensorless
Implementado en `observers.py`.

### A. Filtro Kalman Extendido (EKF) con HFI
- Discretizado a 4 estados mecánicos/eléctricos para estimar posición y velocidad del rotor basándose en la saliencia inductiva.
- Procesa armónicos inducidos por la señal HFI para estimar el ángulo en bajas velocidades.

### B. Observador de Perturbaciones (DOB) Unificado
- Integrado directamente sobre los estados matemáticos del EKF.
- Funciona como un filtro paso-bajos (ancho de banda = 25 Hz) sobre el error de modelo, estimando el torque total de perturbación (`d_hat`), que incluye fricciones, torque de carga externo y rizado de reluctancia (Cogging).

---

## 4. Feedforward de Cogging (Adaptativo en Tiempo Real)
Implementado en `feedforward/cogging_lut.py`. Una de las adiciones más robustas orientadas a la portabilidad hacia un Microcontrolador en C.

### Lógica `RealTimeCoggingLUT`
1. **Arreglo Circular**: Memoria preasignada (ej. 360 bins) que representa una vuelta eléctrica/mecánica completa.
2. **Modo de Aprendizaje (Identificación)**:
   - Durante las primeras `N` vueltas (ej. `turns = 2.0`), el feedforward inyectado es $0$.
   - El sistema captura el torque estimado por el DOB (`d_hat_dob`).
   - **Compensación de Desfase**: Como el DOB tiene un retraso provocado por su filtro ($T_d = \frac{1}{2\pi \cdot 25}$), el algoritmo calcula matemáticamente a qué ángulo físico pertenecía esa perturbación: $\theta_{true} = \theta_m - \omega_m \cdot T_d$.
   - Se actualiza el bin correspondiente mediante un filtro de paso bajo simple (EWMA): `LUT[bin] = (1 - \alpha)*LUT[bin] + \alpha * d_hat`.
3. **Modo de Aplicación (Feedforward)**:
   - Al completar el aprendizaje, el LUT detiene su actualización principal y comienza a inyectar el torque guardado directamente como prealimentación (Feedforward) al PID de velocidad.
   - **Resultado Crítico**: Una vez que el Feedforward de cogging opera, la velocidad se estabiliza. Con el rizado interno cancelado, la señal del DOB (`d_hat_dob`) queda completamente limpia y dedicada **únicamente** a estimar posibles perturbaciones externas desconocidas.

---

## 5. Ruteo Dinámico y Pruebas (`main.py`)
El archivo principal soporta un diccionario de conexión de señales extremadamente flexible, actuando como interruptores lógicos (simulando "Go-To" blocks de Simulink).
Permite configurar el sistema en:
- **Modo Sensado**: Realimentación de encoders ideales (`theta_m_real`, `Wm_real`).
- **Modo Híbrido**: Cualquier combinación de señales reales estimadas en diferentes lazos.
- **Modo Sensorless**: Operación 100% sobre estimaciones del EKF (`theta_m_ekf`, `Wm_ekf`).
