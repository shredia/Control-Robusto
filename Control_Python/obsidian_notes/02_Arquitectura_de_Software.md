# Arquitectura de Software (Simulación en Python)

El entorno de simulación en Python (ubicado en `Control_Python/`) está diseñado de forma modular para facilitar el testeo de algoritmos de control antes de pasarlos a hardware.

## Módulos Principales

### `main.py`
Es el **punto de entrada principal** para una simulación estándar. 
- Define las funciones de referencia de posición, velocidad y el torque de carga simulado.
- Instancia y ejecuta la simulación.
- Exporta los resultados en un formato interactivo HTML (`scope_stepper.html`) que simula el comportamiento de un "Scope" de Simulink.

### `parameters.py`
Contiene los **parámetros físicos y de control** del sistema.
- Parámetros del motor Stepper (Resistencia, Inductancia $L_d$, $L_q$, Inercia, Flujo magnético, Polos).
- Ganancias de los controladores PI.
- Matrices de covarianza $Q$ y $R$ para la sintonización del Filtro de Kalman Extendido (EKF).

### `pmsm_model.py`
Define la **planta (el motor real en simulación)**. 
- Contiene las ecuaciones diferenciales (eléctricas y mecánicas) del Motor Síncrono de Imanes Permanentes (PMSM) / Stepper híbrido.
- Se encarga de calcular las derivadas de las corrientes y la velocidad mecánica.

### `simulation.py`
Gestiona el **ciclo principal de la simulación**.
- Implementa el lazo de integración numérica (por ejemplo, Euler o Runge-Kutta).
- Conecta el controlador (FOC), el estimador (EKF en `observers.py`), el inversor (PWM en `inverter.py`) y la planta (`pmsm_model.py`).
- Guarda el historial de las variables a través del tiempo.

### `observers.py`
Contiene la lógica de los **observadores de estado**.
- Aquí vive la implementación del **EKF (Extended Kalman Filter)**, el cual predice y corrige los estados del motor (corrientes, velocidad, posición y perturbación de carga).

### `run_experiments.py`
Un script auxiliar diseñado para **correr múltiples escenarios de prueba** de forma automatizada y comparar resultados (ej: comportamiento a lazo abierto vs. lazo cerrado).

### `inverter.py`
Simula el **Inversor Trifásico / Bifásico**. Convierte las señales de control de voltaje deseadas en voltajes aplicables a las fases del motor, considerando las limitaciones de modulación.
