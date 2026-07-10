# Cómo Ejecutar y Experimentar

## Pre-requisitos
Para ejecutar las simulaciones en Python necesitas tener instaladas las librerías científicas estándar:
- `numpy`
- `matplotlib`
- `plotly` (para los gráficos interactivos HTML)

Puedes instalarlas usando:
```bash
pip install numpy matplotlib plotly
```

## Ejecución Básica
Para correr una simulación simple y generar los gráficos, ejecuta:
```bash
python main.py
```
**Resultado:** Esto creará (o sobrescribirá) el archivo `outputs/scope_stepper.html`. Ábrelo en cualquier navegador web para explorar interactivamente los resultados de las variables como Posición, Velocidad, Corrientes, Torque, etc.

## Cambiar Parámetros del Motor o Control
Abre el archivo `parameters.py`. Aquí puedes ajustar:
- **`par.J`**: Inercia del motor.
- **`par.B`**: Coeficiente de fricción viscosa.
- **`par.Kp`, `par.Ki`**: Ganancias de los controladores de corriente o velocidad.
- **`par.Q`, `par.R`**: Matrices del EKF. Modifica `Q` para decirle al filtro cuánto confías en el modelo de cada estado, y `R` para cuánto confías en tus mediciones (ruido).

## Modificar el Escenario de Simulación
En `main.py` encontrarás funciones como:
- `get_position_reference(t)`
- `get_speed_reference(t)`
- `get_load_torque(t)`

Puedes modificar estas funciones (por ejemplo, para agregar escalones de carga o perfiles de velocidad variables) y ver cómo reacciona el control EKF ante estas perturbaciones.

## Experimentos Masivos
Si deseas probar varios controladores o escenarios y comparar gráficas de forma automatizada, usa:
```bash
python run_experiments.py
```
Este script está diseñado para ciclar entre diferentes parámetros y guardar los resultados de forma estructurada.
