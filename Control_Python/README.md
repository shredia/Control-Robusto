# Proyecto de control en Python

La entrada recomendada es `run.py`. Las decisiones de arquitectura están en
`project_config.py`; los parámetros físicos y las ganancias numéricas permanecen
en `parameters.py`.

Instalación del entorno de visualización:

```powershell
python -m pip install -r requirements.txt
```

## Ejecuciones habituales

Desde la carpeta `Control_Python`:

```powershell
# Planta con sensores ideales (referencia de comparación)
python run.py --profile sensed

# EKF electromecánico normal, sin HFI, DOB ni LMS
python run.py --profile ekf

# EKF normal con estimación/compensación de perturbación
python run.py --profile ekf-dob

# EKF electromecánico corregido mediante HFI
python run.py --profile ekf-hfi

# Método derivado de los papers: HFI + Kalman [theta, omega, alpha]
python run.py --profile papers

# Método de los papers con DOB y LMS
python run.py --profile papers-dob-lms
```

Para una prueba rápida:

```powershell
python run.py --profile papers --duration 0.05
```

Para generar gráficos:

```powershell
python run.py --profile ekf-hfi --plots
```

## Interruptores

Cada perfil se puede ajustar sin editar código:

```powershell
python run.py --profile papers --dob --lms
python run.py --profile ekf-hfi --no-dob --no-lms
python run.py --profile ekf --lob
```

El programa rechaza configuraciones incoherentes. Por ejemplo, `--lms --no-dob`
falla porque el LMS aprende a partir de la perturbación entregada por el DOB.

## Estructura

- `run.py`: entrada única para simulaciones normales.
- `project_config.py`: perfiles e interruptores de arquitectura.
- `parameters.py`: motor, tiempos de muestreo y ganancias.
- `simulation.py`: coordinación de planta, control y observadores.
- `observers.py`: EKF, HFI/Kalman cinemático, SMO y LOB.
- `control_blocks/`: controladores de posición, velocidad y corriente.
- `feedforward/`: estimador LMS de cogging.
- `pmsm_model.py` e `inverter.py`: planta física.
- `main.py`: visualización histórica; se conserva por compatibilidad.
- `run_experiments.py`: comparación histórica de cuatro casos; se conserva por
  compatibilidad, pero para nuevos experimentos debe usarse `run.py`.

## Significado de los observadores

- `real`: usa las variables reales de la planta.
- `ekf`: EKF electromecánico sin corrección angular HFI.
- `ekf_hfi`: EKF electromecánico con corrección angular HFI.
- `papers`: usa el ángulo HFI como medición de un Kalman cinemático con estados
  `[theta_e, omega_e, alpha_e]`; la velocidad mecánica es `omega_e / P`.

Los resultados principales siempre quedan disponibles como
`theta_m_observer`, `theta_e_observer` y `omega_m_observer`, independientemente
del observador seleccionado.
