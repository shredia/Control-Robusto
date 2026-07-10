# Estructura del Proyecto: Control Robusto de Clinostato

Bienvenido a las notas de documentación del proyecto de **Control Robusto para un Clinostato 2D con Motores Stepper**.

Esta bóveda de Obsidian (o carpeta de notas) está diseñada para que cualquier persona que se integre al proyecto pueda comprender rápidamente **qué** estamos haciendo, **por qué** lo estamos haciendo y **cómo** está estructurado el código fuente.

## Índice de Contenidos

1. [[01_Descripcion_del_Proyecto]] - ¿Qué es un clinostato? ¿Cuál es el problema a resolver?
2. [[02_Arquitectura_de_Software]] - ¿Cómo está estructurada la simulación en Python? (Archivos, módulos y responsabilidades).
3. [[03_Teoria_de_Control_y_Observadores]] - Explicación de los algoritmos utilizados (FOC, EKF, HFI).
4. [[04_Como_Ejecutar_y_Experimentar]] - Guía para correr simulaciones y probar nuevos escenarios.
5. [[05_Formulas_y_Matematicas]] - Detalle matemático de la saturación de voltaje/corriente, separación frecuencial de lazos y cálculo de ganancias (Kp, Ki, Kd, DOB).
6. [[06_EKF_y_HFI_Detalle]] - Crítica técnica y ecuaciones detalladas del funcionamiento interno del EKF, la matriz Jacobiana y el demodulador HFI.

---
**Nota para nuevos desarrolladores:**
El objetivo principal de este proyecto es lograr un control de velocidad estable a bajas revoluciones sin utilizar sensores mecánicos (encoders), estimando todo a partir de las corrientes eléctricas del motor.
