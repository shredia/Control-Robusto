# Teoría de Control y Observadores

Este proyecto combina control clásico con estimación avanzada de estados. A continuación se detallan los tres pilares de nuestra estrategia:

## 1. FOC (Field Oriented Control)
El control orientado a campo nos permite controlar el motor stepper como si fuera un motor DC tradicional.
- Transformamos las corrientes estacionarias de las fases (por ejemplo, bifásicas $\alpha-\beta$ o trifásicas $abc$) a un marco de referencia rotatorio ($d-q$).
- **Eje $d$ (Flujo):** La corriente $i_d$ controla el flujo magnético (comúnmente se regula a 0 para maximizar la eficiencia).
- **Eje $q$ (Torque):** La corriente $i_q$ es proporcional al torque generado. Controlando $i_q$ controlamos la aceleración del sistema.

## 2. EKF (Filtro de Kalman Extendido)
Como no tenemos encoders, estimamos la velocidad y posición a partir de la medición de voltajes y corrientes eléctricas.
El motor stepper es un sistema **no lineal**, por lo que usamos la versión extendida del filtro de Kalman.
El EKF funciona en dos pasos:
1. **Predicción:** Usa el modelo matemático del motor para predecir el próximo estado (corriente, posición, velocidad, torque de carga).
2. **Corrección:** Compara la corriente eléctrica real medida con la corriente predicha y usa el error para ajustar (corregir) todos los estados estimados.

> **Nota sobre perturbaciones:** El EKF en este proyecto está ampliado para estimar un "estado extendido", que corresponde a la perturbación de carga (el torque variable del sistema de poleas). Al estimarlo, el control puede aplicar una pre-alimentación (Feed-forward) para cancelarlo.

## 3. HFI (Inyección de Alta Frecuencia)
**Problema:** A bajas velocidades (las cuales son típicas en un clinostato), el voltaje inducido por el giro del rotor (BEMF) es casi cero. El EKF depende de este BEMF para "observar" la posición. A velocidad nula, el sistema se vuelve **inobservable**.
**Solución (HFI):** Los motores stepper tienen saliencia magnética (las inductancias $L_d$ y $L_q$ son distintas). Al inyectar una señal de voltaje de alta frecuencia, las corrientes de respuesta contendrán información sobre la posición del rotor, sin importar la velocidad. Esta información alimenta al estimador para mantener la sincronía incluso al arrancar desde cero.
