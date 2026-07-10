# Descripción del Proyecto

## ¿Qué estamos haciendo?
Estamos desarrollando e implementando un sistema de **control en lazo cerrado sin sensores (sensorless)** para un **clinostato 2D** accionado por **motores paso a paso (steppers híbridos)**. 

El clinostato es una "máquina de ingravidez" que hace girar sus ejes para promediar temporalmente el vector de gravedad, simulando un entorno de microgravedad útil para experimentos espaciales y biológicos.

## ¿Por qué lo estamos haciendo? (La Problemática)
1. **Transmisión por poleas y perturbaciones:** Para reducir costos y evitar cableado complejo (slip-rings), el sistema usa poleas. Esto introduce un torque resistente variable dependiendo de la posición angular, lo cual actúa como una perturbación mecánica.
2. **Ausencia de sensores (Sensorless):** El diseño mecánico dificulta usar encoders para medir velocidad o posición continuamente.
3. **Pérdida de sincronía a lazo abierto:** Si operamos un motor stepper a lazo abierto a bajas velocidades con perturbaciones de carga variables, se producen vibraciones, desviaciones de velocidad e incluso pérdida de sincronía.

## Nuestra Solución Propuesta
Para lograr un movimiento rotacional constante (velocidad estable) a pesar de las perturbaciones y sin sensores mecánicos, utilizamos:
- **FOC (Field Oriented Control):** Control orientado a campo para regular el torque independizando las corrientes magnéticas y de par.
- **EKF (Extended Kalman Filter):** Un observador de estados que, a partir de señales eléctricas (voltaje y corriente), estima la **velocidad mecánica**, la **posición del rotor** y compensa las **perturbaciones** (como el torque de carga variable o fricción).
- **HFI (High-Frequency Injection):** El EKF pierde observabilidad a bajas velocidades porque la fuerza contraelectromotriz (BEMF) es muy pequeña. Inyectamos alta frecuencia para aprovechar la saliencia magnética del motor y estimar la posición incluso a velocidad nula.
