# Estructura de tesis paso a paso

Esta nota propone una estructura para ordenar la tesis de forma progresiva: primero se presenta el problema físico y experimental, luego el sistema electromecánico, después el modelo matemático, y finalmente las estrategias de control y estimación.

## Idea central de la tesis

La tesis busca estudiar e implementar una estrategia de control sin sensores mecánicos para un clinostato 2D accionado por motores paso a paso. El foco está en mantener una velocidad baja, estable y uniforme pese a perturbaciones mecánicas variables producidas por la transmisión y la carga del sistema.

El hilo conductor recomendado es:

1. Por qué importa controlar bien un clinostato.
2. Qué problema aparece al usar motores stepper y transmisión por poleas.
3. Por qué el control en lazo abierto no basta.
4. Qué variables no se pueden medir directamente.
5. Cómo se modela el motor y el sistema.
6. Cómo se controla el motor mediante FOC y PI.
7. Cómo se estima velocidad/posición usando EKF.
8. Por qué se necesita HFI a baja velocidad.
9. Cómo se valida la propuesta por simulación y/o prototipo.

---

## Capítulo 1 - Introducción general

Este capítulo debe convencer al lector de que el problema existe y vale la pena resolverlo.

### 1.1 Contexto: clinostatos y simulación de microgravedad

Explicar qué es un clinostato y por qué se usa para modificar el promedio temporal del vector de gravedad.

Puntos clave:

- Aplicaciones en investigación espacial, biología, agricultura o estudios de crecimiento.
- Diferencia entre clinostatos de 1, 2 y 3 ejes.
- Necesidad de movimiento lento, continuo y estable.
- Importancia de mantener una velocidad aproximadamente constante.

### 1.2 Estado del arte

Esta sección debe mostrar qué se ha hecho y dónde aparece la brecha.

#### 1.2.1 Clinostatos y configuraciones mecánicas

Exponer las configuraciones típicas:

- Clinostatos con motores montados en ejes móviles.
- Uso de slip-rings para alimentación y sensado.
- Configuraciones con poleas para evitar cableado rotativo.
- Ventajas y desventajas de cada alternativa.

#### 1.2.2 Motores usados en clinostatos

Comparar brevemente:

- Motores DC.
- Motores stepper.
- Motores brushless/PMSM si corresponde como referencia conceptual.

Enfatizar que los motores stepper son atractivos por simplicidad, torque a baja velocidad y control incremental, pero pueden presentar vibraciones, pérdida de pasos y sensibilidad a perturbaciones.

#### 1.2.3 Control de motores stepper

Presentar el paso desde lo simple a lo avanzado:

- Control a lazo abierto.
- Microstepping.
- Control PI/PID.
- Control vectorial/FOC.
- Control sensorless.
- Observadores de estado.
- EKF.
- HFI para baja velocidad.

### 1.3 Problemática

Aquí debe quedar claro el conflicto principal:

El clinostato usa motores fijos y transmisión por poleas, lo que reduce complejidad mecánica, pero genera un torque resistente variable y dificulta instalar sensores mecánicos permanentes.

Subproblemas recomendados:

- Perturbación de torque dependiente de la posición.
- Ausencia de medición directa de velocidad y posición.
- Baja velocidad de operación, donde la fuerza contraelectromotriz entrega poca información.
- Posible pérdida de observabilidad.
- Vibración y rizado de corriente en operación a lazo abierto.

### 1.4 Propuesta de solución

Presentar la solución como una arquitectura integrada:

- Control FOC para desacoplar corrientes en ejes d-q.
- Control PI para corriente y velocidad.
- EKF para estimar estados no medibles.
- HFI para mejorar observabilidad a baja o nula velocidad.
- Compensación/estimación de perturbaciones mecánicas.

### 1.5 Objetivos

#### 1.5.1 Objetivo general

Implementar y evaluar un sistema de control sensorless basado en FOC, EKF y HFI para mejorar la estabilidad de velocidad de un clinostato 2D accionado por motores stepper bajo perturbaciones mecánicas no medidas.

#### 1.5.2 Objetivos específicos

Orden recomendado:

1. Modelar el motor stepper bifásico considerando dinámica eléctrica, mecánica y saliencia.
2. Caracterizar parámetros eléctricos y mecánicos del motor.
3. Modelar la perturbación de torque asociada al clinostato.
4. Diseñar el control FOC con lazos PI de corriente y velocidad.
5. Analizar observabilidad y controlabilidad del sistema.
6. Implementar un EKF para estimar velocidad y posición.
7. Incorporar HFI para mejorar observabilidad a baja velocidad.
8. Comparar desempeño entre lazo abierto, lazo cerrado sensado y lazo cerrado sensorless.

### 1.6 Alcances y limitaciones

Separar claramente qué sí se abordará y qué queda fuera.

Alcances:

- Simulación del motor y del control.
- Estimación de velocidad/posición mediante señales eléctricas.
- Validación en un prototipo o banco experimental acotado.
- Evaluación de seguimiento de velocidad y rechazo de perturbaciones.

Limitaciones:

- Sensibilidad al ruido eléctrico.
- Dependencia del modelo usado por el EKF.
- Recursos computacionales requeridos.
- Validación experimental limitada respecto al sistema final completo.

---

## Capítulo 2 - Marco teórico

Este capítulo debe entregar las herramientas necesarias para entender la metodología, sin adelantar resultados.

### 2.1 Motor stepper bifásico

Explicar:

- Principio de funcionamiento.
- Fases A/B desfasadas 90 grados.
- Paso completo, microstepping y movimiento continuo.
- Ventajas a baja velocidad.
- Problemas de lazo abierto.

### 2.2 Modelo electromecánico del motor

Presentar el modelo como base para control y estimación.

Subsecciones sugeridas:

#### 2.2.1 Dinámica eléctrica en fases A-B

Incluir corrientes, voltajes, resistencia, inductancias y fuerza contraelectromotriz.

#### 2.2.2 Dinámica mecánica

Incluir torque electromagnético, inercia, fricción, torque de carga y torque perturbador.

#### 2.2.3 Relación entre ángulo mecánico y eléctrico

Explicar el rol del número de pares de polos o dientes del rotor.

### 2.3 Saliencia inductiva

Explicar por qué la inductancia cambia con la posición del rotor.

Puntos clave:

- Diferencia entre Ld y Lq.
- Inductancia promedio L0.
- Componente variable L2.
- Importancia de la saliencia para HFI.

### 2.4 Transformada de Park

Mostrar por qué conviene pasar del marco estacionario A-B al marco rotatorio d-q.

Orden recomendado:

1. Sistema bifásico original.
2. Rotación del marco de referencia.
3. Separación entre eje de flujo d y eje de torque q.
4. Ventaja para control PI y FOC.

### 2.5 Control orientado a campo en motores stepper

Explicar FOC como una forma de controlar torque y flujo de manera desacoplada.

Subtítulos sugeridos:

#### 2.5.1 Corriente d como flujo

#### 2.5.2 Corriente q como torque

#### 2.5.3 Desacoples eléctricos

#### 2.5.4 Lazos de corriente y velocidad

### 2.6 Observadores de estado

Introducir la idea antes del EKF.

Explicar:

- Qué es un estado.
- Qué significa estimar una variable no medida.
- Por qué no basta medir solo corrientes.
- Qué rol tiene el modelo matemático.

### 2.7 Filtro de Kalman Extendido

Presentar el EKF como observador para sistemas no lineales.

Subsecciones sugeridas:

#### 2.7.1 Predicción del modelo

#### 2.7.2 Corrección con mediciones

#### 2.7.3 Matrices de ruido Q y R

#### 2.7.4 Estados estimados: corriente, velocidad, posición y perturbación

### 2.8 Observabilidad y controlabilidad

Explicar por qué el sistema puede ser controlable pero no completamente observable en ciertas condiciones.

Puntos clave:

- Rango de matriz de observabilidad.
- Rango de matriz de controlabilidad.
- Pérdida de información a baja velocidad.
- Dependencia del punto de operación.

### 2.9 High Frequency Injection

Presentar HFI como solución al problema de baja observabilidad.

Orden recomendado:

1. A baja velocidad hay poca fuerza contraelectromotriz.
2. Se inyecta una señal de alta frecuencia en el eje d.
3. La saliencia produce una respuesta dependiente de la posición.
4. Esa respuesta permite corregir la estimación angular.

---

## Capítulo 3 - Metodología

Este capítulo debe explicar qué se hizo, en qué orden y con qué criterios.

### 3.1 Caracterización del motor

Describir el motor usado: stepper bifásico 17HS3401S.

Subsecciones sugeridas:

#### 3.1.1 Parámetros de datasheet

Incluir tabla con:

- Voltaje nominal.
- Corriente por fase.
- Resistencia.
- Inductancia.
- Torque de mantenimiento.
- Torque de detención.
- Inercia.
- Ángulo de paso.
- Número de fases.

#### 3.1.2 Parámetros derivados

Explicar cálculo de:

- Pasos por revolución.
- Constante de torque/fuerza.
- L0.
- Ld y Lq si aplica.
- Número de pares de polos o relación eléctrica-mecánica.

#### 3.1.3 Medición experimental de inductancia

Exponer el procedimiento de medición con tester LCR.

Puntos clave:

- Medición de La y Lb para distintas posiciones.
- Identificación de valores máximos y mínimos.
- Cálculo de Ld, Lq, L0 y L2.
- Validación contra datasheet.

### 3.2 Modelado del sistema

Separar el modelo en bloques.

#### 3.2.1 Modelo eléctrico A-B

#### 3.2.2 Modelo mecánico

#### 3.2.3 Modelo en coordenadas d-q

#### 3.2.4 Modelo de perturbación de torque

Aquí conviene explicar cómo se representa el torque resistente variable del clinostato.

### 3.3 Diseño del control PI/FOC

Mostrar la arquitectura de control.

Subsecciones sugeridas:

#### 3.3.1 Lazo de corriente d

Objetivo: regular Id, usualmente cerca de cero o del valor definido por estrategia.

#### 3.3.2 Lazo de corriente q

Objetivo: regular Iq porque está asociada al torque.

#### 3.3.3 Lazo de velocidad

Objetivo: convertir error de velocidad en referencia de torque/corriente q.

#### 3.3.4 Desacoples del modelo

Explicar términos compensados para tratar los ejes como sistemas SISO.

#### 3.3.5 Saturaciones y anti-windup

Si se implementó en simulación/Python, conviene incluirlo porque fortalece la tesis desde el punto de vista práctico.

### 3.4 Diseño del EKF

Explicar el estimador paso a paso.

Subsecciones sugeridas:

#### 3.4.1 Estados estimados

Ejemplo:

- Id.
- Iq.
- Velocidad mecánica.
- Posición mecánica o eléctrica.
- Perturbación, si se incluye como estado aumentado.

#### 3.4.2 Entradas del EKF

- Vd.
- Vq.
- Parámetros del motor.

#### 3.4.3 Mediciones disponibles

- Corriente de fase o corrientes d-q.
- Voltajes aplicados o referencias de voltaje.

#### 3.4.4 Predicción

#### 3.4.5 Corrección

#### 3.4.6 Sintonización de Q y R

### 3.5 Análisis de observabilidad y controlabilidad

Mostrar el análisis linealizado.

Orden sugerido:

1. Definir vector de estados.
2. Definir vector de entradas.
3. Definir salidas medidas.
4. Calcular matriz A.
5. Calcular matriz B.
6. Construir matriz de observabilidad.
7. Construir matriz de controlabilidad.
8. Interpretar rangos.

Resultado esperado según la tesis actual:

- Observabilidad parcial en el punto de operación analizado.
- Controlabilidad completa en el punto de operación analizado.
- Necesidad de HFI para mejorar estimación a baja velocidad.

### 3.6 Implementación de HFI

Explicar cómo se agrega la señal de alta frecuencia.

Subsecciones sugeridas:

#### 3.6.1 Señal inyectada

#### 3.6.2 Respuesta esperada en corriente

#### 3.6.3 Detección de error angular

#### 3.6.4 Integración con EKF o lazo de estimación

### 3.7 Compensación de perturbaciones

Esta sección puede conectar con el torque variable del clinostato.

Opciones de subtítulos:

#### 3.7.1 Observador de perturbaciones

#### 3.7.2 Compensación feedforward

#### 3.7.3 Tabla/LUT de cogging o torque periódico

#### 3.7.4 Evaluación del rechazo de perturbación

### 3.8 Escenarios de simulación y validación

Orden recomendado:

1. Motor a lazo abierto sin perturbación.
2. Motor a lazo abierto con perturbación.
3. Control FOC con medición ideal.
4. Control FOC con EKF.
5. Control FOC + EKF + HFI a baja velocidad.
6. Control con compensación de perturbación.

Métricas sugeridas:

- Error de velocidad.
- Rizado de velocidad.
- Error de posición/ángulo estimado.
- Rizado de corriente.
- Esfuerzo de control.
- Tiempo de establecimiento.
- Robustez frente a cambios de carga.

---

## Capítulo 4 - Resultados esperados o resultados

Si aún no están todos los resultados, este capítulo puede planificarse desde ahora.

### 4.1 Validación del modelo del motor

Comparar comportamiento simulado con lo esperado según datasheet o mediciones.

### 4.2 Desempeño en lazo abierto

Mostrar por qué el lazo abierto no es suficiente.

### 4.3 Desempeño del control FOC sensado

Usarlo como caso base ideal.

### 4.4 Desempeño del EKF

Evaluar estimación de velocidad y posición.

### 4.5 Desempeño con HFI a baja velocidad

Mostrar mejora de observabilidad.

### 4.6 Rechazo de perturbaciones

Comparar antes/después de compensación.

### 4.7 Comparación global

Tabla recomendada:

| Caso | Error de velocidad | Rizado | Error de estimación | Esfuerzo de control | Comentario |
| --- | --- | --- | --- | --- | --- |
| Lazo abierto | Alto | Alto | No aplica | Bajo | Sensible a perturbaciones |
| FOC sensado | Bajo | Bajo | No aplica | Medio | Caso ideal |
| FOC + EKF | Medio/Bajo | Medio | Depende de velocidad | Medio | Pierde calidad a baja velocidad |
| FOC + EKF + HFI | Bajo | Bajo | Mejorado | Medio/Alto | Mejor opción sensorless |

---

## Capítulo 5 - Conclusiones

Las conclusiones deberían responder directamente a los objetivos.

### 5.1 Cumplimiento de objetivos

Relacionar cada objetivo específico con lo logrado.

### 5.2 Aportes principales

Posibles aportes:

- Modelo del motor stepper con saliencia para control sensorless.
- Integración de FOC, EKF y HFI en un clinostato 2D.
- Análisis de observabilidad a baja velocidad.
- Estrategia para enfrentar perturbaciones mecánicas sin sensores adicionales.

### 5.3 Limitaciones observadas

Retomar ruido, dependencia del modelo y validación experimental limitada.

### 5.4 Trabajo futuro

Ideas:

- Implementación embebida en microcontrolador.
- Validación experimental de larga duración.
- Identificación más precisa de fricción y torque de carga.
- Optimización automática de parámetros del EKF.
- Comparación con otros observadores.

---

## Anexos recomendados

### Anexo A - Datasheet del motor

Guardar especificaciones completas del 17HS3401S.

### Anexo B - Mediciones de inductancia

Incluir imágenes y tabla de mediciones.

### Anexo C - Parámetros de simulación

Incluir valores usados en MATLAB/Simulink y Python.

### Anexo D - Código o pseudocódigo del EKF

Evitar pegar código excesivo en el cuerpo principal de la tesis.

### Anexo E - Diagramas de control

Incluir diagramas completos de FOC, PI, EKF, HFI y compensación.

---

## Orden narrativo recomendado

Para que la exposición sea clara, conviene escribir cada parte con esta lógica:

1. Problema físico: el clinostato necesita velocidad estable.
2. Restricción mecánica: se usan poleas y motores fijos para evitar slip-rings.
3. Consecuencia: aparece torque variable y no hay medición directa del eje.
4. Restricción de control: el lazo abierto no compensa perturbaciones.
5. Modelo: se representa el motor stepper en A-B y luego en d-q.
6. Control: FOC permite regular torque mediante Iq.
7. Estimación: EKF reemplaza sensores mecánicos usando señales eléctricas.
8. Dificultad: a baja velocidad se pierde observabilidad.
9. Solución complementaria: HFI usa saliencia para recuperar información angular.
10. Validación: se comparan escenarios y métricas.

---

## Notas relacionadas

- [[Notas_Migracion_Obsidian]]
- [[03_Teoria_de_Control_y_Observadores]]
- [[05_Formulas_y_Matematicas]]
- [[06_EKF_y_HFI_Detalle]]
- [[FOC_Control]]
- [[Observers]]
- [[PMSM_Model]]
