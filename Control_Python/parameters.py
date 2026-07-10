import numpy as np

# ======================================================================
# CONFIGURACIÓN DE REALIMENTACIÓN (FEEDBACK SELECTION)
# ======================================================================
# Define de dónde se extraen las variables para cerrar cada lazo de control:
# Opciones disponibles: 
#   - 'real' : Usa los sensores físicos/medición directa de la planta
#   - 'ekf'  : Usa los estados estimados del Filtro de Kalman Extendido (EKF + HFI)
pos_feedback_source = 'real'      # Lazo de posición: 'real' o 'ekf'
speed_feedback_source = 'real'    # Lazo de velocidad: 'real' o 'ekf'
current_feedback_source = 'real'  # Lazo de corrientes (Park): 'real' o 'ekf'

# ======================================================================
# 1. PARÁMETROS FÍSICOS DEL MOTOR STEPPER BIFÁSICO (CON SALIENCIA)
# ======================================================================
Step_angle = 1.8            # pasos en grados del motor
N_phases = 2                # numero de phases
N_steps = 360 / Step_angle  # numero de pasos del motor (200)
N_teeths = N_steps / N_phases # Cantidad de dientes del rotor (100)
P = N_teeths / 2            # Número de pares de polos (50)

R_s = 2.5                   # Resistencia de estator (Ohm)
L_d = 3.66e-3               # Inductancia de eje d (Henry)
L_q = 9.61e-3               # Inductancia de eje q (Henry)

# Inercia y fricción interna nominales (utilizadas por el controlador/observador)
J_internal = 54e-7          # Inercia del motor (kg*m^2)
B_internal = 1e-4           # Fricción viscosa (N*m*s/rad)

# Carga externa (para simular desajuste de parámetros)
J_external = J_internal * 1.0
B_external = B_internal * 1.0

# Parámetros REALES de la planta física
J_real = J_internal + J_external  # 108e-7 kg*m^2
B_real = B_internal + B_external  # 2e-4 N*m*s/rad

# Variables nominales del modelo del controlador
J_var = J_internal
B_var = B_internal

# Parámetros para el torque de cogging/detent
Tdm = 18e-3                 # Torque de detent (N*m)
Phi = np.pi / 2.0           # Fase inicial del cogging (rad)

# ======================================================================
# 2. CONDICIONES INICIALES DEL MOTOR
# ======================================================================
i_d0 = 0.0                  # Corriente eje d inicial (A)
i_q0 = 0.0                  # Corriente eje q inicial (A)
omega_m0 = 0.0              # Velocidad mecánica angular inicial (rad/s)
theta_m0 = 0.0              # Posición mecánica angular inicial (rad)

# ======================================================================
# 3. PARÁMETROS DE SIMULACIÓN Y TIEMPOS DE MUESTREO (MULTITASA)
# ======================================================================
t_end = 10.0                # Tiempo total de simulación (s)
dt_sim = 1e-5               # Paso de simulación continua (RK4) (10 us)

frecuency_simulation = 10e3
f_carrier = 20e3
sample_time = 1.0 / frecuency_simulation

f_current = f_carrier / 2.0
f_ekf = f_current
Ts_ekf = 1.0 / f_ekf

Ts_current = 1.0 / f_current
Ts_Wm = Ts_current * 10.0
Ts_pos = Ts_Wm * 10.0

dt_current = Ts_current
dt_speed = Ts_Wm
dt_pos = Ts_pos
dt_observer = Ts_current

f_wm = 1.0 / Ts_Wm
Ts_DO = Ts_Wm

# Voltaje del bus DC y límites de salida
V_dc = 24.0                 # Voltaje de alimentación del bus DC (V)
V_max = V_dc                # En inversor bifásico, la tensión máxima por fase es V_dc
I_nom = 1.7                 # Corriente nominal del motor (A) (acorde a InitFcn.m)
i_max = I_nom * np.sqrt(2)  # Corriente máxima permitida (A)
Thold = 392e-3              # Holding torque (N*m)

Kt = Thold / (np.sqrt(2) * I_nom)
psi_m = Kt / P
Ke = Kt
omega_max = 30.0            # Velocidad mecánica máxima de consigna (rad/s)

# ======================================================================
# 4. SINTONIZACIÓN DE GANANCIAS EN CASCADA POR ASIGNACIÓN DE POLOS
# ======================================================================
fbw_d = 500.0               # Hz (BW corriente d)
fbw_q = 500.0               # Hz (BW corriente q)

wd_d = 2.0 * np.pi * fbw_d    # rad/s
wd_q = 2.0 * np.pi * fbw_q    # rad/s 

fbw_Wm = 50.0               # Hz (Bw Wm)

shi_d = 0.707
shi_q = 0.707

# Ganancias de corrientes d y q para motor con saliencia (L_d != L_q)
Kp_d_salient = 2.0 * shi_d * wd_d * L_d - R_s
Ki_id = (wd_d ** 2) * L_d
Kp_id = Kp_d_salient
Kd_id = 0.0

Kp_q_salient = 2.0 * shi_q * wd_q * L_q - R_s
Ki_iq = (wd_q ** 2) * L_q
Kp_iq = Kp_q_salient
Kd_iq = 0.0

# Sintonización de velocidad usando parámetros del modelo (J_var, B_var)
shi_w = 1.0
wn_w = 2.0 * np.pi * fbw_Wm

Kp_w = (2.0 * shi_w * wn_w * J_var - B_var) / Ke
Ki_w = (wn_w ** 2) * J_var / Ke
Kd_w = 0.0
Tf_w = 10.0 / 1000.0

Kp_speed = Kp_w
Ki_speed = Ki_w
Kd_speed = Kd_w

# Ganancias de posición (control PD)
Kp_pos = 1.0                # Ganancia Proporcional (rigidez)
Ki_pos = 0.0                # Ganancia Integral
Kd_pos = 0.05               # Ganancia Derivativa (amortiguamiento ante cogging)

# ======================================================================
# 5. PARÁMETROS DEL OBSERVADOR (SMO, LOB & EKF+HFI+DOB)
# ======================================================================
smo_l_gain = 80.0           # Ganancia de conmutación del SMO
smo_f_cutoff = 300.0        # Frecuencia de corte del filtro pasabajos (Hz)
smo_omega_pll = 250.0       # Frecuencia natural del PLL (rad/s)

lob_p_obs = 50.0            # Polo del observador de torque de carga (rad/s)

# Sintonización del Observador de Perturbaciones (DOB de 3er orden para EKF)
fbw_DO = 50.0                  # Hz (Ancho de banda del DOB ampliado para captar cogging)
wn_do = 2.0 * np.pi * fbw_DO   # rad/s
BJ_hat = B_var / J_var

B0 = 3.0 * wn_do - BJ_hat
B1 = 3.0 * (wn_do ** 2) - B0 * BJ_hat
B2 = (wn_do ** 3) * (J_var / P)

beta0 = B0
beta1 = B1
beta2 = B2

# Parámetros del inyector HFI y demodulación
f_h = 800.0                  # Frecuencia de inyección en Hz
wh = 2.0 * np.pi * f_h        # Frecuencia en rad/s
Amplitud_HFI = 2.0            # Voltios de inyección (aumentado para mayor SNR)
HFI_enable = 1                # 1: Habilitado, 0: Deshabilitado

# Coeficientes del filtro BPF HFI (Butterworth @ Fs=10kHz)
b0_hfi = 0.030469
b1_hfi = 0.0
b2_hfi = -0.030469
a1_hfi = -1.822695
a2_hfi = 0.939062
