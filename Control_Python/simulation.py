import numpy as np
import parameters as par
from project_config import ProjectConfig
from pmsm_model import StepperMotor2Phase
from inverter import TwoPhaseInverter
def wrapPi(x):
    """
    Acota ángulos periódicos entre -pi y pi.
    """
    return np.arctan2(np.sin(x), np.cos(x))

def park(i_alpha, i_beta, theta_e):
    """
    Transformada de Park.
    """
    cos_t = np.cos(theta_e)
    sin_t = np.sin(theta_e)
    i_d =  i_alpha * cos_t + i_beta * sin_t
    i_q = -i_alpha * sin_t + i_beta * cos_t
    return i_d, i_q
from control_blocks.PID_posicion import PIDPositionHold as PID_position_hold
from control_blocks.PID_velocidad import PISpeedRIMTPA as PI_speed_RI_MTPA
from control_blocks.PID_corrientes import PIDCorrientes as PI_dq_vectorial_fixed
from observers import EKF_HFI_Observer, HfiKinematicKalman, LoadTorqueObserver
from feedforward.lms_cogging import RealTimeLMSCogging

class PMSMSimulation:
    """
    Coordinador de simulación multitasa en tiempo mixto (continuo-discreto)
    para un Motor Stepper Bifásico con control en cascada de 3 niveles y observadores.
    Instancia y mapea directamente los bloques de control equivalentes a tus S-Functions.
    """
    def __init__(self, t_end=None, dt_sim=None, dt_pos=None, dt_speed=None, dt_current=None,
                 dt_observer=None, inverter_mode='average',
                 signals_routing=None, config=None):
        
        # Configuración de enrutamiento de señales
        self.config = (config or ProjectConfig(inverter_mode=inverter_mode)).validate()
        self.signals_routing = signals_routing or self._default_signal_routing()
        
        # Cargar valores por defecto desde parameters.py
        self.t_end = t_end if t_end is not None else par.t_end
        self.dt_sim = dt_sim if dt_sim is not None else par.dt_sim
        self.dt_pos = dt_pos if dt_pos is not None else par.dt_pos
        self.dt_speed = dt_speed if dt_speed is not None else par.dt_speed
        self.dt_current = dt_current if dt_current is not None else par.dt_current
        self.dt_observer = dt_observer if dt_observer is not None else par.dt_observer
        self.inverter_mode = self.config.inverter_mode
        
        # Inicializar el motor paso a paso con las condiciones físicas REALES (plant)
        self.motor = StepperMotor2Phase(
            R_s=par.R_s, L_d=par.L_d, L_q=par.L_q, psi_m=par.psi_m,
            J=par.J_real, B=par.B_real, P=par.P,
            Tdm=par.Tdm, N_steps=par.N_steps, Phi=par.Phi
        )
        self.motor.set_states([par.i_d0, par.i_q0, par.omega_m0, par.theta_m0])
        
        # 5. Instanciación del estimador LMS de Cogging (Aprende durante 2 vueltas)
        self.cogging_lut = RealTimeLMSCogging(N_steps=par.N_steps, fbw_DO=25.0, learning_rate=0.01, turns_to_learn=2.0)
        
        # Inicializar inversor bifásico
        self.inverter = TwoPhaseInverter(V_dc=par.V_dc, f_sw=1.0/self.dt_current)
        
        # Inicializar los tres bloques de control independientes (estilo Simulink)
        self.pos_controller = PID_position_hold()
        self.speed_controller = PI_speed_RI_MTPA()
        self.current_controller = PI_dq_vectorial_fixed()
        
        # Inicializar observador de torque de carga clásico (LOB)
        self.lob = LoadTorqueObserver(
            J=par.J_var, B=par.B_var, dt_obs=self.dt_observer, p_obs=par.lob_p_obs
        )

        # Inicializar el observador híbrido EKF+HFI+DOB de tu modelo de Simulink
        self.ekf_hfi = EKF_HFI_Observer(
            Ke=par.Ke, Kt=par.Kt, R=par.R_s, J=par.J_var, Ld=par.L_d, Lq=par.L_q,
            Nr=par.P, N_steps=par.N_steps, B=par.B_var, Tdtm=par.Tdm,
            beta0=par.beta0, beta1=par.beta1, beta2=par.beta2,
            wh=par.wh, b0=par.b0_hfi, b1=par.b1_hfi, b2=par.b2_hfi,
            a1=par.a1_hfi, a2=par.a2_hfi, Ts=self.dt_observer
        )

        self.hfi_speed_kf = HfiKinematicKalman(
            pole_pairs=par.P,
            dt=self.dt_observer,
            q_jerk=par.hfi_kf_q_jerk,
            sigma_theta=par.hfi_kf_sigma_theta,
        )

        # Las opciones se fijan aquí; los scripts no deben mutar internamente los bloques.
        self.ekf_hfi.use_bpf = self.config.enable_bpf
        self.ekf_hfi.force_hfi_off = not self.config.enable_hfi
        q_tx = 1e-8 if self.config.enable_dob else 1e-12
        self.ekf_hfi.Qk_override = np.diag([1e-4, 1e-4, 1e-7, 1e-10, q_tx])

        # Referencias intermedias retenidas para ZOH
        self.omega_ref_held = 0.0
        self.i_q_ref_held = 0.0
        self.i_d_ref_held = 0.0

    def _default_signal_routing(self):
        source = self.config.feedback_source
        if source == "papers":
            theta_m, omega_m, theta_e = "theta_m_hfi_kf", "Wm_hfi_kf", "theta_e_hfi_kf"
        elif source == "ekf":
            theta_m, omega_m, theta_e = "theta_m_ekf", "Wm_ekf", "theta_e_ekf"
        else:
            theta_m, omega_m, theta_e = "theta_m_real", "Wm_real", "theta_e_real"
        return {
            "pos_theta_m": theta_m,
            "pos_Wm": omega_m,
            "speed_Wm_ref": "omega_ref_ext",
            "speed_Wm": omega_m,
            "curr_Wm": omega_m,
            "curr_Theta_e": theta_e,
        }

    def run(self, get_theta_ref, get_load_torque, get_omega_ref=None):
        """
        Ejecuta la simulación completa del sistema.
        """
        t_arr = np.arange(0.0, self.t_end, self.dt_sim)
        N = len(t_arr)
        
        # Estructura para registrar datos (historial)
        history = {
            't': t_arr,
            'theta_ref': np.zeros(N), 'theta_m': np.zeros(N), 'theta_m_est': np.zeros(N),
            'omega_ref': np.zeros(N), 'omega_m': np.zeros(N), 'omega_m_est': np.zeros(N),
            'theta_m_observer': np.zeros(N), 'omega_m_observer': np.zeros(N),
            'theta_e_observer': np.zeros(N),
            'theta_e': np.zeros(N), 'theta_e_est': np.zeros(N),
            'theta_hfi_e': np.zeros(N), 'amp_hfi': np.zeros(N), 'threshold_hfi': np.zeros(N), 'hfi_valid': np.zeros(N),
            'theta_e_hfi_kf': np.zeros(N), 'omega_e_hfi_kf': np.zeros(N),
            'alpha_e_hfi_kf': np.zeros(N), 'omega_m_hfi_kf': np.zeros(N),
            'innovation_hfi_kf': np.zeros(N),
            'e_w': np.zeros(N), 'e_theta_hfi': np.zeros(N), 'e_theta_ekf': np.zeros(N),
            'K_w': np.zeros(N), 'K_tx': np.zeros(N),
            'i_a': np.zeros(N), 'i_b': np.zeros(N),
            'i_a_hfi': np.zeros(N), 'i_b_hfi': np.zeros(N),
            'i_a_fund': np.zeros(N), 'i_b_fund': np.zeros(N),
            'i_d': np.zeros(N), 'i_q': np.zeros(N),
            'i_d_ref': np.zeros(N), 'i_q_ref': np.zeros(N),
            'i_d_est': np.zeros(N), 'i_q_est': np.zeros(N),
            'i_a_est': np.zeros(N), 'i_b_est': np.zeros(N),
            'torque_e': np.zeros(N), 'torque_l': np.zeros(N),
            'torque_l_est': np.zeros(N), 'd_hat_dob': np.zeros(N),
            'T_cogging_est': np.zeros(N), 'T_cogging_real': np.zeros(N),
            'v_d': np.zeros(N), 'v_q': np.zeros(N)
        }
        
        # Con consignas iniciales retenidas (ZOH)
        v_a_cmd, v_b_cmd, v_d_cmd, v_q_cmd = 0.0, 0.0, 0.0, 0.0
        
        T_L_est = 0.0
        d_hat_dob = 0.0
        x_est = np.zeros(5)
        amp_hfi_val = 0.0
        theta_hfi_e_val = 0.0
        theta_e_hfi_kf = 0.0
        omega_e_hfi_kf = 0.0
        alpha_e_hfi_kf = 0.0
        omega_m_hfi_kf = 0.0
        
        # Calcular los pasos enteros equivalentes
        steps_pos = int(round(self.dt_pos / self.dt_sim))
        steps_speed = int(round(self.dt_speed / self.dt_sim))
        steps_current = int(round(self.dt_current / self.dt_sim))
        steps_observer = int(round(self.dt_observer / self.dt_sim))
        
        for k, t in enumerate(t_arr):
            # 1. Referencias continuas del sistema
            theta_ref = get_theta_ref(t)
            
            # Obtener theta_m actual para perturbaciones dependientes de posición
            theta_m_current = self.motor.theta_m
            
            # Soporte de compatibilidad hacia atrás si get_load_torque no acepta 2 argumentos
            try:
                T_L = get_load_torque(t, theta_m_current)
            except TypeError:
                T_L = get_load_torque(t)
            
            # 2. Corrientes físicas del motor paso a paso (fases A y B)
            i_d_real = self.motor.i_d
            i_q_real = self.motor.i_q
            theta_e_real = self.motor.theta_e
            
            # Park inversa usando el theta_e real para obtener corrientes físicas de fase
            cos_t = np.cos(theta_e_real)
            sin_t = np.sin(theta_e_real)
            i_a = i_d_real * cos_t - i_q_real * sin_t
            i_b = i_d_real * sin_t + i_q_real * cos_t
            
            # 3. Lazo Discreto del Observador (LOB y EKF+HFI+DOB síncronos)
            if k % steps_observer == 0:
                if self.config.enable_lob:
                    T_L_est = self.lob.step(self.motor.torque_e, self.motor.omega_m)
                x_est, d_hat_dob, amp_hfi_val, theta_hfi_e_val, threshold_hfi_val, K_w_val, K_tx_val = self.ekf_hfi.step(U=[v_d_cmd, v_q_cmd], X=[i_a, i_b])
                hfi_measurement_valid = amp_hfi_val >= threshold_hfi_val
                theta_e_hfi_kf, omega_e_hfi_kf, alpha_e_hfi_kf, omega_m_hfi_kf = self.hfi_speed_kf.step(
                    theta_hfi_e_val, measurement_valid=hfi_measurement_valid
                )
                
            # =========================================================================
            # PANEL DE CONEXIONES Y PARÁMETROS (ESTILO SIMULINK)
            # Puedes modificar libremente cualquier señal de entrada o ganancia aquí.
            # =========================================================================
            
            # A. SEÑALES DISPONIBLES EN EL SISTEMA
            signals_avail = {
                'theta_m_real': self.motor.theta_m,
                'Wm_real': self.motor.omega_m,
                'theta_e_real': self.motor.theta_e,
                'theta_m_ekf': x_est[3],
                'Wm_ekf': x_est[2],
                'theta_e_ekf': wrapPi(x_est[3] * par.P),
                'omega_ref_ext': get_omega_ref(t) if get_omega_ref else 0.0,
                'Wm_hfi_kf': omega_m_hfi_kf,
                'theta_e_hfi_kf': theta_e_hfi_kf,
                'theta_m_hfi_kf': theta_e_hfi_kf / par.P,
            }
            
            # B. EJECUCIÓN SECUENCIAL DE LOS BLOQUES DE CONTROL (Cascada)
            
            # --- Bloque 1: PID Posición ---
            if k % steps_pos == 0:
                pos_theta_target = theta_ref
                pos_theta_m      = signals_avail.get(self.signals_routing.get('pos_theta_m', 'theta_m_real'), signals_avail['theta_m_real'])
                pos_Wm           = signals_avail.get(self.signals_routing.get('pos_Wm', 'Wm_real'), signals_avail['Wm_real'])
                pos_Ts           = self.dt_pos
                pos_Kp           = par.Kp_pos
                pos_Ki           = par.Ki_pos
                pos_Kd           = par.Kd_pos
                
                # Para ajustar dinámicamente Wtraj_max según la velocidad nominal
                pos_Wtraj_max    = abs(get_omega_ref(t)) if get_omega_ref else 12.0
                
                self.omega_ref_held, _, _ = self.pos_controller.step(
                    theta_target=pos_theta_target,
                    theta_m=pos_theta_m,
                    Wm=pos_Wm,
                    Ts=pos_Ts,
                    Kp_pos=pos_Kp,
                    Ki_pos=pos_Ki,
                    Kd_pos=pos_Kd,
                    Wtraj_max=pos_Wtraj_max
                )
            
            # Actualizar señales disponibles con la salida fresca del lazo de posición
            signals_avail['Wm_ref_from_pos'] = self.omega_ref_held
            
            # Feedforward LUT de cogging adaptativo
            feedforward_theta_m = signals_avail.get(self.signals_routing.get('pos_theta_m', 'theta_m_real'), signals_avail['theta_m_real'])
            feedforward_omega_m = signals_avail.get(self.signals_routing.get('speed_Wm', 'Wm_real'), signals_avail['Wm_real'])
            
            # El LUT aprende usando d_hat_dob. Retorna > 0 solo después de aprender.
            dob_for_control = d_hat_dob if self.config.enable_dob else 0.0
            if self.config.enable_lms:
                T_cogging_est = self.cogging_lut.step(
                    feedforward_theta_m, feedforward_omega_m, dob_for_control
                )
            else:
                T_cogging_est = 0.0
            
            # Tx es el residuo lento. Tcomp será T_cogging_est + Tx.
            Tx = dob_for_control
            
            # --- Bloque 2: PID Velocidad (PI + RI + MTPA) ---
            if k % steps_speed == 0:
                # Se lee la referencia de velocidad, por defecto conectada a la salida del lazo de posición
                speed_Wm_ref     = signals_avail.get(self.signals_routing.get('speed_Wm_ref', 'Wm_ref_from_pos'), signals_avail['Wm_ref_from_pos'])
                self.speed_Wm_ref_held = speed_Wm_ref
                speed_Wm         = signals_avail.get(self.signals_routing.get('speed_Wm', 'Wm_real'), signals_avail['Wm_real'])
                speed_Ts         = self.dt_speed
                speed_Kp         = par.Kp_speed
                speed_Ki         = par.Ki_speed
                speed_Kd         = par.Kd_speed
                speed_Tf         = 10.0 / 1000.0
                speed_Tl_tdm     = Tx + T_cogging_est
                speed_B          = par.B_var
                speed_Imax       = par.i_max
                speed_RI_enable  = 0
                speed_a_ri       = 0.0
                speed_b_ri       = 0.0
                speed_c_ri       = 0.0
                speed_d_ri       = 0.0
                speed_Kr_norm    = 0.0
                speed_Ld         = par.L_d
                speed_Lq         = par.L_q
                speed_Ke         = par.Ke
                speed_P          = par.P
                
                self.i_q_ref_held, self.i_d_ref_held, _, _, _ = self.speed_controller.step(
                    Wm_ref=speed_Wm_ref,
                    Wm=speed_Wm,
                    Ts=speed_Ts,
                    Kp_w=speed_Kp,
                    Ki_w=speed_Ki,
                    Kd_w=speed_Kd,
                    Tf_w=speed_Tf,
                    Tl_tdm=speed_Tl_tdm,
                    B=speed_B,
                    Imax=speed_Imax,
                    RI_enable=speed_RI_enable,
                    a_ri=speed_a_ri,
                    b_ri=speed_b_ri,
                    c_ri=speed_c_ri,
                    d_ri=speed_d_ri,
                    Kr_norm=speed_Kr_norm,
                    Ld=speed_Ld,
                    Lq=speed_Lq,
                    Ke=speed_Ke,
                    P=speed_P
                )
                
            # Actualizar señales disponibles con salidas frescas del lazo de velocidad
            signals_avail['Iq_ref_from_speed'] = self.i_q_ref_held
            signals_avail['Id_ref_from_speed'] = self.i_d_ref_held
            
            # --- Bloque 3: PID Corrientes dq (PI + Desacople + HFI) ---
            if k % steps_current == 0:
                curr_Wm          = signals_avail.get(self.signals_routing.get('curr_Wm', 'Wm_real'), signals_avail['Wm_real'])
                curr_X           = [i_a, i_b]
                curr_Theta_e     = signals_avail.get(self.signals_routing.get('curr_Theta_e', 'theta_e_real'), signals_avail['theta_e_real'])
                curr_Ts          = self.dt_current
                curr_Vdc         = par.V_dc
                curr_Kp_d        = par.Kp_id
                curr_Ki_d        = par.Ki_id
                curr_Kp_q        = par.Kp_iq
                curr_Ki_q        = par.Ki_iq
                curr_Psi         = par.psi_m
                curr_Ld          = par.L_d
                curr_Lq          = par.L_q
                curr_P           = par.P
                curr_Iq_ref_req  = signals_avail.get(self.signals_routing.get('curr_Iq_ref', 'Iq_ref_from_speed'), signals_avail['Iq_ref_from_speed'])
                curr_Amp_HFI     = par.Amplitud_HFI
                curr_wh          = par.wh
                curr_HFI_enable  = int(self.config.enable_hfi)
                curr_Id_ref      = signals_avail.get(self.signals_routing.get('curr_Id_ref', 'Id_ref_from_speed'), signals_avail['Id_ref_from_speed'])
                
                # Compensación del torque de cogging ideal DESACTIVADA
                real_theta_m = signals_avail.get(self.signals_routing.get('pos_theta_m', 'theta_m_real'), signals_avail['theta_m_real'])
                T_cogging_real = par.Tdm * np.sin(par.N_steps * real_theta_m + par.Phi)
                curr_Iq_comp = 0.0 # Se desactiva la trampa ideal para que el LMS y DOB trabajen
                
                curr_Kt          = par.Kt
                curr_I_max       = par.i_max
                
                Vd, Vq, _, _, Valpha, Vbetha, _, _, _, _ = self.current_controller.step(
                    Wm=curr_Wm,
                    X=curr_X,
                    Theta_e=curr_Theta_e,
                    Ts=curr_Ts,
                    Vdc=curr_Vdc,
                    Kp_d=curr_Kp_d,
                    Ki_d=curr_Ki_d,
                    Kp_q=curr_Kp_q,
                    Ki_q=curr_Ki_q,
                    Psi=curr_Psi,
                    Ld=curr_Ld,
                    Lq=curr_Lq,
                    P=curr_P,
                    Iq_ref_req=curr_Iq_ref_req,
                    Amplitud_HFI=curr_Amp_HFI,
                    wh=curr_wh,
                    HFI_enable=curr_HFI_enable,
                    Id_ref=curr_Id_ref,
                    Iq_comp=curr_Iq_comp,
                    Kt=curr_Kt,
                    I_max=curr_I_max
                )
                v_a_cmd, v_b_cmd = Valpha, Vbetha
                v_d_cmd, v_q_cmd = Vd, Vq
            # 7. Inversor Bifásico (tiempo continuo dt_sim con consignas ZOH)
            v_a, v_b = self.inverter.step(v_a_cmd, v_b_cmd, t, mode=self.inverter_mode)
            
            # Transformar voltajes del inversor al plano d-q físico real del motor
            v_d, v_q = park(v_a, v_b, theta_e_real)
            
            # 8. Integración física del motor paso a paso (dt_sim)
            self.motor.step_rk4(v_d, v_q, T_L, self.dt_sim)
            
            # 9. Registro de datos para telemetría
            history['theta_ref'][k] = theta_ref
            history['theta_m'][k] = self.motor.theta_m
            history['theta_m_est'][k] = x_est[3]
            history['omega_ref'][k] = getattr(self, 'speed_Wm_ref_held', self.omega_ref_held)
            history['omega_m'][k] = self.motor.omega_m
            history['omega_m_est'][k] = x_est[2]
            history['theta_e'][k] = wrapPi(theta_e_real)
            history['theta_e_est'][k] = wrapPi(x_est[3] * par.P)
            history['theta_hfi_e'][k] = theta_hfi_e_val
            history['theta_e_hfi_kf'][k] = theta_e_hfi_kf
            history['omega_e_hfi_kf'][k] = omega_e_hfi_kf
            history['alpha_e_hfi_kf'][k] = alpha_e_hfi_kf
            history['omega_m_hfi_kf'][k] = omega_m_hfi_kf
            history['innovation_hfi_kf'][k] = self.hfi_speed_kf.innovation

            if self.config.observer == "papers":
                history['theta_m_observer'][k] = theta_e_hfi_kf / par.P
                history['omega_m_observer'][k] = omega_m_hfi_kf
                history['theta_e_observer'][k] = theta_e_hfi_kf
            elif self.config.observer in ("ekf", "ekf_hfi"):
                history['theta_m_observer'][k] = x_est[3]
                history['omega_m_observer'][k] = x_est[2]
                history['theta_e_observer'][k] = wrapPi(x_est[3] * par.P)
            else:
                history['theta_m_observer'][k] = self.motor.theta_m
                history['omega_m_observer'][k] = self.motor.omega_m
                history['theta_e_observer'][k] = wrapPi(self.motor.theta_e)
            history['amp_hfi'][k] = amp_hfi_val
            history['threshold_hfi'][k] = threshold_hfi_val
            history['hfi_valid'][k] = 1.0 if amp_hfi_val >= threshold_hfi_val else 0.0
            history['K_w'][k] = K_w_val
            history['K_tx'][k] = K_tx_val
            
            # Cálculo de errores
            history['e_w'][k] = x_est[2] - self.motor.omega_m
            history['e_theta_hfi'][k] = wrapPi(theta_hfi_e_val - theta_e_real)
            history['e_theta_ekf'][k] = wrapPi(x_est[3] * par.P - theta_e_real)
            
            history['i_a'][k] = i_a
            history['i_b'][k] = i_b
            
            # Registrar componentes BPF
            history['i_a_fund'][k] = self.ekf_hfi.Ia_fund
            history['i_b_fund'][k] = self.ekf_hfi.Ib_fund
            history['i_a_hfi'][k] = i_a - self.ekf_hfi.Ia_fund
            history['i_b_hfi'][k] = i_b - self.ekf_hfi.Ib_fund
            
            # Cálculo de Id e Iq medidas usando el ángulo estimado (para comparación justa con EKF)
            theta_e_est_unwrapped = x_est[3] * par.P
            cos_est = np.cos(theta_e_est_unwrapped)
            sin_est = np.sin(theta_e_est_unwrapped)
            history['i_d'][k] = i_a * cos_est + i_b * sin_est
            history['i_q'][k] = -i_a * sin_est + i_b * cos_est
            
            history['i_d_ref'][k] = self.i_d_ref_held
            history['i_q_ref'][k] = self.i_q_ref_held
            history['i_d_est'][k] = x_est[0]
            history['i_q_est'][k] = x_est[1]
            
            # Reconstrucción de Ia_hat e Ib_hat a partir de Id_est e Iq_est usando el mismo ángulo
            history['i_a_est'][k] = x_est[0] * cos_est - x_est[1] * sin_est
            history['i_b_est'][k] = x_est[0] * sin_est + x_est[1] * cos_est
            
            history['torque_e'][k] = self.motor.torque_e
            history['torque_l'][k] = T_L
            history['torque_l_est'][k] = T_L_est
            history['v_d'][k] = v_d_cmd
            history['v_q'][k] = v_q_cmd
            history['d_hat_dob'][k] = d_hat_dob
            history['T_cogging_est'][k] = T_cogging_est
            history['T_cogging_real'][k] = par.Tdm * np.sin(par.N_steps * self.motor.theta_m + par.Phi)
            
        # 10. Cálculo de RMSE final
        rmse_w = np.sqrt(np.mean(history['e_w']**2))
        rmse_w_hfi_kf = np.sqrt(np.mean((history['omega_m_hfi_kf'] - history['omega_m'])**2))
        rmse_w_observer = np.sqrt(np.mean((history['omega_m_observer'] - history['omega_m'])**2))
        e_theta_observer = wrapPi(history['theta_e_observer'] - history['theta_e'])
        rmse_theta_observer = np.sqrt(np.mean(e_theta_observer**2))
        
        # e_theta_ekf ya está calculado usando wrapPi(theta_e_est - theta_e_real) en el bucle
        rmse_theta = np.sqrt(np.mean(history['e_theta_ekf']**2))
        rmse_theta_deg = np.degrees(rmse_theta)
        rmse_theta_m = rmse_theta / par.P
        
        rmse_tx = np.sqrt(np.mean((history['d_hat_dob'] - history['torque_l'])**2))
        hfi_activation_ratio = np.mean(history['hfi_valid'])
        
        history['rmse_w'] = rmse_w
        history['rmse_w_hfi_kf'] = rmse_w_hfi_kf
        history['rmse_w_observer'] = rmse_w_observer
        history['rmse_theta_observer'] = rmse_theta_observer
        history['observer_name'] = self.config.observer
        history['config'] = self.config
        history['rmse_theta'] = rmse_theta
        history['rmse_theta_deg'] = rmse_theta_deg
        history['rmse_theta_m'] = rmse_theta_m
        history['rmse_tx'] = rmse_tx
        history['hfi_activation_ratio'] = hfi_activation_ratio
        
        print(f"\n[Métricas de Estimación]")
        print(f" - RMSE Wm       : {rmse_w:.4f} rad/s")
        print(f" - RMSE Wm HFI-KF: {rmse_w_hfi_kf:.4f} rad/s")
        print(f" - RMSE seleccionado: {rmse_w_observer:.4f} rad/s ({self.config.observer})")
        print(f" - RMSE Theta_e  : {rmse_theta:.4f} rad_e ({rmse_theta_deg:.2f}°_e)")
        print(f" - RMSE Theta_m  : {rmse_theta_m:.4f} rad_m")
        print(f" - RMSE Tx       : {rmse_tx:.4f} Nm")
        print(f" - Activación HFI: {hfi_activation_ratio*100:.1f} %")
            
        return history
