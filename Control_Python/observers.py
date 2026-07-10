import numpy as np

def wrapPi(x):
    """
    Acota ángulos periódicos entre -pi y pi.
    Equivalente a wrap_pi del bloque de MATLAB.
    """
    return np.arctan2(np.sin(x), np.cos(x))


class SlidingModeObserver:
    """
    Observador de Modos Deslizantes (SMO) en Tiempo Discreto 
    para la estimación de la posición y velocidad eléctrica de un Motor Stepper Bifásico.
    Incluye un PLL (Phase-Locked Loop) discreto integrado para estimación libre de retraso.
    Las entradas en coordenadas de estator alfa-beta corresponden directamente a las fase A y B.
    """
    def __init__(self, R_s, L_s, psi_m, P, dt_obs, l_gain=120.0, f_cutoff=200.0, omega_pll=150.0):
        self.R_s = R_s          # Resistencia de estator (Ohm)
        self.L_s = L_s          # Inductancia de estator promedio (H)
        self.psi_m = psi_m      # Flujo de imanes (Wb)
        self.P = P              # Pares de polos (dientes del rotor)
        self.dt = dt_obs        # Periodo de muestreo del observador (s)
        
        # Parámetros del SMO
        self.l = l_gain         # Ganancia del modo deslizante
        self.f_c = f_cutoff     # Frecuencia de corte del filtro pasa bajos (Hz)
        self.tau = 1.0 / (2.0 * np.pi * self.f_c)  # Constante de tiempo del filtro
        
        # Estados internos del SMO
        self.i_alpha_est = 0.0
        self.i_beta_est = 0.0
        self.e_alpha_est = 0.0  # Back-EMF filtrada alfa (fase A)
        self.e_beta_est = 0.0   # Back-EMF filtrada beta (fase B)
        
        # Parámetros del PLL
        zeta = 1.0  # Críticamente amortiguado
        self.K_p_pll = 2.0 * zeta * omega_pll
        self.K_i_pll = omega_pll ** 2
        
        # Estados internos del PLL
        self.theta_e_est = 0.0    # Posición eléctrica estimada (rad)
        self.omega_e_est = 0.0    # Velocidad eléctrica estimada (rad/s)
        self.pll_integrator = 0.0

    def step(self, v_alpha, v_beta, i_alpha, i_beta):
        """
        Ejecuta un paso de estimación del SMO + PLL.
        """
        err_i_alpha = self.i_alpha_est - i_alpha
        err_i_beta = self.i_beta_est - i_beta
        
        z_alpha = self.l * np.sign(err_i_alpha)
        z_beta = self.l * np.sign(err_i_beta)
        
        di_alpha = (-self.R_s * self.i_alpha_est + v_alpha - z_alpha) / self.L_s
        di_beta = (-self.R_s * self.i_beta_est + v_beta - z_beta) / self.L_s
        
        self.i_alpha_est += self.dt * di_alpha
        self.i_beta_est += self.dt * di_beta
        
        de_alpha = (z_alpha - self.e_alpha_est) / self.tau
        de_beta = (z_beta - self.e_beta_est) / self.tau
        
        self.e_alpha_est += self.dt * de_alpha
        self.e_beta_est += self.dt * de_beta
        
        epsilon = -self.e_alpha_est * np.cos(self.theta_e_est) - self.e_beta_est * np.sin(self.theta_e_est)
        
        self.pll_integrator += self.K_i_pll * self.dt * epsilon
        self.omega_e_est = self.K_p_pll * epsilon + self.pll_integrator
        
        self.theta_e_est += self.dt * self.omega_e_est
        self.theta_e_est = self.theta_e_est % (2.0 * np.pi)
        
        omega_m_est = self.omega_e_est / self.P
        
        return self.theta_e_est, omega_m_est


class LoadTorqueObserver:
    """
    Observador de Torque de Carga (Luenberger Observer) en tiempo discreto.
    """
    def __init__(self, J, B, dt_obs, p_obs=30.0):
        self.J = J              # Inercia del rotor (kg*m^2)
        self.B = B              # Fricción viscosa (N*m*s/rad)
        self.dt = dt_obs        # Periodo de muestreo (s)
        
        self.L1 = 2.0 * p_obs - (self.B / self.J)
        self.L2 = self.J * (p_obs ** 2)
        
        self.omega_m_est = 0.0
        self.T_L_est = 0.0

    def step(self, torque_e, omega_m):
        err_omega = omega_m - self.omega_m_est
        d_omega = (torque_e - self.T_L_est - self.B * self.omega_m_est) / self.J + self.L1 * err_omega
        d_TL = -self.L2 * err_omega
        
        self.omega_m_est += self.dt * d_omega
        self.T_L_est += self.dt * d_TL
        
        return self.T_L_est


class EKF_HFI_Observer:
    """
    Filtro de Kalman Extendido (EKF) con Inyección de Alta Frecuencia (HFI)
    y Observador de Perturbaciones (DOB) integrado.
    Alineado 100% con la S-Function ekf_dq_4states_final de tu modelo de Simulink.
    """
    def __init__(self, Ke, Kt, R, J, Ld, Lq, Nr, N_steps, B, Tdtm, beta0, beta1, beta2, wh, b0, b1, b2, a1, a2, Ts, use_bpf=True):
        self.Ke = Ke
        self.Kt = Kt
        self.R = R
        self.J = J
        self.Ld = Ld
        self.Lq = Lq
        self.Nr = Nr            # Pares de polos (dientes del rotor)
        self.N_steps = N_steps  # Pasos por revolución (200)
        self.B = B
        self.Tdtm = Tdtm
        
        # Ganancias DOB
        self.beta0 = beta0
        self.beta1 = beta1
        self.beta2 = beta2
        
        # Parámetros HFI BPF y frecuencia
        self.wh = wh
        self.b0 = b0
        self.b1 = b1
        self.b2 = b2
        self.a1 = a1
        self.a2 = a2
        
        self.Ts = Ts
        
        # Inicialización de variables persistentes (estados internos)
        self.x_hat = np.zeros(5)
        self.P = np.diag([1e-4, 1e-4, 1e-2, 1e-6, 1e-4])
        self.d_do = 0.0
        self.theta_do = 0.0
        self.omega_do = 0.0
        self.k = 0
        self.hfi_phase_acc = 0.0
        self.use_bpf = use_bpf
        
        # Filtros de banda de inyección HFI (BPF)
        self.xa1 = 0.0
        self.xa2 = 0.0
        self.ya1 = 0.0
        self.ya2 = 0.0
        
        self.xb1 = 0.0
        self.xb2 = 0.0
        self.yb1 = 0.0
        self.yb2 = 0.0
        
        self.t_h = 0.0
        self.eps_f = 0.0
        self.amp_hfi = 0.0
        self.theta_hfi_e = 0.0
        
        self.demod_s = 0.0
        self.demod_c = 0.0
        self.amp_noise_mean = 0.005
        self.amp_noise_var = 1e-6
        
        self.Ia_fund = 0.0
        self.Ib_fund = 0.0
        
        self.theta_hfi_int = 0.0
        self.omega_hfi_int = 0.0
        self.use_hfi = True

    def step(self, U, X):
        """
        Ejecuta un paso de estimación EKF + HFI + DOB.
        
        Parámetros:
        - U: Voltajes de fase fundamentales aplicados [Vd_fund, Vq_fund] (V)
        - X: Corrientes de fase medidas [i_alpha, i_beta] (A)
        
        Retorna:
        - x_out: Vector de estados estimados [i_d, i_q, omega_m, theta_m, Tx]
        - d_hat: Perturbación de torque estimada por el DOB (N*m)
        - amp_hfi: Amplitud detectada de la inyección HFI
        - theta_hfi_e: Ángulo eléctrico estimado por HFI
        """
        Ia_raw = X[0]
        Ib_raw = X[1]
        
        self.k += 1
        
        # 1. Filtro BPF para extracción quirúrgica de HFI (mejor que LPF puro para evitar phase-lag)
        if self.use_bpf:
            Iah = self.b0 * Ia_raw + self.b1 * self.xa1 + self.b2 * self.xa2 - self.a1 * self.ya1 - self.a2 * self.ya2
            self.xa2 = self.xa1; self.xa1 = Ia_raw; self.ya2 = self.ya1; self.ya1 = Iah
    
            Ibh = self.b0 * Ib_raw + self.b1 * self.xb1 + self.b2 * self.xb2 - self.a1 * self.yb1 - self.a2 * self.yb2
            self.xb2 = self.xb1; self.xb1 = Ib_raw; self.yb2 = self.yb1; self.yb1 = Ibh
        else:
            Iah = 0.0
            Ibh = 0.0
            
        # Componente HFI
        Ia_hfi = Iah
        Ib_hfi = Ibh
        
        # Componente fundamental (sin el ruido de inyección)
        self.Ia_fund = Ia_raw - Iah
        self.Ib_fund = Ib_raw - Ibh
        
        Ia_ekf = self.Ia_fund
        Ib_ekf = self.Ib_fund

        # 2. Procesamiento HFI
        theta_e_ekf_prev = wrapPi(self.Nr * self.x_hat[3])
        ch_prev = np.cos(theta_e_ekf_prev)
        sh_prev = np.sin(theta_e_ekf_prev)
        
        Idh_prev = Ia_hfi * ch_prev + Ib_hfi * sh_prev
        Iqh_prev = -Ia_hfi * sh_prev + Ib_hfi * ch_prev
        
        self.t_h += self.Ts
        
        # Demodulación I/Q ortogonal para amplitud y fase
        demod_s_raw = Idh_prev * np.sin(self.wh * self.t_h)
        demod_c_raw = Idh_prev * np.cos(self.wh * self.t_h)
        
        # LPF para demodulación (I/Q)
        f_lpf_iq = 50.0
        alpha_iq = self.Ts / (self.Ts + 1.0 / (2.0 * np.pi * f_lpf_iq))
        self.demod_s += alpha_iq * (demod_s_raw - self.demod_s)
        self.demod_c += alpha_iq * (demod_c_raw - self.demod_c)
        
        self.amp_hfi = 2.0 * np.sqrt(self.demod_s**2 + self.demod_c**2)
        
        # Demodulación para el error angular (eje q)
        eps_raw = Iqh_prev * np.sin(self.wh * self.t_h)
        f_lpf_hfi = 50.0
        alpha_lpf = self.Ts / (self.Ts + 1.0 / (2.0 * np.pi * f_lpf_hfi))
        self.eps_f += alpha_lpf * (eps_raw - self.eps_f)

        # Actualización de umbral de ruido de HFI cuando está inactivo
        if not self.use_hfi:
            alpha_th = self.Ts / (self.Ts + 1.0 / (2.0 * np.pi * 1.0)) # 1 Hz
            self.amp_noise_mean += alpha_th * (self.amp_hfi - self.amp_noise_mean)
            self.amp_noise_var += alpha_th * ((self.amp_hfi - self.amp_noise_mean)**2 - self.amp_noise_var)
            
        threshold_hfi = self.amp_noise_mean + 3.0 * np.sqrt(self.amp_noise_var)

        # Normalización del error HFI
        DeltaL = self.Lq - self.Ld
        Lavg = 0.5 * (self.Ld + self.Lq)
        if abs(DeltaL) < 1e-12:
            eps_norm = 0.0
        else:
            saliency = DeltaL / Lavg
            eps_norm = self.eps_f / abs(saliency)

        K_hfi_pll = 5e-3
        theta_delay = 0.0
        # CORRECCIÓN DE SIGNO: eps_f es positivo cuando theta_ekf atrasa a theta_real.
        # Por lo tanto, debemos SUMAR (K_hfi_pll * eps_norm) para corregir el atraso.
        theta_hfi_raw = wrapPi(theta_e_ekf_prev + K_hfi_pll * eps_norm + theta_delay)
        
        # 3. Detección de ambigüedad de polaridad (N/S)
        if abs(wrapPi(theta_hfi_raw - theta_e_ekf_prev)) > np.pi / 2.0:
            theta_hfi_raw = wrapPi(theta_hfi_raw + np.pi)
            
        self.theta_hfi_e = theta_hfi_raw
        
        # 4. Voltajes a dq (ya vienen en dq fundamental, no requieren Park!)
        Vd_ekf = U[0]
        Vq_ekf = U[1]

        # 5. EKF - Predicción
        Id = self.x_hat[0]
        Iq = self.x_hat[1]
        Wm = self.x_hat[2]
        Tx = self.x_hat[4]
        
        dId = (Vd_ekf - self.R * Id + self.Nr * Wm * self.Lq * Iq) / self.Ld
        dIq = (Vq_ekf - self.R * Iq - self.Nr * Wm * self.Ld * Id - self.Ke * Wm) / self.Lq

        Te = self.Kt * Iq - self.Nr * (self.Ld - self.Lq) * Id * Iq
        dWm = (Te - Tx - self.B * Wm) / self.J
        dTh = Wm
        dTx = 0.0

        x_pred = self.x_hat + self.Ts * np.array([dId, dIq, dWm, dTh, dTx])

        # Jacobiano
        A = np.zeros((5, 5))
        A[0, 0] = -self.R / self.Ld
        A[0, 1] = self.Nr * Wm * self.Lq / self.Ld
        A[0, 2] = self.Nr * self.Lq * Iq / self.Ld
        
        A[1, 0] = -self.Nr * Wm * self.Ld / self.Lq
        A[1, 1] = -self.R / self.Lq
        A[1, 2] = -(self.Nr * self.Ld * Id + self.Ke) / self.Lq
        
        A[2, 0] = -self.Nr * (self.Ld - self.Lq) * Iq / self.J
        A[2, 1] = (self.Kt - self.Nr * (self.Ld - self.Lq) * Id) / self.J
        A[2, 2] = -self.B / self.J
        A[2, 4] = -1.0 / self.J
        
        A[3, 2] = 1.0

        Phi = np.eye(5) + A * self.Ts
        
        # Qk Override si se establecen externamente para los experimentos
        if hasattr(self, 'Qk_override') and self.Qk_override is not None:
            Qk = self.Qk_override
        else:
            Qk = np.diag([1e-4, 1e-4, 1e-7, 1e-10, 1e-8])
            
        P_pred = Phi @ self.P @ Phi.T + Qk

        # 6. EKF - Primera Corrección (Corrientes Ia, Ib fundamentales)
        theta_e_pred = self.Nr * x_pred[3]
        c_pred = np.cos(theta_e_pred)
        s_pred = np.sin(theta_e_pred)
        
        Id_p = x_pred[0]
        Iq_p = x_pred[1]
        
        Ia_pred = Id_p * c_pred - Iq_p * s_pred
        Ib_pred = Id_p * s_pred + Iq_p * c_pred
        h_curr = np.array([Ia_pred, Ib_pred])
        
        Ck = np.array([
            [c_pred, -s_pred, 0.0, self.Nr * (-Id_p * s_pred - Iq_p * c_pred), 0.0],
            [s_pred,  c_pred, 0.0, self.Nr * (Id_p * c_pred - Iq_p * s_pred),  0.0]
        ])
        
        Rk_curr = np.diag([1e-3, 1e-3])
        z_curr = np.array([Ia_ekf, Ib_ekf])
        innov_curr = z_curr - h_curr
        
        Sk_curr = Ck @ P_pred @ Ck.T + Rk_curr
        K_curr = (P_pred @ Ck.T) @ np.linalg.solve(Sk_curr, np.eye(2))
        
        x_corr1 = x_pred + K_curr @ innov_curr
        
        I5 = np.eye(5)
        I_KC = I5 - K_curr @ Ck
        P_corr1 = I_KC @ P_pred @ I_KC.T + K_curr @ Rk_curr @ K_curr.T
        P_corr1 = 0.5 * (P_corr1 + P_corr1.T)

        # Inicializar K_hfi en caso de que no se use
        K_hfi = np.zeros((5, 1))

        # 7. EKF - Segunda Corrección (HFI adaptativo)
        hfi_valid = self.amp_hfi >= threshold_hfi
        
        force_hfi_off = hasattr(self, 'force_hfi_off') and self.force_hfi_off
        
        # Acumulador de fase para trigger HFI adaptativo
        self.hfi_phase_acc += self.wh * self.Ts
        hfi_trigger = False
        if self.hfi_phase_acc >= 2.0 * np.pi:
            self.hfi_phase_acc -= 2.0 * np.pi
            hfi_trigger = True
            
        if not force_hfi_off and hfi_valid and hfi_trigger:
            self.use_hfi = True
        else:
            self.use_hfi = False

        if self.use_hfi:
            H_hfi = np.array([[0.0, 0.0, 0.0, self.Nr, 0.0]])
            innov_hfi = np.array([wrapPi(self.theta_hfi_e - self.Nr * x_corr1[3])])
            
            # Mantener R_hfi según instrucciones
            sigma_theta_hfi = 0.15
            
            R_hfi = np.array([[sigma_theta_hfi**2]])
            
            Sk_hfi = H_hfi @ P_corr1 @ H_hfi.T + R_hfi
            K_hfi = (P_corr1 @ H_hfi.T) / Sk_hfi[0, 0]
            
            x_corr2 = x_corr1 + K_hfi @ innov_hfi
            
            I_KH = I5 - K_hfi @ H_hfi
            P_corr2 = I_KH @ P_corr1 @ I_KH.T + K_hfi @ R_hfi @ K_hfi.T
            P_corr2 = 0.5 * (P_corr2 + P_corr2.T)
        else:
            x_corr2 = x_corr1
            P_corr2 = P_corr1

        # 8. Protecciones (Anti-Windup y Anti-NaN)
        # Limitar la velocidad y el torque para evitar explosiones
        Wm_max = 500.0
        Tx_max = 10.0
        if x_corr2[2] > Wm_max: x_corr2[2] = Wm_max
        if x_corr2[2] < -Wm_max: x_corr2[2] = -Wm_max
        if x_corr2[4] > Tx_max: x_corr2[4] = Tx_max
        if x_corr2[4] < -Tx_max: x_corr2[4] = -Tx_max
        
        # Si algo se volvió NaN o Infinito, resetear el filtro
        if not np.all(np.isfinite(x_corr2)) or not np.all(np.isfinite(P_corr2)):
            print(f"[EKF Warning] NaN o Inf detectado! Reseteando filtro a t={self.t_h:.3f}s")
            self.reset()
            x_corr2 = self.x_hat
            P_corr2 = self.P
        else:
            self.x_hat = x_corr2
            self.P = P_corr2

        d_hat = self.x_hat[4]
        x_out = self.x_hat.copy()
        
        K_w = np.linalg.norm(K_curr[2, :]) + (abs(K_hfi[2, 0]) if self.use_hfi else 0.0)
        K_tx = np.linalg.norm(K_curr[4, :]) + (abs(K_hfi[4, 0]) if self.use_hfi else 0.0)
        
        return x_out, d_hat, self.amp_hfi, self.theta_hfi_e, threshold_hfi, K_w, K_tx

    def reset(self):
        self.x_hat = np.zeros(5)
        self.P = np.diag([1e-4, 1e-4, 1e-2, 1e-6, 1e-4])
        self.Ia_fund = 0.0
        self.Ib_fund = 0.0
        self.t_h = 0.0
        self.eps_f = 0.0
        self.amp_hfi = 0.0
        self.theta_hfi_e = 0.0
        self.use_hfi = True
