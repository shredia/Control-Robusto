import numpy as np

class PISpeedRIMTPA:
    """
    Lazo de Velocidad y MTPA: PID_Wm(Wm_ref, Wm, Ts, Kp_w, Ki_w, Kd_w, Tf_w, Tl_tdm, B, Imax, RI_enable, a_ri, b_ri, c_ri, d_ri, Kr_norm, Ld, Lq, Ke, P).
    Basado en el bloque PI_speed_RI_MTPA de Simulink.
    """
    def __init__(self):
        # Variables persistentes (estados del bloque)
        self.int_w = 0.0
        self.Wm_prev = None
        self.dWm_filt_prev = 0.0
        
        # Estados del controlador resonante (cogging ripple suppressor)
        self.y1_ri = 0.0
        self.y2_ri = 0.0
        self.e1_ri = 0.0
        self.e2_ri = 0.0

    def step(self, Wm_ref, Wm, Ts, Kp_w, Ki_w, Kd_w, Tf_w, Tl_tdm, B, Imax, RI_enable, a_ri, b_ri, c_ri, d_ri, Kr_norm, Ld, Lq, Ke, P):
        """
        Ejecuta un paso discreto del control de velocidad mecánico con MTPA y resonadores.
        
        Parámetros:
        - Wm_ref: Velocidad mecánica de referencia (rad/s)
        - Wm: Velocidad mecánica realimentada (rad/s)
        - Ts: Período de muestreo (s)
        - Kp_w, Ki_w, Kd_w: Ganancias del regulador PID
        - Tf_w: Constante de tiempo del filtro derivativo (s)
        - Tl_tdm: Torque de carga de compensación feedforward (N*m)
        - B: Fricción viscosa nominal (N*m*s/rad)
        - Imax: Límite de corriente máxima del estator (A)
        - RI_enable: Flag para activar resonador RI (1: Habilitado, 0: Deshabilitado)
        - a_ri, b_ri, c_ri, d_ri, Kr_norm: Coeficientes del resonador discreto
        - Ld, Lq: Inductancias dq nominales (H)
        - Ke: Constante de torque Kt (N*m/A)
        - P: Número de pares de polos (dientes)
        
        Retorna:
        - Iq_ref: Corriente de referencia q óptima (A)
        - Id_ref: Corriente de referencia d óptima (A)
        - Iq_cogging: Salida de corriente de compensación del resonador (A)
        - Iq_PI: Salida de corriente del bloque PI/PID (A)
        - error_Wm: Error de velocidad (rad/s)
        """
        if self.Wm_prev is None:
            self.Wm_prev = Wm

        Ts_eff = max(Ts, 1e-9)
        Imax_eff = max(abs(Imax), 1e-6)
        Kt_pm = Ke
        Krel = P * (Ld - Lq)

        # 1. Error de velocidad
        error_Wm = Wm_ref - Wm

        # 2. Derivada filtrada sobre velocidad medida
        dWm_raw = (Wm - self.Wm_prev) / Ts_eff
        if Tf_w > 0.0:
            alpha_d = Ts_eff / (Tf_w + Ts_eff)
            dWm_filt = self.dWm_filt_prev + alpha_d * (dWm_raw - self.dWm_filt_prev)
        else:
            dWm_filt = dWm_raw
        
        D_w = -Kd_w * dWm_filt

        # 3. Feedforward de carga y roce viscoso
        T_ff = Tl_tdm + B * Wm_ref
        Iq_ff = T_ff / Kt_pm if abs(Kt_pm) > 1e-9 else 0.0

        # 4. Resonador para ripple / cogging
        if RI_enable != 0:
            Iq_RI_new = c_ri * self.y1_ri - d_ri * self.y2_ri + Kr_norm * (error_Wm - a_ri * self.e1_ri + b_ri * self.e2_ri)
        else:
            Iq_RI_new = 0.0

        # 5. Integrador candidato
        int_w_candidate = self.int_w + Ki_w * Ts_eff * error_Wm

        # 6. Control de velocidad + MTPA + límite vectorial (Dos pasadas para Anti-windup condicional)
        Iq_ref = 0.0
        Id_ref = 0.0
        Iq_PI = 0.0
        
        for pass_idx in range(2):
            Iq_PI = Kp_w * error_Wm + int_w_candidate + D_w
            Iq_ref_raw = Iq_ff + Iq_PI + Iq_RI_new

            # Petición de usuario: Id_ref = 0, Iq_ref = proveniente del PID de velocidad
            Id_ref = 0.0
            Iq_ref = max(-Imax_eff, min(Imax_eff, Iq_ref_raw))
            sat = abs(Iq_ref_raw) > Imax_eff

            # Protección final de magnitud vectorial (Círculo de Corriente)
            Iref_mag = np.hypot(Id_ref, Iq_ref)
            if Iref_mag > Imax_eff:
                scale = Imax_eff / Iref_mag
                Id_ref *= scale
                Iq_ref *= scale
                sat = True

            # Anti-windup condicional del PI de velocidad
            int_increment = int_w_candidate - self.int_w
            integrator_worsens_sat = (Iq_ref_raw * int_increment) > 0
            
            if pass_idx == 0 and sat and integrator_worsens_sat:
                int_w_candidate = self.int_w
            else:
                break

        # Guarda integrador final
        self.int_w = int_w_candidate

        # 7. Actualización de memorias
        self.Wm_prev = Wm
        self.dWm_filt_prev = dWm_filt

        if RI_enable != 0:
            self.y2_ri = self.y1_ri
            self.y1_ri = Iq_RI_new
            self.e2_ri = self.e1_ri
            self.e1_ri = error_Wm
        else:
            self.y1_ri = 0.0
            self.y2_ri = 0.0
            self.e1_ri = 0.0
            self.e2_ri = 0.0

        Iq_cogging = Iq_RI_new

        return Iq_ref, Id_ref, Iq_cogging, Iq_PI, error_Wm

    def reset(self):
        self.int_w = 0.0
        self.Wm_prev = None
        self.dWm_filt_prev = 0.0
        self.y1_ri = 0.0
        self.y2_ri = 0.0
        self.e1_ri = 0.0
        self.e2_ri = 0.0
