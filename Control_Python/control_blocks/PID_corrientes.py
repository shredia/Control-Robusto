import numpy as np

class PIDCorrientes:
    """
    Lazo de Corrientes dq: PID_Idq(Wm, X, Theta_e, Ts, Vdc, Kp_d, Ki_d, Kp_q, Ki_q, Psi, Ld, Lq, P, Iq_ref_req, Amplitud_HFI, wh, HFI_enable, Id_ref, Iq_comp, Kt, I_max).
    Basado en la S-Function de Simulink PI_dq_vectorial_fixed.
    """
    def __init__(self):
        # Variables persistentes (estados del bloque)
        self.v_int_d = 0.0
        self.v_int_q = 0.0
        self.t_hfi = 0.0

    def step(self, Wm, X, Theta_e, Ts, Vdc, Kp_d, Ki_d, Kp_q, Ki_q, Psi, Ld, Lq, P, Iq_ref_req, Amplitud_HFI, wh, HFI_enable, Id_ref, Iq_comp, Kt, I_max):
        """
        Ejecuta un paso discreto del lazo de corrientes dq con desacoplamiento y HFI.
        
        Parámetros:
        - Wm: Velocidad mecánica realimentada (rad/s)
        - X: Vector de corrientes estatóricas bifásicas [Ia_real, Ib_real] (A)
        - Theta_e: Ángulo eléctrico realimentado (rad)
        - Ts: Período de muestreo (s)
        - Vdc: Voltaje del bus DC (V)
        - Kp_d, Ki_d: Ganancias del eje d (PI)
        - Kp_q, Ki_q: Ganancias del eje q (PI)
        - Psi: Flujo de imanes nominal (Wb)
        - Ld, Lq: Inductancias dq nominales (H)
        - P: Número de pares de polos (dientes)
        - Iq_ref_req: Referencia de corriente q (A)
        - Amplitud_HFI: Voltaje de inyección de alta frecuencia (V)
        - wh: Frecuencia de inyección HFI (rad/s)
        - HFI_enable: Flag de inyección HFI (1: Habilitado, 0: Deshabilitado)
        - Id_ref: Referencia de corriente d (A) (MTPA)
        - Iq_comp: Corriente de compensación adicional (A)
        - Kt: Constante de torque (N*m/A)
        - I_max: Corriente máxima permitida (A)
        
        Retorna:
        - Vd, Vq: Tensiones de consigna fundamental d y q (V)
        - Id_out, Iq_out: Corrientes medidas transformadas en d y q (A)
        - Valpha, Vbetha: Tensiones de consigna físicas de fase A y B (V)
        - Vd_hfi, Vq_hfi: Tensiones dq incluyendo inyección HFI (V)
        - err_d, err_q: Errores de corriente en d y q (A)
        """
        # 1. Protecciones numéricas
        Ts_eff = max(Ts, 1e-9)
        Vdc_eff = max(abs(Vdc), 1e-6)
        Imax_eff = max(abs(I_max), 1e-6)

        # 2. Corrientes alpha-beta (A y B) a dq (Park)
        Ialpha = X[0]
        Ibeta = X[1]
        cos_th = np.cos(Theta_e)
        sin_th = np.sin(Theta_e)

        Id = Ialpha * cos_th + Ibeta * sin_th
        Iq = -Ialpha * sin_th + Ibeta * cos_th

        # 3. Referencias de corriente
        Id_ref_limited = max(-Imax_eff, min(Imax_eff, Id_ref))
        Iq_limit = np.sqrt(max(0.0, Imax_eff**2 - Id_ref_limited**2))
        Iq_ref_total = Iq_ref_req + Iq_comp
        Iq_ref = max(-Iq_limit, min(Iq_limit, Iq_ref_total))

        # 4. Errores de corriente
        err_d = Id_ref_limited - Id
        err_q = Iq_ref - Iq

        # Velocidad eléctrica
        We = P * Wm

        # 5. Desacople usando corrientes reales
        V_cross_d = -We * Lq * Iq
        V_cross_q = We * (Ld * Id + Psi)

        # 6. PI fundamental dq
        Vd_pre = Kp_d * err_d + self.v_int_d + V_cross_d
        Vq_pre = Kp_q * err_q + self.v_int_q + V_cross_q

        # 7. HFI y margen de tensión disponible
        self.t_hfi += Ts_eff
        if HFI_enable != 0:
            v_hfi = Amplitud_HFI * np.cos(wh * self.t_hfi)
        else:
            v_hfi = 0.0

        # Límite vectorial de tensión con SVPWM (o H-bridges en plano de Clarke)
        # Para inversores bifásicos independientes, la tensión máxima por fase es Vdc_eff
        # Sin embargo, para mantener compatibilidad exacta con tu bloque usamos Vdc_eff / sqrt(3)
        Vmax = Vdc_eff / np.sqrt(3.0)

        if HFI_enable != 0:
            Vfund_max = max(0.0, Vmax - abs(Amplitud_HFI))
        else:
            Vfund_max = Vmax

        # 8. Saturación vectorial de tensión fundamental
        Vmag_pre = np.hypot(Vd_pre, Vq_pre)
        if Vmag_pre > Vfund_max and Vmag_pre > 1e-12:
            scale_v = Vfund_max / Vmag_pre
            Vd = Vd_pre * scale_v
            Vq = Vq_pre * scale_v
        else:
            Vd = Vd_pre
            Vq = Vq_pre

        # 9. Anti-windup por back-calculation (de tu S-Function)
        Kaw_d = Ki_d / abs(Kp_d) if abs(Kp_d) > 1e-9 else 0.0
        Kaw_q = Ki_q / abs(Kp_q) if abs(Kp_q) > 1e-9 else 0.0

        self.v_int_d += Ts_eff * (Ki_d * err_d + Kaw_d * (Vd - Vd_pre))
        self.v_int_q += Ts_eff * (Ki_q * err_q + Kaw_q * (Vq - Vq_pre))

        # 10. Inyección HFI en eje d
        Vd_hfi = Vd + v_hfi
        Vq_hfi = Vq

        # Saturación final de seguridad
        Vmag_hfi = np.hypot(Vd_hfi, Vq_hfi)
        if Vmag_hfi > Vmax and Vmag_hfi > 1e-12:
            scale_hfi = Vmax / Vmag_hfi
            Vd_hfi *= scale_hfi
            Vq_hfi *= scale_hfi

        # 11. Transformada inversa de Park
        Valpha = Vd_hfi * cos_th - Vq_hfi * sin_th
        Vbetha = Vd_hfi * sin_th + Vq_hfi * cos_th

        # Salidas de monitoreo
        Id_out = Id
        Iq_out = Iq

        return Vd, Vq, Id_out, Iq_out, Valpha, Vbetha, Vd_hfi, Vq_hfi, err_d, err_q

    def reset(self):
        self.v_int_d = 0.0
        self.v_int_q = 0.0
        self.t_hfi = 0.0
