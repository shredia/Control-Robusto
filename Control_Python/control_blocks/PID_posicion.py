import numpy as np

def wrapPi(x):
    """
    Acota ángulos periódicos entre -pi y pi.
    Equivalente a wrapToPi de MATLAB: wrapPi(x) = atan2(sin(x), cos(x))
    """
    return np.arctan2(np.sin(x), np.cos(x))

class PIDPositionHold:
    """
    Lazo de Posición Mecánico: PID_posicion(theta_target, theta_m, Wm, Ts, Kp_pos, Ki_pos, Kd_pos).
    Basado en el bloque PID_position_hold de Simulink.
    """
    def __init__(self):
        # Variables persistentes (estados del bloque)
        self.theta_ref_int = None
        self.integral_pos = 0.0

    def step(self, theta_target, theta_m, Wm, Ts, Kp_pos, Ki_pos, Kd_pos, Wtraj_max=12.0):
        """
        Ejecuta un paso discreto de control de posición con generador de trayectorias.
        
        Parámetros:
        - theta_target: Consigna / Referencia externa de posición (rad)
        - theta_m: Posición mecánica realimentada (rad)
        - Wm: Velocidad mecánica realimentada (rad/s)
        - Ts: Período de muestreo del lazo (s)
        - Kp_pos, Ki_pos, Kd_pos: Ganancias del regulador
        
        Retorna:
        - Wm_ref_pos: Velocidad de referencia generada (rad/s)
        - theta_m_ref: Salida del generador de trayectoria interna (rad)
        - e_theta: Error de posición angular acotado (rad)
        """
        if self.theta_ref_int is None:
            self.theta_ref_int = theta_m
            self.integral_pos = 0.0

        # Parámetros del bloque de Simulink
        Wmax_pos = 15.0      # Límite de velocidad de salida del lazo de posición [rad/s]
        Kaw = 20.0           # Ganancia anti-windup [1/s]

        # 1. Generador de trayectoria (mueve la referencia interna suavemente hacia theta_target)
        e_traj = wrapPi(theta_target - self.theta_ref_int)
        dtheta_max = Wtraj_max * Ts
        dtheta = min(max(e_traj, -dtheta_max), dtheta_max)
        self.theta_ref_int = wrapPi(self.theta_ref_int + dtheta)

        # Velocidad de la trayectoria generada
        Wm_traj_ref = dtheta / Ts
        theta_m_ref = self.theta_ref_int

        # 2. Error de posición angular acotado en [-pi, pi)
        e_theta = wrapPi(theta_m_ref - theta_m)

        # 3. Control PID de posición usando velocidad realimentada en vez de derivar e_theta
        Wm_unsat = Wm_traj_ref + Kp_pos * e_theta + self.integral_pos + Kd_pos * (Wm_traj_ref - Wm)

        # 4. Saturación
        Wm_ref_pos = min(max(Wm_unsat, -Wmax_pos), Wmax_pos)

        # 5. Anti-Windup por Back-Calculation
        self.integral_pos += Ts * (Ki_pos * e_theta + Kaw * (Wm_ref_pos - Wm_unsat))

        return Wm_ref_pos, theta_m_ref, e_theta

    def reset(self):
        self.theta_ref_int = None
        self.integral_pos = 0.0
