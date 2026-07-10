import numpy as np

class StepperMotor2Phase:
    """
    Modelo continuo de un Motor Paso a Paso Bifásico (Stepper Motor) con saliencia,
    que incluye torque de reluctancia y torque de cogging/detent.
    Alineado matemáticamente con la S-Function Simulink PMSM_salient_dq.m.
    """
    def __init__(self, R_s=2.5, L_d=3.66e-3, L_q=9.61e-3, psi_m=0.00326, J=1.08e-5, B=2.0e-4, P=50, Tdm=18e-3, N_steps=200, Phi=np.pi/2.0):
        # Parámetros del motor (Reales de la planta física)
        self.R_s = R_s        # Resistencia de estator (Ohm)
        self.L_d = L_d        # Inductancia de eje d (Henry)
        self.L_q = L_q        # Inductancia de eje q (Henry)
        self.psi_m = psi_m    # Flujo de imanes permanentes (Weber)
        self.J = J            # Inercia real del rotor + carga (kg*m^2)
        self.B = B            # Coeficiente de fricción viscosa real (N*m*s/rad)
        self.P = P            # Pares de polos (dientes del rotor)
        
        # Parámetros de cogging / detent
        self.Tdm = Tdm        # Amplitud del torque de detent (N*m)
        self.N_steps = N_steps# Pasos por revolución (frecuencia espacial de cogging)
        self.Phi = Phi        # Fase de cogging (rad)

        # Variables de estado iniciales
        self.i_d = 0.0        # Corriente de eje d (A)
        self.i_q = 0.0        # Corriente de eje q (A)
        self.omega_m = 0.0    # Velocidad mecánica angular (rad/s)
        self.theta_m = 0.0    # Posición mecánica angular (rad)

    @property
    def omega_e(self):
        """Velocidad eléctrica angular (rad/s)."""
        return self.P * self.omega_m

    @property
    def theta_e(self):
        """Posición eléctrica angular (rad) acotada en [0, 2*pi)."""
        return (self.P * self.theta_m) % (2.0 * np.pi)

    @property
    def torque_e(self):
        """
        Torque electromagnético generado por el motor stepper bifásico (N*m).
        Consistente con la S-Function PMSM_salient_dq.m de Simulink:
        Te = P * psi_m * i_q - P * (L_d - L_q) * i_d * i_q
        """
        return self.P * self.psi_m * self.i_q - self.P * (self.L_d - self.L_q) * self.i_d * self.i_q

    def get_states(self):
        """Retorna el vector de estados actual."""
        return np.array([self.i_d, self.i_q, self.omega_m, self.theta_m])

    def set_states(self, states):
        """Establece el vector de estados."""
        self.i_d, self.i_q, self.omega_m, self.theta_m = states

    def derivatives(self, states, v_d, v_q, T_L):
        """
        Calcula las derivadas dx/dt de los estados del motor.
        
        Estados:
        - states[0]: i_d
        - states[1]: i_q
        - states[2]: omega_m
        - states[3]: theta_m
        """
        i_d, i_q, omega_m, theta_m = states
        omega_e = self.P * omega_m

        # Ecuaciones eléctricas del estator en el marco d-q
        di_d = (v_d - self.R_s * i_d + omega_e * self.L_q * i_q) / self.L_d
        di_q = (v_q - self.R_s * i_q - omega_e * self.L_d * i_d - omega_e * self.psi_m) / self.L_q

        # Torque electromagnético
        T_e = self.P * self.psi_m * i_q - self.P * (self.L_d - self.L_q) * i_d * i_q
        
        # Torque de cogging / detent
        T_detent = self.Tdm * np.sin(self.N_steps * theta_m + self.Phi)

        # Dinámica mecánica
        # d_omega_m/dt = (Te - T_L - B_real*omega - T_detent) / J_real
        d_omega_m = (T_e - T_L - self.B * omega_m - T_detent) / self.J
        d_theta_m = omega_m

        return np.array([di_d, di_q, d_omega_m, d_theta_m])

    def step_rk4(self, v_d, v_q, T_L, dt):
        """
        Integra los estados del motor un paso 'dt' hacia adelante 
        utilizando Runge-Kutta de 4° orden (RK4).
        """
        x = self.get_states()

        k1 = self.derivatives(x, v_d, v_q, T_L)
        k2 = self.derivatives(x + 0.5 * dt * k1, v_d, v_q, T_L)
        k3 = self.derivatives(x + 0.5 * dt * k2, v_d, v_q, T_L)
        k4 = self.derivatives(x + dt * k3, v_d, v_q, T_L)

        x_new = x + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        
        # El ángulo mecánico theta_m se deja crecer sin límites para el seguimiento de posición acumulada.
        # La posición eléctrica theta_e se limita internamente a [0, 2*pi) mediante su @property.

        self.set_states(x_new)
        return self.get_states()

# Alias para compatibilidad con simuladores existentes
PMSMMotor = StepperMotor2Phase
