import numpy as np
import os

class RealTimeLMSCogging:
    """
    Filtro Adaptativo LMS en tiempo real para aprendizaje y feedforward de cogging.
    Estima los primeros 2 armónicos espaciales del torque de cogging.
    Diseñado para ser portado a C en un microcontrolador.
    """
    def __init__(self, N_steps=200, fbw_DO=25.0, learning_rate=0.01, turns_to_learn=2.0):
        self.N_steps = N_steps
        self.fbw_DO = fbw_DO
        self.Td = 1.0 / (2.0 * np.pi * fbw_DO) if fbw_DO > 0 else 0.0
        self.mu = learning_rate
        self.turns_to_learn = turns_to_learn
        
        # Pesos del filtro LMS: [w1_sin1, w2_cos1, w3_sin2, w4_cos2]
        self.W = np.zeros(4)
        
        self.accumulated_theta = 0.0
        self.last_theta_m = None
        self.mode = 'learn'
        # Filtro Pasa-Altas (HPF) para aislar el rizado de d_hat
        self.d_hat_lpf = 0.0
        self.alpha_lpf = 0.02  # Constante de tiempo para LPF (~2% de actualización)

    def step(self, theta_m, omega_m, d_hat):
        """
        Debe ser llamado en cada paso de control de velocidad.
        Retorna el torque de feedforward de cogging a inyectar estimado por el LMS.
        """
        if self.last_theta_m is not None:
            delta_theta = theta_m - self.last_theta_m
            # Manejar wrap-around natural si el ángulo viene del EKF (-pi a pi)
            if delta_theta > np.pi:
                delta_theta -= 2.0 * np.pi
            elif delta_theta < -np.pi:
                delta_theta += 2.0 * np.pi
            self.accumulated_theta += abs(delta_theta)
            
        self.last_theta_m = theta_m
        turns = self.accumulated_theta / (2.0 * np.pi)
        
        # Compensar el retraso del filtro DOB para alinear en fase
        theta_true = theta_m - omega_m * self.Td
        
        # Vector de excitación X para 2 armónicos
        # Armónico fundamental: N_steps * theta
        # Segundo armónico: 2 * N_steps * theta
        arg1 = self.N_steps * theta_true
        arg2 = 2.0 * self.N_steps * theta_true
        
        X = np.array([
            np.sin(arg1),
            np.cos(arg1),
            np.sin(arg2),
            np.cos(arg2)
        ])
        
        # Estimación actual del torque de cogging
        y_hat = np.dot(self.W, X)
        
        # 1. Aislar la componente alterna (AC / Rizado) del Disturbio
        # Usamos un filtro Pasa-Bajos simple (EWMA) para encontrar el DC (escalones)
        self.d_hat_lpf = self.alpha_lpf * d_hat + (1.0 - self.alpha_lpf) * self.d_hat_lpf
        # La componente AC (rizado) es la diferencia
        d_hat_hpf = d_hat - self.d_hat_lpf
        
        # 2. ¡EL SECRETO DEL FEEDFORWARD ADAPTATIVO EN LAZO CERRADO!
        # d_hat (Tx del EKF) ya es el *residuo* o perturbación no modelada.
        # Si el LMS cancela el cogging, Tx no tendrá rizado (d_hat_hpf = 0).
        # Por lo tanto, el error para el LMS *es* directamente d_hat_hpf.
        error = d_hat_hpf
        
        # 3. Actualización de pesos LMS: W(k+1) = W(k) + mu * error * X
        self.W += self.mu * error * X
        
        # Aplicación: Extraer el torque estimado por el modelo de Fourier (LMS)
        # Usar theta_m actual (no retrasado) para el feedforward inyectado
        arg1_out = self.N_steps * theta_m
        arg2_out = 2.0 * self.N_steps * theta_m
        
        X_out = np.array([
            np.sin(arg1_out),
            np.cos(arg1_out),
            np.sin(arg2_out),
            np.cos(arg2_out)
        ])
        
        self.mode = 'learn_and_apply'
        return np.dot(self.W, X_out)
