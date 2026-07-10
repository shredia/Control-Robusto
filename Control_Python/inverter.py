import numpy as np

class TwoPhaseInverter:
    """
    Modelo de inversor bifásico (compuesto por 2 puentes H completos independientes).
    Soporta:
    1. Modelo Promedio (Average Model): Aplica los voltajes de referencia directamente,
       limitándolos al voltaje del bus DC disponible ([-Vdc, Vdc]).
    2. Modelo de Conmutación (Switching Model): Simula la conmutación de los interruptores
       de los puentes H comparando las referencias con una portadora triangular de frecuencia f_sw.
    """
    def __init__(self, V_dc=48.0, f_sw=10000.0):
        self.V_dc = V_dc       # Voltaje del bus DC (V)
        self.f_sw = f_sw       # Frecuencia de conmutación PWM (Hz)

    def calculate_duty_cycles(self, v_a_ref, v_b_ref):
        """
        Calcula los Ciclos de Trabajo (Duty Cycles) bipolarizados entre 0 y 1.
        Para un puente H completo, un ciclo de trabajo d = 0.5 produce un voltaje medio de 0 V,
        d = 1.0 produce Vdc y d = 0.0 produce -Vdc.
        """
        # Limitar las referencias al bus de alimentación
        v_a_sat = np.clip(v_a_ref, -self.V_dc, self.V_dc)
        v_b_sat = np.clip(v_b_ref, -self.V_dc, self.V_dc)
        
        # Mapear de [-Vdc, Vdc] a [0, 1]
        d_a = (v_a_sat / (2.0 * self.V_dc)) + 0.5
        d_b = (v_b_sat / (2.0 * self.V_dc)) + 0.5
        
        return v_a_sat, v_b_sat, d_a, d_b

    def get_voltage_switching(self, d_a, d_b, t):
        """
        Simula la conmutación bipolar real comparando los duty cycles con una portadora triangular.
        Retorna los voltajes de fase instantáneos v_a, v_b aplicados al estator.
        """
        # Periodo de conmutación
        T_sw = 1.0 / self.f_sw
        
        # Tiempo normalizado dentro del periodo portador [0, T_sw]
        t_norm = t % T_sw
        
        # Portadora triangular simétrica entre 0 y 1
        if t_norm < (T_sw / 2.0):
            carrier = t_norm / (T_sw / 2.0)
        else:
            carrier = (T_sw - t_norm) / (T_sw / 2.0)
            
        # Conmutación del puente H (Bipolar): 1 si el ciclo de trabajo supera a la portadora, -1 si no.
        S_a = 1.0 if d_a > carrier else -1.0
        S_b = 1.0 if d_b > carrier else -1.0
        
        # Voltajes aplicados a las fases A y B
        v_a = S_a * self.V_dc
        v_b = S_b * self.V_dc
        
        return v_a, v_b

    def step(self, v_a_ref, v_b_ref, t, mode='average'):
        """
        Ejecuta el inversor bifásico y retorna los voltajes de fase v_a, v_b aplicados al motor.
        
        Modos:
        - 'average': Modulación promedio (voltajes limitados directamente a [-Vdc, Vdc]).
        - 'switching': Simulación por conmutación física comparando con portadora a f_sw.
        """
        v_a_sat, v_b_sat, d_a, d_b = self.calculate_duty_cycles(v_a_ref, v_b_ref)
        
        if mode == 'average':
            return v_a_sat, v_b_sat
        elif mode == 'switching':
            return self.get_voltage_switching(d_a, d_b, t)
        else:
            raise ValueError(f"Modo de inversor bifásico no válido: {mode}")
