import numpy as np
import matplotlib.pyplot as plt
import parameters as par
import os
from simulation import PMSMSimulation

# Crear carpeta de salidas si no existe
OUTPUT_DIR = "outputs"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def get_position_reference(t):
    """
    Genera una rampa de posición mecánica (rad) para lograr
    una velocidad angular constante de 10 rad/s.
    theta(t) = 10 * t
    """
    return 10.0 * t

def get_speed_reference(t):
    """
    Genera la velocidad de referencia (rad/s) para prealimentación.
    """
    return 40.0

def get_load_torque(t):
    """
    Genera la perturbación de torque de carga externa (N*m).
    Escalón de 0.05 Nm en t = 0.5s para probar el rechazo de perturbaciones (DOB).
    """
    if t >= 0.5:
        return 0.05
    return 0.0

def export_to_interactive_scope(history, filename="scope_stepper.html"):
    """
    Genera un visor dinámico interactivo en formato HTML que funciona como
    un "Scope de Simulink". Permite hacer zoom, desplazar ejes, leer coordenadas
    al pasar el cursor y ocultar/mostrar señales dinámicamente haciendo clic en la leyenda.
    """
    filepath = os.path.join(OUTPUT_DIR, filename)
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        
        # Submuestreo para que el HTML no pese 300 MB y Python no colapse (Errno 22)
        step = 10 
        t = history['t'][::step]
        h = {k: v[::step] if isinstance(v, np.ndarray) else v for k, v in history.items()}
        
        # Crear subgráficos en cuadrícula de 6 filas x 2 columnas (12 gráficos totales)
        fig = make_subplots(
            rows=6, cols=2, 
            shared_xaxes=True,
            subplot_titles=(
                "Posición Mecánica (rad)", 
                "Velocidad Mecánica (rad/s)",
                "Corriente Eje d (Flujo) (A)", 
                "Corriente Eje q (Torque) (A)",
                "Corriente Fase A (A)",
                "Corriente Fase B (A)",
                "Ángulos Eléctricos (rad)",
                "Torques y Disturbios (N*m)",
                "Errores de Estimación (Vel y Ángulo)",
                "Amplitud HFI y Umbral",
                "Torque Cogging (Real vs LMS Est)"
            )
        )
        
        # Colores consistentes
        c_ref = '#e74c3c'  # Rojo
        c_real = '#2980b9' # Azul
        c_est = '#27ae60'  # Verde
        
        # 1. Posición (Fila 1, Col 1)
        fig.add_trace(go.Scatter(x=t, y=h['theta_ref'], name='Ref. Posición (θ_ref)', line=dict(color=c_ref, dash='dash')), row=1, col=1)
        fig.add_trace(go.Scatter(x=t, y=h['theta_m'], name='Real Posición (θ_m)', line=dict(color=c_real, width=2)), row=1, col=1)
        fig.add_trace(go.Scatter(x=t, y=h['theta_m_est'], name='Est. Posición EKF (θ_est)', line=dict(color=c_est, width=1.5, dash='dot')), row=1, col=1)
        
        # 2. Velocidad (Fila 1, Col 2)
        fig.add_trace(go.Scatter(x=t, y=h['omega_ref'], name='Ref. Velocidad (ω_ref)', line=dict(color=c_ref, dash='dash')), row=1, col=2)
        fig.add_trace(go.Scatter(x=t, y=h['omega_m'], name='Real Velocidad (ω_m)', line=dict(color=c_real, width=2)), row=1, col=2)
        fig.add_trace(go.Scatter(x=t, y=h['omega_m_est'], name='Est. Velocidad EKF (ω_est)', line=dict(color=c_est, width=1.5, dash='dot')), row=1, col=2)
        
        # 3. Corriente Id (Fila 2, Col 1)
        fig.add_trace(go.Scatter(x=t, y=h['i_d_ref'], name='Ref. Corriente Id', line=dict(color='#95a5a6', dash='dash')), row=2, col=1)
        fig.add_trace(go.Scatter(x=t, y=h['i_d'], name='Real Corriente Id', line=dict(color='#16a085')), row=2, col=1)
        fig.add_trace(go.Scatter(x=t, y=h['i_d_est'], name='Est. Corriente Id', line=dict(color='#f39c12', dash='dot')), row=2, col=1)

        # 4. Corriente Iq (Fila 2, Col 2)
        fig.add_trace(go.Scatter(x=t, y=h['i_q_ref'], name='Ref. Corriente Iq', line=dict(color='#e67e22', dash='dash')), row=2, col=2)
        fig.add_trace(go.Scatter(x=t, y=h['i_q'], name='Real Corriente Iq', line=dict(color='#9b59b6')), row=2, col=2)
        fig.add_trace(go.Scatter(x=t, y=h['i_q_est'], name='Est. Corriente Iq', line=dict(color=c_est, dash='dot')), row=2, col=2)
        
        # 5. Corriente Fase A (Fila 3, Col 1)
        fig.add_trace(go.Scatter(x=t, y=h['i_a'], name='i_a (Fase A)', line=dict(color='#1abc9c')), row=3, col=1)
        fig.add_trace(go.Scatter(x=t, y=h['i_a_est'], name='Est. i_a', line=dict(color=c_est, dash='dot')), row=3, col=1)
        
        # 6. Corriente Fase B (Fila 3, Col 2)
        fig.add_trace(go.Scatter(x=t, y=h['i_b'], name='i_b (Fase B)', line=dict(color='#e74c3c')), row=3, col=2)
        fig.add_trace(go.Scatter(x=t, y=h['i_b_est'], name='Est. i_b', line=dict(color='#f39c12', dash='dot')), row=3, col=2)

        # 7. Ángulos Eléctricos (Fila 4, Col 1)
        fig.add_trace(go.Scatter(x=t, y=h['theta_e'], name='Real θ_e', line=dict(color=c_real)), row=4, col=1)
        fig.add_trace(go.Scatter(x=t, y=h['theta_e_est'], name='Est. EKF θ_e', line=dict(color=c_est, dash='dot')), row=4, col=1)
        fig.add_trace(go.Scatter(x=t, y=h['theta_hfi_e'], name='Est. HFI θ_hfi', line=dict(color='#f39c12', dash='dash')), row=4, col=1)

        # 8. Torques (Fila 4, Col 2)
        fig.add_trace(go.Scatter(x=t, y=h['torque_e'], name='Torque Eléctrico (Te)', line=dict(color=c_real)), row=4, col=2)
        fig.add_trace(go.Scatter(x=t, y=h['torque_l'], name='Torque Carga (TL)', line=dict(color='#f39c12', dash='dash')), row=4, col=2)
        fig.add_trace(go.Scatter(x=t, y=h['torque_l_est'], name='Torque Carga Est.', line=dict(color='#8e44ad', dash='dot')), row=4, col=2)
        fig.add_trace(go.Scatter(x=t, y=h['d_hat_dob'], name='Disturbio Est.', line=dict(color=c_est, dash='dot')), row=4, col=2)
        
        # 9. Errores de Estimación (Fila 5, Col 1)
        fig.add_trace(go.Scatter(x=t, y=h['e_w'], name='Error Wm', line=dict(color='#c0392b')), row=5, col=1)
        fig.add_trace(go.Scatter(x=t, y=h['e_theta_ekf'], name='Error θ_ekf', line=dict(color=c_est)), row=5, col=1)
        fig.add_trace(go.Scatter(x=t, y=h['e_theta_hfi'], name='Error θ_hfi', line=dict(color='#f39c12', dash='dash')), row=5, col=1)
        
        # 10. Amplitud HFI (Fila 5, Col 2)
        fig.add_trace(go.Scatter(x=t, y=h['amp_hfi'], name='Amplitud HFI', line=dict(color='#8e44ad')), row=5, col=2)
        fig.add_trace(go.Scatter(x=t, y=h['threshold_hfi'], name='Umbral R_hfi', line=dict(color='gray', dash='dot')), row=5, col=2)

        # 11. Cogging (Fila 6, Col 1)
        fig.add_trace(go.Scatter(x=t, y=h['T_cogging_real'], name='T_cogging Real', line=dict(color=c_real)), row=6, col=1)
        fig.add_trace(go.Scatter(x=t, y=h['T_cogging_est'], name='T_cogging LMS Est', line=dict(color=c_est, dash='dot')), row=6, col=1)
        
        # Diseño general y configuración del eje temporal
        fig.update_layout(
            title_text="Visor Interactivo de Simulación (Scope) - Sensorless EKF+HFI+DOB",
            height=1250,
            showlegend=True,
            hovermode="x unified",
            template="plotly_white"
        )
        
        fig.update_xaxes(title_text="Tiempo (s)", row=5, col=1)
        fig.update_xaxes(title_text="Tiempo (s)", row=5, col=2)
        
        fig.write_html(filepath)
        print(f"\n[ÉXITO] Visor dinámico interconectado guardado como: '{filepath}'")
        print(f">> Puedes abrirlo desde la carpeta '{OUTPUT_DIR}' en tu navegador para hacer zoom.")
        
    except ImportError:
        print(f"\n[INFO] Instala plotly (pip install plotly) para guardar el visor interactivo en '{filepath}'.")

def plot_bpf_zoom(history, title_prefix=""):
    t = history['t']
    mask = (t >= 0.30) & (t <= 0.305)
    t_zoom = t[mask]
    
    if len(t_zoom) == 0:
        return
        
    fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle(f"Zoom BPF (0.30s - 0.305s) - {title_prefix}", fontsize=14, fontweight='bold')
    
    axs[0].plot(t_zoom, history['i_a'][mask], label='Ia_raw', color='blue', alpha=0.5, linewidth=2)
    axs[0].plot(t_zoom, history['i_a_fund'][mask], label='Ia_fund', color='green', linewidth=1.5)
    axs[0].plot(t_zoom, history['i_a_hfi'][mask], label='Ia_hfi', color='red', linewidth=1)
    axs[0].set_title("Separación BPF - Fase A")
    axs[0].set_ylabel("Corriente (A)")
    axs[0].grid(True)
    axs[0].legend()
    
    axs[1].plot(t_zoom, history['i_b'][mask], label='Ib_raw', color='blue', alpha=0.5, linewidth=2)
    axs[1].plot(t_zoom, history['i_b_fund'][mask], label='Ib_fund', color='green', linewidth=1.5)
    axs[1].plot(t_zoom, history['i_b_hfi'][mask], label='Ib_hfi', color='red', linewidth=1)
    axs[1].set_title("Separación BPF - Fase B")
    axs[1].set_xlabel("Tiempo (s)")
    axs[1].set_ylabel("Corriente (A)")
    axs[1].grid(True)
    axs[1].legend()
    
    plt.tight_layout()
    filename = f"sim_stepper_bpf_zoom_{title_prefix.lower().replace(' ', '_')}.png"
    filepath = os.path.join(OUTPUT_DIR, filename)
    plt.savefig(filepath, dpi=300)
    print(f"Gráfico BPF zoom guardado como: {filepath}")
    plt.close()

def plot_results(history, title_prefix=""):
    """Genera gráficos detallados y estéticos de la simulación en Matplotlib, incluyendo estimaciones."""
    t = history['t']
    
    # Crear figura con 12 subgráficos
    fig, axs = plt.subplots(12, 1, figsize=(10, 36), sharex=True)
    fig.suptitle(f"Simulación Stepper Bifásico - {title_prefix}", fontsize=14, fontweight='bold')
    
    # Colores elegantes
    c_ref = '#e74c3c'      # Rojo
    c_real = '#2980b9'     # Azul
    c_est = '#27ae60'      # Verde para estimación
    c_phase = ['#1abc9c', '#9b59b6'] # Verde azulado, Púrpura para fases A y B
    c_err = '#c0392b'      # Rojo oscuro para error
    c_hfi = '#f39c12'      # Naranja
    
    # 1. Gráfico de Posición
    axs[0].plot(t, history['theta_ref'], '--', color=c_ref, label='Referencia ($\theta_{ref}$)', linewidth=1.5)
    axs[0].plot(t, history['theta_m'], color=c_real, label='Real ($\theta_m$)', linewidth=2.0)
    axs[0].plot(t, history['theta_m_est'], ':', color=c_est, label='Estimada EKF ($\theta_{est}$)', linewidth=1.8)
    axs[0].set_ylabel('Posición (rad)', fontweight='bold')
    axs[0].grid(True, linestyle=':', alpha=0.6)
    axs[0].legend(loc='upper left')
    axs[0].set_title('Seguimiento de Posición Mecánica', fontsize=11, loc='left')
    
    # 2. Gráfico de Velocidad
    axs[1].plot(t, history['omega_ref'], '--', color=c_ref, label='Ref (Lazo Pos)', linewidth=1.5)
    axs[1].plot(t, history['omega_m'], color=c_real, label='Real ($\omega_m$)', linewidth=2.0)
    axs[1].plot(t, history['omega_m_est'], ':', color=c_est, label='Estimada EKF ($\omega_{est}$)', linewidth=1.8)
    axs[1].axhline(y=10.0, color='#7f8c8d', linestyle=':', label='Nominal (10 rad/s)', alpha=0.8)
    axs[1].set_ylabel('Velocidad (rad/s)', fontweight='bold')
    axs[1].grid(True, linestyle=':', alpha=0.6)
    axs[1].legend(loc='lower right')
    axs[1].set_title('Respuesta de Velocidad Mecánica', fontsize=11, loc='left')
    
    # 3. Gráfico de Corriente Id
    axs[2].plot(t, history['i_d_ref'], '--', color='#7f8c8d', label='$i_d$ Ref', linewidth=1.2)
    axs[2].plot(t, history['i_d'], color='#16a085', label='$i_d$ Real (Medida con EKF angle)', linewidth=1.5)
    axs[2].plot(t, history['i_d_est'], ':', color=c_hfi, label='$i_d$ Est (EKF)', linewidth=1.5)
    axs[2].set_ylabel('Corriente $i_d$ (A)', fontweight='bold')
    axs[2].grid(True, linestyle=':', alpha=0.6)
    axs[2].legend(loc='upper right')
    axs[2].set_title('Corriente eje d (Flujo)', fontsize=11, loc='left')
    
    # 4. Gráfico de Corriente Iq
    axs[3].plot(t, history['i_q_ref'], '--', color=c_ref, label='$i_q$ Ref', linewidth=1.2)
    axs[3].plot(t, history['i_q'], color=c_real, label='$i_q$ Real (Medida con EKF angle)', linewidth=1.5)
    axs[3].plot(t, history['i_q_est'], ':', color=c_est, label='$i_q$ Est (EKF)', linewidth=1.5)
    axs[3].set_ylabel('Corriente $i_q$ (A)', fontweight='bold')
    axs[3].grid(True, linestyle=':', alpha=0.6)
    axs[3].legend(loc='upper right')
    axs[3].set_title('Corriente eje q (Torque)', fontsize=11, loc='left')
    
    # 5. Gráfico de Corriente Ia
    axs[4].plot(t, history['i_a'], color=c_phase[0], label='$i_a$ (Fase A)', linewidth=1.2)
    axs[4].plot(t, history['i_a_est'], ':', color=c_est, label='$i_a$ Est', linewidth=1.0)
    axs[4].set_ylabel('Corriente $i_a$ (A)', fontweight='bold')
    axs[4].grid(True, linestyle=':', alpha=0.6)
    axs[4].legend(loc='upper right')
    axs[4].set_title('Corriente Estatórica Fase A', fontsize=11, loc='left')
    
    # 6. Gráfico de Corriente Ib
    axs[5].plot(t, history['i_b'], color=c_phase[1], label='$i_b$ (Fase B)', linewidth=1.2)
    axs[5].plot(t, history['i_b_est'], ':', color=c_hfi, label='$i_b$ Est', linewidth=1.0)
    axs[5].set_ylabel('Corriente $i_b$ (A)', fontweight='bold')
    axs[5].grid(True, linestyle=':', alpha=0.6)
    axs[5].legend(loc='upper right')
    axs[5].set_title('Corriente Estatórica Fase B', fontsize=11, loc='left')
    
    # 7. Gráfico de Ángulos Eléctricos
    axs[6].plot(t, history['theta_e'], color=c_real, label='Real ($\theta_e$)', linewidth=2.0)
    axs[6].plot(t, history['theta_e_est'], ':', color=c_est, label='Estimado EKF ($\theta_{e\_est}$)', linewidth=1.8)
    axs[6].plot(t, history['theta_hfi_e'], '--', color=c_hfi, label='Estimado HFI ($\theta_{e\_hfi}$)', linewidth=1.8)
    axs[6].set_ylabel('Ángulo (rad)', fontweight='bold')
    axs[6].grid(True, linestyle=':', alpha=0.6)
    axs[6].legend(loc='lower right')
    axs[6].set_title('Seguimiento de Ángulo Eléctrico (HFI + EKF)', fontsize=11, loc='left')
    
    # 8. Gráfico de Error de Velocidad y Errores Angulares
    axs[7].plot(t, history['e_w'], color=c_err, label='$e_W$ (Error Vel)', linewidth=1.5)
    axs[7].set_ylabel('Error Wm (rad/s)', fontweight='bold')
    axs[7].grid(True, linestyle=':', alpha=0.6)
    axs[7].legend(loc='upper left')
    axs[7].set_title('Error de Estimación de Velocidad', fontsize=11, loc='left')
    
    # 9. Gráfico de Error de Ángulo Eléctrico
    axs[8].plot(t, history['e_theta_ekf'], color=c_est, label='$e_{\\theta\_ekf}$', linewidth=1.5)
    axs[8].plot(t, history['e_theta_hfi'], '--', color=c_hfi, label='$e_{\\theta\_hfi}$', linewidth=1.5)
    axs[8].set_ylabel('Error Ángulo (rad)', fontweight='bold')
    axs[8].grid(True, linestyle=':', alpha=0.6)
    axs[8].legend(loc='upper left')
    axs[8].set_title('Error de Estimación Angular Eléctrica', fontsize=11, loc='left')
    
    # 10. Gráfico de Amplitud HFI
    axs[9].plot(t, history['amp_hfi'], color=c_phase[1], label='Amplitud HFI', linewidth=1.5)
    axs[9].plot(t, history['threshold_hfi'], color='gray', linestyle=':', label='Umbral R_hfi')
    axs[9].set_ylabel('Amplitud (A)', fontweight='bold')
    axs[9].grid(True, linestyle=':', alpha=0.6)
    axs[9].legend(loc='upper left')
    axs[9].set_title('Amplitud Detectada HFI (Idh)', fontsize=11, loc='left')
    
    # 11. Gráfico de Torque de Carga y Perturbación
    axs[10].plot(t, history['torque_l'], color=c_real, label='$T_x$ Real (Carga)', linewidth=2.0)
    axs[10].plot(t, history['d_hat_dob'], '--', color=c_est, label='$T_x$ Est (DOB)', linewidth=1.5)
    axs[10].set_ylabel('Torque (Nm)', fontweight='bold')
    axs[10].set_xlabel('Tiempo (s)', fontweight='bold')
    axs[10].grid(True, linestyle=':', alpha=0.6)
    axs[10].legend(loc='upper left')
    axs[10].set_title('Estimación de Perturbación (Tx)', fontsize=11, loc='left')
    
    # 12. Gráfico de Cogging Real vs Estimado
    axs[11].plot(t, history['T_cogging_real'], color=c_real, label='$T_{cogging}$ Real', linewidth=2.0)
    axs[11].plot(t, history['T_cogging_est'], '--', color=c_est, label='$T_{cogging}$ Est (LMS)', linewidth=1.5)
    axs[11].set_ylabel('Cogging (Nm)', fontweight='bold')
    axs[11].set_xlabel('Tiempo (s)', fontweight='bold')
    axs[11].grid(True, linestyle=':', alpha=0.6)
    axs[11].legend(loc='upper left')
    axs[11].set_title('Estimación Adaptativa de Cogging (LMS)', fontsize=11, loc='left')
    
    plt.tight_layout()
    
    filename = f"sim_stepper_{title_prefix.lower().replace(' ', '_')}.png"
    filepath = os.path.join(OUTPUT_DIR, filename)
    plt.savefig(filepath, dpi=300)
    print(f"Gráfico guardado como: {filepath}")
    plt.close()


if __name__ == "__main__":
    print("======================================================================")
    print(" SIMULACIÓN DE MOTOR STEPPER BIFÁSICO EN TIEMPO MIXTO (CASCADA 3 LAZOS) ")
    print("======================================================================")
    print(f"Parámetros temporales configurados (parameters.py):")
    print(f" - Simulación planta continua (dt_sim): {par.dt_sim*1e6:.1f} us")
    print(f" - Lazo discreto corriente (dt_current): {par.dt_current*1e6:.1f} us ({1/par.dt_current*1e-3:.1f} kHz)")
    print(f" - Lazo discreto velocidad (dt_speed): {par.dt_speed*1e6:.1f} us ({1/par.dt_speed*1e-3:.1f} kHz)")
    print(f" - Lazo discreto posición (dt_pos): {par.dt_pos*1e6:.1f} us ({1/par.dt_pos:.1f} Hz)")
    print(f" - Lazo del observador (dt_observer): {par.dt_observer*1e6:.1f} us ({1/par.dt_observer*1e-3:.1f} kHz)")
    print("----------------------------------------------------------------------")
    
    # =========================================================================
    # PANEL DE ENRUTAMIENTO DE SEÑALES (TIPO GOTO SIMULINK)
    # =========================================================================
    # Modifica el valor de cada entrada para cambiar la señal que recibe cada bloque.
    # Señales disponibles: 
    #  - Posición: 'theta_m_real', 'theta_m_ekf'
    #  - Velocidad: 'Wm_real', 'Wm_ekf'
    #  - Ángulo Eléctrico: 'theta_e_real', 'theta_e_ekf'
    # =========================================================================
    signals = {
        # 1. Lazo de Posición
        'pos_theta_m': 'theta_m_real',
        'pos_Wm': 'Wm_real',
        
        # 2. Lazo de Velocidad
        'speed_Wm_ref': 'omega_ref_ext', # <--- Se ignora el lazo de posición, referencia fija
        'speed_Wm': 'Wm_real',
        
        # 3. Lazo de Corriente
        'curr_Wm': 'Wm_real',
        'curr_Theta_e': 'theta_e_real'
    }
    
    print(f"\nConexión de Señales en los Lazos:")
    print(f" - Lazo de Posición:  {signals['pos_theta_m']} / {signals['pos_Wm']}")
    print(f" - Lazo de Velocidad: {signals['speed_Wm']}")
    print(f" - Lazo de Corriente: {signals['curr_Theta_e']} / {signals['curr_Wm']}")
    
    # Determinar prefijo de título para los archivos
    p_fb_type = "ekf" if "ekf" in signals['pos_theta_m'] else "real"
    s_fb_type = "ekf" if "ekf" in signals['speed_Wm'] else "real"
    c_fb_type = "ekf" if "ekf" in signals['curr_Theta_e'] else "real"
    
    if p_fb_type == 'ekf' and s_fb_type == 'ekf' and c_fb_type == 'ekf':
        title_prefix = "Sensorless"
    elif p_fb_type == 'real' and s_fb_type == 'real' and c_fb_type == 'real':
        title_prefix = "Sensado"
    else:
        title_prefix = f"Hibrido_P_{p_fb_type}_S_{s_fb_type}_C_{c_fb_type}"
        
    print(f"\nEjecutando Simulación: Lazo Cerrado en modo {title_prefix}...")
    sim = PMSMSimulation(
        inverter_mode='average',
        signals_routing=signals
    )
    
    history = sim.run(
        get_theta_ref=get_position_reference,
        get_load_torque=get_load_torque,
        get_omega_ref=get_speed_reference
    )
    
    # 1. Gráficos de Matplotlib
    plot_bpf_zoom(history, title_prefix=title_prefix)
    plot_results(history, title_prefix=title_prefix)
    
    # 2. Visor Interactivo (Scope)
    export_to_interactive_scope(history, filename="scope_stepper.html")
    
    print("\n¡Simulación completada con éxito!")
    print("======================================================================")
