import os
import numpy as np
import matplotlib.pyplot as plt
import parameters as par
from simulation import PMSMSimulation

OUTPUT_DIR = "outputs"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def get_position_reference(t):
    return 10.0 * t

def get_speed_reference(t):
    if t < 0.1:
        return 0.0
    elif t < 0.3:
        return 50.0 * (t - 0.1)
    else:
        return 10.0

def get_load_torque(t, theta_m=0.0):
    # Escalón de torque de carga
    T_L = 0.02 if t >= 0.5 else 0.0
    # Torque gravitacional senoidal
    T_grav = 0.01 * np.sin(theta_m)
    # Roce o desbalance variable lento
    T_fric = 0.005 * np.sin(2.0 * np.pi * 0.5 * t)
    return T_L + T_grav + T_fric

def run_case(case_name, use_bpf, use_hfi, qk_tx, title):
    print(f"\n--- Ejecutando {case_name}: {title} ---")
    signals = {
        'pos_theta_m': 'theta_m_real',
        'pos_Wm': 'Wm_real',
        'speed_Wm_ref': 'omega_ref_ext',
        'speed_Wm': 'Wm_real',
        'curr_Wm': 'Wm_real',
        'curr_Theta_e': 'theta_e_real'
    }
    
    sim = PMSMSimulation(
        inverter_mode='average',
        signals_routing=signals
    )
    
    # Configurar observador
    sim.ekf_hfi.use_bpf = use_bpf
    sim.ekf_hfi.Qk_override = np.diag([1e-4, 1e-4, 1e-7, 1e-10, qk_tx])
    
    if not use_hfi:
        sim.ekf_hfi.force_hfi_off = True
        sim.current_controller.HFI_enable = 0
    else:
        sim.ekf_hfi.force_hfi_off = False
        sim.current_controller.HFI_enable = 1
        
    history = sim.run(
        get_theta_ref=get_position_reference,
        get_load_torque=get_load_torque,
        get_omega_ref=get_speed_reference
    )
    
    return history

def print_table(histories):
    print("\n" + "="*85)
    print(" TABLA COMPARATIVA DE CASOS EKF + HFI + DOB ")
    print("="*85)
    print(f"{'Caso':<8} | {'RMSE Vel (rad/s)':<16} | {'RMSE Ángulo (°_e)':<17} | {'RMSE Tx (Nm)':<13} | {'HFI Activación (%)':<18}")
    print("-" * 85)
    
    names = ["Caso A", "Caso B", "Caso C", "Caso D"]
    for i, h in enumerate(histories):
        print(f"{names[i]:<8} | {h['rmse_w']:<16.4f} | {h['rmse_theta_deg']:<17.2f} | {h['rmse_tx']:<13.4f} | {h['hfi_activation_ratio']*100:<18.1f}")
    print("="*85 + "\n")
    
    print("Observaciones de estabilidad:")
    print(" - Caso A (Base): Tx asume 0, no compensa perturbaciones, el error angular suele ser elevado si no hay BPF porque el ruido afecta el estado.")
    print(" - Caso B (+BPF): Al filtrar componentes de alta frecuencia, mejora la limpieza del modelo fundamental, pero sin DOB no sigue perturbaciones.")
    print(" - Caso C (+DOB): El observador logra seguir la perturbación lenta (gravedad, roce y escalón), mejorando radicalmente la velocidad.")
    print(" - Caso D (+HFI): Añade corrección sensoless, activándose solo cuando SNR es seguro, acotando el error acumulado.")

def plot_comparisons(histories):
    cA, cB, cC, cD = histories
    t = cD['t']
    
    fig, axs = plt.subplots(4, 2, figsize=(16, 14), sharex=True)
    fig.suptitle('Comparación EKF: Progresión de Arquitecturas', fontsize=16, fontweight='bold')
    
    # 1. Velocidad Real vs Est
    axs[0,0].plot(t, cD['omega_m'], 'k', label='Real')
    axs[0,0].plot(t, cB['omega_m_est'], 'g--', label='Est Caso B')
    axs[0,0].plot(t, cC['omega_m_est'], 'b--', label='Est Caso C (DOB)')
    axs[0,0].plot(t, cD['omega_m_est'], 'm:', label='Est Caso D (HFI)')
    axs[0,0].set_title('Velocidad Mecánica (rad/s)')
    axs[0,0].legend()
    axs[0,0].grid()
    
    # 2. Error Velocidad
    axs[0,1].plot(t, cA['e_w'], 'r', alpha=0.6, label=f"Caso A (RMSE={cA['rmse_w']:.2f})")
    axs[0,1].plot(t, cB['e_w'], 'g', alpha=0.6, label=f"Caso B (RMSE={cB['rmse_w']:.2f})")
    axs[0,1].plot(t, cC['e_w'], 'b', label=f"Caso C (RMSE={cC['rmse_w']:.2f})")
    axs[0,1].plot(t, cD['e_w'], 'm--', label=f"Caso D (RMSE={cD['rmse_w']:.2f})")
    axs[0,1].set_title('Error de Velocidad e_w (rad/s)')
    axs[0,1].legend()
    axs[0,1].grid()
    
    # 3. Ángulos Eléctricos
    axs[1,0].plot(t, cD['theta_e'], 'k', label='Real')
    axs[1,0].plot(t, cD['theta_e_est'], 'b--', label='Est EKF (Caso D)')
    axs[1,0].plot(t, cD['theta_hfi_e'], 'm:', label='Est HFI (Caso D)')
    axs[1,0].set_title('Ángulo Eléctrico (rad)')
    axs[1,0].legend()
    axs[1,0].grid()
    
    # 4. Error Angular Circular
    axs[1,1].plot(t, cA['e_theta_ekf'], 'r', alpha=0.6, label=f"Caso A (RMSE={cA['rmse_theta']:.2f} rad)")
    axs[1,1].plot(t, cC['e_theta_ekf'], 'b', label=f"Caso C (RMSE={cC['rmse_theta']:.2f} rad)")
    axs[1,1].plot(t, cD['e_theta_ekf'], 'm--', label=f"Caso D (RMSE={cD['rmse_theta']:.2f} rad)")
    axs[1,1].set_title('Error Angular Eléctrico circular (rad)')
    axs[1,1].legend()
    axs[1,1].grid()
    
    # 5. Torque Perturbador Tx
    axs[2,0].plot(t, cD['torque_l'], 'k', label='Tx Real (Perturbaciones)')
    axs[2,0].plot(t, cA['d_hat_dob'], 'r--', label=f"Est Caso A (RMSE={cA['rmse_tx']:.2f})")
    axs[2,0].plot(t, cC['d_hat_dob'], 'b', label=f"Est Caso C (RMSE={cC['rmse_tx']:.2f})")
    axs[2,0].plot(t, cD['d_hat_dob'], 'm--', label=f"Est Caso D (RMSE={cD['rmse_tx']:.2f})")
    axs[2,0].set_title('Estimación de Torque Perturbador Tx (Nm)')
    axs[2,0].legend()
    axs[2,0].grid()
    
    # 6. Ganancias Kalman Kw, Ktx (Caso D)
    axs[2,1].plot(t, cD['K_w'], 'c', label='Kw (Caso D)')
    axs[2,1].plot(t, cD['K_tx'], 'orange', label='Ktx (Caso D)')
    axs[2,1].set_title('Ganancias de Kalman (Caso D)')
    axs[2,1].set_yscale('log')
    axs[2,1].legend()
    axs[2,1].grid()
    
    # 7. Amplitud HFI y Activación
    axs[3,0].plot(t, cD['amp_hfi'], 'purple', label='Amp HFI Caso D')
    axs[3,0].plot(t, cD['threshold_hfi'], 'k:', label='Threshold')
    axs[3,0].set_title('Amplitud HFI y Umbral Adaptativo (A)')
    axs[3,0].set_xlabel('Tiempo (s)')
    axs[3,0].legend()
    axs[3,0].grid()
    
    # 8. HFI Valid
    axs[3,1].fill_between(t, 0, cD['hfi_valid'], color='purple', alpha=0.3, label='HFI Active (Caso D)')
    axs[3,1].set_title('Estado de Validación HFI (hfi_valid)')
    axs[3,1].set_xlabel('Tiempo (s)')
    axs[3,1].set_yticks([0, 1])
    axs[3,1].legend()
    axs[3,1].grid()
    
    plt.tight_layout()
    filepath = os.path.join(OUTPUT_DIR, "comparacion_ekf_casos.png")
    plt.savefig(filepath, dpi=300)
    print(f"Gráfico comparativo guardado en {filepath}")
    plt.close()

if __name__ == "__main__":
    print("=========================================================")
    print(" EXPERIMENTOS DE ARQUITECTURA EKF + HFI + DOB ")
    print("=========================================================")
    
    # Caso A: EKF base, sin HFI y sin DOB.
    hA = run_case("Caso A", use_bpf=False, use_hfi=False, qk_tx=1e-12, title="Base (Sin BPF, Sin DOB, Sin HFI)")
    
    # Caso B: EKF + separación BPF, sin actualización HFI y con Tx rígido.
    hB = run_case("Caso B", use_bpf=True, use_hfi=False, qk_tx=1e-12, title="+ BPF (Sin DOB, Sin HFI)")
    
    # Caso C: EKF + separación BPF + DOB.
    hC = run_case("Caso C", use_bpf=True, use_hfi=False, qk_tx=1e-8, title="+ BPF + DOB (Sin HFI)")
    
    # Caso D: EKF + separación BPF + DOB + HFI validado.
    hD = run_case("Caso D", use_bpf=True, use_hfi=True, qk_tx=1e-8, title="+ BPF + DOB + HFI Adaptativo")
    
    print_table([hA, hB, hC, hD])
    plot_comparisons([hA, hB, hC, hD])
    
    print("Experimentos finalizados.")
