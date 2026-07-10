function PMSM_salient(block)
    setup(block);
end

function setup(block)

    % Entradas: Va, Vb, Tl_ext
    block.NumInputPorts  = 3;

    % Salidas:
    % Ia, Ib, Wm, We, Th_m, Th_e, Te, Tdetent, Tnet
    block.NumOutputPorts = 10;

    block.SetPreCompInpPortInfoToDynamic;
    block.SetPreCompOutPortInfoToDynamic;

    % Va y Vb no afectan directamente las salidas.
    % Tl_ext sí afecta directamente Tnet.
    block.InputPort(1).DirectFeedthrough = false;
    block.InputPort(2).DirectFeedthrough = false;
    block.InputPort(3).DirectFeedthrough = true;

    block.NumDialogPrms = 1;
    block.NumContStates = 4;
    block.SampleTimes   = [0 0];

    block.RegBlockMethod('InitializeConditions', @InitConditions);
    block.RegBlockMethod('Outputs',              @Outputs);
    block.RegBlockMethod('Derivatives',          @Derivatives);

end

function InitConditions(block)

    % Estados: Ia, Ib, Wm, Th_m
    block.ContStates.Data = [0; 0; 0; 0];

end

function [Te, L_mat, dL_dTh, Tdetent, Tload_pos, Tnet,Tneg] = ...
    calculate_physics(p, Ia, Ib, Wm, Th_m, Tl_ext)

    Th_e   = p.P * Th_m;
    Th_e_2 = 2 * Th_e;

    L0 = (p.Ld + p.Lq)/2;
    L2 = (p.Ld - p.Lq)/2;

    % Matriz de inductancias en coordenadas a-b
    L_mat = [L0 + L2*cos(Th_e_2),  L2*sin(Th_e_2);
             L2*sin(Th_e_2),      L0 - L2*cos(Th_e_2)];

    % Derivada respecto del ángulo eléctrico
    dL_dTh = [-2*L2*sin(Th_e_2),   2*L2*cos(Th_e_2);
                2*L2*cos(Th_e_2), 2*L2*sin(Th_e_2)];

    I_vec = [Ia; Ib];

    % Torque de reluctancia
    Te_reluc = 0.5 * I_vec' * (p.P * dL_dTh) * I_vec;

    % Torque síncrono de imanes
    Te_sync = p.Ke * (Ib*cos(Th_e) - Ia*sin(Th_e));

    % Torque electromagnético total
    Te = Te_sync + Te_reluc;

    % Torque interno de detent/cogging
    Tdetent = p.Tdm * sin(p.Nr*Th_m + p.Phi);

    % Carga dependiente de posición mecánica.
    % Ejemplo: gravedad, excentricidad, poleas, desbalance, etc.
    Tload_pos = p.Tload_amp * ...
                sin(p.Nload*Th_m + p.Phi_load);

    % Suma firmada de todos los torques mecánicos

    Tneg = + Tl_ext ...
         + Tload_pos ...
         + p.B_real*Wm ...
         + Tdetent;
     Tnet = Te ...
         - Tneg;
end

function Outputs(block)

    p = block.DialogPrm(1).Data;

    x    = block.ContStates.Data;
    Ia   = x(1);
    Ib   = x(2);
    Wm   = x(3);
    Th_m = x(4);

    Tl_ext = block.InputPort(3).Data;

    [Te, ~, ~, Tdetent, ~, Tnet,Tneg] = ...
        calculate_physics(p, Ia, Ib, Wm, Th_m, Tl_ext);

    block.OutputPort(1).Data = Ia;
    block.OutputPort(2).Data = Ib;
    block.OutputPort(3).Data = Wm;
    block.OutputPort(4).Data = p.P * Wm;
    block.OutputPort(5).Data = Th_m;
    block.OutputPort(6).Data = p.P * Th_m;
    block.OutputPort(7).Data = Te;
    block.OutputPort(8).Data = Tdetent;
    block.OutputPort(9).Data = Tnet;
    block.OutputPort(10).Data = Tneg;

end

function Derivatives(block)

    p = block.DialogPrm(1).Data;

    Va     = block.InputPort(1).Data;
    Vb     = block.InputPort(2).Data;
    Tl_ext = block.InputPort(3).Data;

    x    = block.ContStates.Data;
    Ia   = x(1);
    Ib   = x(2);
    Wm   = x(3);
    Th_m = x(4);

    I_vec = [Ia; Ib];
    Th_e  = p.P * Th_m;

    [~, L_mat, dL_dTh, ~, ~, Tnet, ~] = ...
        calculate_physics(p, Ia, Ib, Wm, Th_m, Tl_ext);

    % Back-EMF
    E_vec = [-p.Ke * Wm * sin(Th_e);
              p.Ke * Wm * cos(Th_e)];

    % dL/dt = dL/dtheta_e * omega_e
    dL_dt = dL_dTh * (p.P * Wm);

    % Dinámica eléctrica
    rhs_elec = [Va; Vb] ...
             - p.R * I_vec ...
             - dL_dt * I_vec ...
             - E_vec;

    dI_vec = L_mat \ rhs_elec;

    % Dinámica mecánica
    dWm = Tnet / p.J_real;
    dTh = Wm;

    block.Derivatives.Data = [dI_vec(1); dI_vec(2); dWm; dTh];

end