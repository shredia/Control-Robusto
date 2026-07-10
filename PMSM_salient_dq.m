function PMSM_salient_dq(block)
    setup(block);
end

function setup(block)

    % Entradas: Vd, Vq, Tl_ext
    block.NumInputPorts  = 3;

    % Salidas:
    % Id, Iq, Wm, We, Th_m, Th_e, Te, Tdet, Tnet
    block.NumOutputPorts = 10;

    block.SetPreCompInpPortInfoToDynamic;
    block.SetPreCompOutPortInfoToDynamic;

    for k = 1:block.NumInputPorts
        block.InputPort(k).Dimensions = 1;
    end

    % Vd y Vq no afectan directamente las salidas instantáneas.
    % Tl_ext sí afecta Tnet.
    block.InputPort(1).DirectFeedthrough = false;
    block.InputPort(2).DirectFeedthrough = false;
    block.InputPort(3).DirectFeedthrough = true;

    for k = 1:block.NumOutputPorts
        block.OutputPort(k).Dimensions = 1;
    end

    block.NumDialogPrms = 1;
    block.NumContStates = 4;   % [Id Iq Wm Th_m]
    block.SampleTimes   = [0 0];

    block.RegBlockMethod('InitializeConditions', @InitConditions);
    block.RegBlockMethod('Outputs',              @Outputs);
    block.RegBlockMethod('Derivatives',          @Derivatives);

end

function InitConditions(block)

    % Estados: [Id; Iq; Wm; Th_m]
    block.ContStates.Data = [0; 0; 0; 0];

end

function [Te, Tdet, Tload_pos, Tvisc, Tnet,Tneg] = ...
    calculate_mechanics(p, Id, Iq, Wm, Th_m, Tl_ext)

    % Torque electromagnético:
    % Te_PM  = Ke*Iq
    % Te_rel = P*(Lq-Ld)*Id*Iq
    Te = p.Ke*Iq - p.P*(p.Ld - p.Lq)*Id*Iq;

    % Torque de detent/cogging interno del motor
    Tdet = p.Tdm * sin(p.Nr*Th_m + p.Phi);

    % Carga dependiente de la posición mecánica.
    % Puede representar gravedad, desbalance, poleas,
    % transmisión, excentricidad, etc.
    Tload_pos = p.Tload_amp * ...
                sin(p.Nload*Th_m + p.Phi_load);

    % Fricción viscosa
    Tvisc = p.B_real * Wm;

    % Suma algebraica de torques sobre el eje.
    % Este torque determina directamente la aceleración.

    Tneg =  + Tl_ext ...
         + Tload_pos ...
         + Tvisc ...
         + Tdet;

     Tnet =  Te - Tneg;

end

function Outputs(block)

    p = block.DialogPrm(1).Data;
    x = block.ContStates.Data;

    Id   = x(1);
    Iq   = x(2);
    Wm   = x(3);
    Th_m = x(4);

    Th_e = p.P * Th_m;
    We   = p.P * Wm;

    Tl_ext = block.InputPort(3).Data;

    [Te, Tdet, ~, ~, Tnet,Tneg] = ...
        calculate_mechanics(p, Id, Iq, Wm, Th_m, Tl_ext);

    block.OutputPort(1).Data = Id;
    block.OutputPort(2).Data = Iq;
    block.OutputPort(3).Data = Wm;
    block.OutputPort(4).Data = We;
    block.OutputPort(5).Data = Th_m;
    block.OutputPort(6).Data = Th_e;
    block.OutputPort(7).Data = Te;
    block.OutputPort(8).Data = Tdet;
    block.OutputPort(9).Data = Tnet;
    block.OutputPort(10).Data = Tneg;

end

function Derivatives(block)

    p = block.DialogPrm(1).Data;

    Vd     = block.InputPort(1).Data;
    Vq     = block.InputPort(2).Data;
    Tl_ext = block.InputPort(3).Data;

    x = block.ContStates.Data;

    Id   = x(1);
    Iq   = x(2);
    Wm   = x(3);
    Th_m = x(4);

    We = p.P * Wm;

    % Flujo de los imanes permanentes
    lambda_m = p.Ke / p.P;

    % =========================
    % Dinámica eléctrica dq
    % =========================
    dId = (Vd - p.R*Id + We*p.Lq*Iq) / p.Ld;

    dIq = (Vq - p.R*Iq ...
          - We*(p.Ld*Id + lambda_m)) / p.Lq;

    % =========================
    % Dinámica mecánica
    % =========================
    [~, ~, ~, ~, Tnet] = ...
        calculate_mechanics(p, Id, Iq, Wm, Th_m, Tl_ext);

    dWm  = Tnet / p.J_real;
    dThm = Wm;

    block.Derivatives.Data = [dId; dIq; dWm; dThm];

end