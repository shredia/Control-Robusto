function PMSM_salient_dq(block)
    setup(block);
end

function setup(block)
    block.NumInputPorts  = 3; % Vd, Vq, Tl
    block.NumOutputPorts = 8; % Id, Iq, Wm, We, Th_m, Th_e, Te,Tdtm

    block.SetPreCompInpPortInfoToDynamic;
    block.SetPreCompOutPortInfoToDynamic;

    for k = 1:block.NumInputPorts
        block.InputPort(k).Dimensions        = 1;
        block.InputPort(k).DirectFeedthrough = true;
    end

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
    % [Id; Iq; Wm; Th_m]
    block.ContStates.Data = [0; 0; 0; 0];
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

    % Para modelo bifásico dq:
    % Te = Te_PM + Te_rel
    % Te_PM  = Ke * Iq
    % Te_rel = P*(Ld-Lq)*Id*Iq
    Te = p.Ke*Iq - p.P*(p.Ld - p.Lq)*Id*Iq;
    Tdtm = p.Tdm * sin(p.Nr*Th_m + p.Phi); %%Revisar si es N_steps Ó Nr (200 valles energeóticos o por pares de polos)


    block.OutputPort(1).Data = Id;
    block.OutputPort(2).Data = Iq;
    block.OutputPort(3).Data = Wm;
    block.OutputPort(4).Data = We;
    block.OutputPort(5).Data = Th_m;
    block.OutputPort(6).Data = Th_e;
    block.OutputPort(7).Data = Te;
    block.OutputPort(8).Data = Tdtm;
end

function Derivatives(block)
    p = block.DialogPrm(1).Data;

    Vd = block.InputPort(1).Data;
    Vq = block.InputPort(2).Data;
    Tl = block.InputPort(3).Data;

    x = block.ContStates.Data;

    Id   = x(1);
    Iq   = x(2);
    Wm   = x(3);
    Th_m = x(4);

    We = p.P * Wm;

    % =========================
    % Modelo eléctrico dq
    % =========================
    % Vd = R*Id + Ld*dId - We*Lq*Iq
    % Vq = R*Iq + Lq*dIq + We*(Ld*Id + lambda_m)
    %
    % Con tu notación:
    % Ke = P*lambda_m  => lambda_m = Ke/P
    %
    % Entonces:
    lambda_m = p.Ke / p.P;

    dId = (Vd - p.R*Id + We*p.Lq*Iq) / p.Ld;
    dIq = (Vq - p.R*Iq - We*(p.Ld*Id + lambda_m)) / p.Lq;

    % =========================
    % Torque electromagnético
    % =========================
    Te = p.Ke*Iq + p.P*(p.Lq - p.Ld)*Id*Iq;

    % =========================
    % Torque detent
    % =========================
    % Recomendado: usar función periódica de ángulo mecánico.
    % Ejemplo simple:
    % Tdet = p.Tdm * sin(p.Nr * Th_m + p.phi_tdm);
    %
    % Si no tienes Nr o phi_tdm, usa phi=0 y define Nr en params.
    

    Tdet = p.Tdm * sin(p.Nr*Th_m + p.Phi); %%Revisar si es N_steps Ó Nr (200 valles energeóticos o por pares de polos)

    % =========================
    % Mecánica
    % =========================
    W_eps = 0.1;

    Tl = Tl*tanh(Wm/W_eps);
    dWm  = (Te - Tl - p.B_real*Wm - Tdet) / p.J_real;
    dThm = Wm;

    block.Derivatives.Data = [dId; dIq; dWm; dThm];
end