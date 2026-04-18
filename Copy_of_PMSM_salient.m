    function PMSM_salient(block)
        setup(block);
    end
    
    function setup(block)
        block.NumInputPorts  = 3; % Va, Vb, Tl
        block.NumOutputPorts = 7; % Ia, Ib, Wm, We, Th_m, Th_e, Te 
        
        block.SetPreCompInpPortInfoToDynamic;
        block.SetPreCompOutPortInfoToDynamic;
    
        block.NumDialogPrms = 1; 
        block.NumContStates = 4; 
        block.SampleTimes = [0 0]; 
    
        block.RegBlockMethod('InitializeConditions', @InitConditions);
        block.RegBlockMethod('Outputs',              @Outputs);
        block.RegBlockMethod('Derivatives',          @Derivatives);
    end
    
    function InitConditions(block)
        block.ContStates.Data = [0; 0; 0; 0]; 
    end
    
    function Outputs(block)
        p = block.DialogPrm(1).Data;
        states = block.ContStates.Data;
        
        Ia = states(1);
        Ib = states(2);
        Wm = states(3);
        Th_m = states(4);
        
        % --- Re-calculamos variables dependientes para las salidas ---
        Th_e = p.P * Th_m;
        Th_e_2 = Th_e*2;
        We = p.P * Wm;
        L2 = (p.Ld - p.Lq)/2;
    
        % Cálculo de Torque (Debe ser idéntico al de Derivatives)
        % Te = Torque_Imán + Torque_Reluctancia
        Te = p.Ke*(Ib*cos(Th_e) - Ia*sin(Th_e)) + ...
             p.P*L2*(2*Ia*Ib*cos(Th_e_2) - (Ia^2 - Ib^2)*sin(Th_e_2));
    
        % Asignación a los 7 puertos
        block.OutputPort(1).Data = Ia;
        block.OutputPort(2).Data = Ib;
        block.OutputPort(3).Data = Wm;
        block.OutputPort(4).Data = We;
        block.OutputPort(5).Data = Th_m;
        block.OutputPort(6).Data = Th_e;
        block.OutputPort(7).Data = Te; % <--- Salida de Torque disponible
    end
    
    function Derivatives(block)
        params = block.DialogPrm(1).Data;
        Va = block.InputPort(1).Data;
        Vb = block.InputPort(2).Data;
        Tl = block.InputPort(3).Data;
        
        Ia = block.ContStates.Data(1);
        Ib = block.ContStates.Data(2);
        Wm = block.ContStates.Data(3);
        Th_m = block.ContStates.Data(4);
    
        Th_e = params.P * Th_m;
        Th_e_2 = 2*Th_e;
        L0 = (params.Ld + params.Lq)/2;
        L2 = (params.Ld - params.Lq)/2;
    
        % Matriz de Inductancia
        La  = L0 + L2*cos(Th_e_2);
        Lb  = L0 - L2*cos(Th_e_2);
        Lab = L2*sin(Th_e_2);
        Ld = L0 + L2;
        Lq = L0 - L2;
        detL = La * Lb - Lab^2;
    
        % Variaciones temporales de L
        dLa_dt  = -2*params.P*Wm * L2*sin(Th_e_2);
        dLb_dt  =  2*params.P*Wm * L2*sin(Th_e_2);
        dLab_dt =  2*params.P*Wm * L2*cos(Th_e_2);
    
        % Back-EMF
        Ea = -params.Ke * Wm * sin(Th_e);
        Eb =  params.Ke * Wm * cos(Th_e);
    
        % Voltajes Netos
        Va_net = Va - params.R*Ia - Ea - (Ia*dLa_dt + Ib*dLab_dt);
        Vb_net = Vb - params.R*Ib - Eb - (Ib*dLb_dt + Ia*dLab_dt);
    
        % Derivadas de Corriente
        dIa = ( Lb*Va_net - Lab*Vb_net) / detL;
        dIb = ( La*Vb_net - Lab*Va_net) / detL;
    
        % Torque
        Te = params.Ke*(Ib*cos(Th_e) - Ia*sin(Th_e)) + ...
             params.P*L2*(2*Ia*Ib*cos(Th_e_2) - (Ia^2 - Ib^2)*sin(Th_e_2));
        
        % Mecánica
        dWm = (Te - Tl - params.B_real*Wm - params.Tdm*sin(2*Th_e_2)) / params.J_real;
        dTh = Wm;
    
        block.Derivatives.Data = [dIa; dIb; dWm; dTh];
    end