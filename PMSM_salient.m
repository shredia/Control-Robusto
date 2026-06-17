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
    
    function [Te, L_mat, dL_dTh] = calculate_physics(p, Ia, Ib, Wm, Th_m)
        % Función auxiliar interna para consistencia total entre Outputs y Derivatives
        Th_e = p.P * Th_m;
        Th_e_2 = 2 * Th_e;
        L0 = (p.Ld + p.Lq)/2;
        L2 = (p.Ld - p.Lq)/2;
        
        % 1. Matriz de Inductancia L(theta)
        L_mat = [ L0 + L2*cos(Th_e_2),      L2*sin(Th_e_2); ...
                  L2*sin(Th_e_2),      L0 - L2*cos(Th_e_2) ];
              
        % 2. Derivada de L respecto a theta eléctrica (dL/dTh_e)
        dL_dTh = [ -2*L2*sin(Th_e_2),   2*L2*cos(Th_e_2); ...
                    2*L2*cos(Th_e_2),   2*L2*sin(Th_e_2) ];
                
        % 3. Cálculo de Torque Matricial
        I_vec = [Ia; Ib];
        % Torque de reluctancia: 0.5 * I' * (dL/dTh_m) * I
        % Recordar: dL/dTh_m = P * dL/dTh_e
        Te_reluc = 0.5 * I_vec' * (p.P * dL_dTh) * I_vec;
        % Torque síncrono (imanes)
        Te_sync  = p.Ke * (Ib*cos(Th_e) - Ia*sin(Th_e));
        
        Te = Te_sync + Te_reluc;
    end
    
    function Outputs(block)
        p = block.DialogPrm(1).Data;
        x = block.ContStates.Data;
        Ia = x(1); Ib = x(2); Wm = x(3); Th_m = x(4);
        
        [Te, ~, ~] = calculate_physics(p, Ia, Ib, Wm, Th_m);
        
        block.OutputPort(1).Data = Ia;
        block.OutputPort(2).Data = Ib;
        block.OutputPort(3).Data = Wm;
        block.OutputPort(4).Data = p.P * Wm;
        block.OutputPort(5).Data = Th_m;
        block.OutputPort(6).Data = p.P * Th_m;
        block.OutputPort(7).Data = Te; 
    end
    
    function Derivatives(block)
        p  = block.DialogPrm(1).Data;
        Va = block.InputPort(1).Data;
        Vb = block.InputPort(2).Data;
        Tl = block.InputPort(3).Data;
        
        x = block.ContStates.Data;
        Ia = x(1); Ib = x(2); Wm = x(3); Th_m = x(4);
        I_vec = [Ia; Ib];
        Th_e = p.P * Th_m;
    
        % Obtener matrices consistentes
        [Te, L_mat, dL_dTh] = calculate_physics(p, Ia, Ib, Wm, Th_m);
    
        % --- Dinámica Eléctrica ---
        % Vector de Back-EMF (E)
        E_vec = [ -p.Ke * Wm * sin(Th_e); ...
                   p.Ke * Wm * cos(Th_e) ];
        
        % Variación temporal de inductancia: dL/dt = (dL/dTh_e * We)
        dL_dt = dL_dTh * (p.P * Wm);
        
        % Ecuación: L*dI/dt = V - R*I - (dL/dt)*I - E
        rhs_elec = [Va; Vb] - p.R * I_vec - dL_dt * I_vec - E_vec;
        dI_vec = L_mat \ rhs_elec; 
        
        % --- Dinámica Mecánica ---
        % Te_total incluye reluctancia y síncrono. Sumamos Detent Torque (Tdm)
        T_detent = p.Tdm * sin(p.Nr* Th_m + p.Phi); % Usualmente 4*P o P según el motor
        dWm = (Te - Tl - p.B_real * Wm - T_detent) / p.J_real;
        dTh = Wm;
    
        block.Derivatives.Data = [dI_vec(1); dI_vec(2); dWm; dTh];
    end