fe_ref = Wm_ref.data(1)*P/(2*pi)
Kmax = 5; %%Primero K Armónicos

t_inicial = 0.2;
t_final  = 0.8;
Fmax = 2000;
plot = 1;
out_Ia_M1 = fft_analyzer(Ia_planta_M1, 'Name',"Ia_M1", 'MaxPeaks',20, 'MinPeakFrac',0.05,'Fref',fe_ref,'Fmax',Fmax,'Kmax',Kmax,'TimeWindow',[t_inicial t_final],'DoPlot',plot);
out_Ib_M1 = fft_analyzer(Ib_planta_M1, 'Name',"Ib_M1", 'MaxPeaks',20, 'MinPeakFrac',0.05,'Fref',fe_ref,'Fmax',Fmax,'Kmax',Kmax,'TimeWindow',[t_inicial t_final],'DoPlot',plot);

out_Ia_M2 = fft_analyzer(Ia_planta_M2, 'Name',"Ia_M2", 'MaxPeaks',20, 'MinPeakFrac',0.05,'Fref',fe_ref,'Fmax',Fmax,'Kmax',Kmax,'TimeWindow',[t_inicial t_final],'DoPlot',plot);
out_Ib_M2 = fft_analyzer(Ib_planta_M2, 'Name',"Ib_M2", 'MaxPeaks',20, 'MinPeakFrac',0.05,'Fref',fe_ref,'Fmax',Fmax,'Kmax',Kmax,'TimeWindow',[t_inicial t_final],'DoPlot',plot);
%%En rad/s
out_We_M1 = fft_analyzer(We_planta_M1, 'Name',"We_M1", 'MaxPeaks',20, 'MinPeakFrac',0.05,'Fref',fe_ref,'Fmax',Fmax,'Kmax',Kmax,'TimeWindow',[t_inicial t_final],'DoPlot',plot);
out_We_M2 = fft_analyzer(We_planta_M2, 'Name',"We_M2", 'MaxPeaks',20, 'MinPeakFrac',0.05,'Fref',fe_ref,'Fmax',Fmax,'Kmax',Kmax,'TimeWindow',[t_inicial t_final],'DoPlot',plot);

fe_M1 = We_planta_M1.mean/(2*pi);
fe_M2 = We_planta_M2.mean/(2*pi);
%%En rad/s
out_Wm_M1 = fft_analyzer(Wm_planta_M1, 'Name',"Wm_M1", 'MaxPeaks',20, 'MinPeakFrac',0.05,'Fref',fe_ref,'Fmax',Fmax,'Kmax',Kmax,'TimeWindow',[t_inicial t_final],'DoPlot',plot);
out_Wm_M2 = fft_analyzer(Wm_planta_M2, 'Name',"Wm_M2", 'MaxPeaks',20, 'MinPeakFrac',0.05,'Fref',fe_ref,'Fmax',Fmax,'Kmax',Kmax,'TimeWindow',[t_inicial t_final],'DoPlot',plot);

out_Te_M1 = fft_analyzer(Te_planta_M1, 'Name',"Te_M1", 'MaxPeaks',20, 'MinPeakFrac',0.05,'Fref',fe_ref,'Fmax',Fmax,'Kmax',Kmax,'TimeWindow',[t_inicial t_final],'DoPlot',plot);
out_Te_M2 = fft_analyzer(Te_planta_M2, 'Name',"Te_M2", 'MaxPeaks',20, 'MinPeakFrac',0.05,'Fref',fe_ref,'Fmax',Fmax,'Kmax',Kmax,'TimeWindow',[t_inicial t_final],'DoPlot',plot);

Tsteady = 0.20;      % segundos finales para calcular fe_ref desde We
bwRel = 0.02;        % banda relativa +/-2% alrededor del target
bwMinHz = 2;         % banda mínima absoluta (Hz)
Kmax = 8;            % armónicos para THD (2..Kmax)

fprintf("fe_ref M1 = %.3f Hz | fe_ref M2 = %.3f Hz\n", fe_M1, fe_M2);


fprintf("THD Ia M1 = %.2f %%\n", 100*out_Ia_M1.THD);
fprintf("THD Ib M1 = %.2f %%\n", 100*out_Ib_M1.THD);
fprintf("THD Ia M2 = %.2f %%\n", 100*out_Ia_M2.THD);
fprintf("THD Ib M2 = %.2f %%\n", 100*out_Ia_M2.THD);



T_Ia = table_fund_and_3harm( ...
    {out_Ia_M1, out_Ib_M1,out_Ia_M2, out_Ib_M2}, ...
    {'Ia M1','Ib M1','Ia M2','Ib M2'});

disp(T_Ia)


fprintf("\n---------------------------------------------Análisis Te ---------------------------------------------\n")
T_Te = table_fund_and_3harm( ...
    {out_Te_M1, out_Te_M2}, ...
    {'Te M1','Te M2'});

disp(T_Te)

fprintf("mean Te M1 = %.2f \n", out_Te_M1.mean);
fprintf("mean Te M2 = %.2f \n", out_Te_M2.mean);
fprintf("ripple P-P Te M1 = %.2f \n", out_Te_M1.ripple_pp);
fprintf("ripple P-P Te M2 = %.2f \n", out_Te_M2.ripple_pp);
fprintf("ripple rms Te M1 = %.2f \n", out_Te_M1.ripple_rms_ac);
fprintf("ripple rms Te M2 = %.2f \n", out_Te_M2.ripple_rms_ac);
fprintf("Max M1= %.4f en t = %.6f s\n", out_Te_M1.extrema.max, out_Te_M1.extrema.t_max);
fprintf("Max M2= %.4f en t = %.6f s\n", out_Te_M2.extrema.max, out_Te_M1.extrema.t_max);
fprintf("Min M1= %.4f en t = %.6f s\n", out_Te_M1.extrema.min, out_Te_M1.extrema.t_min);
fprintf("Min M2= %.4f en t = %.6f s\n", out_Te_M2.extrema.min, out_Te_M1.extrema.t_min);
fprintf("P-P raw M1= %.4f\n", out_Te_M1.extrema.pp_raw);
fprintf("P-P raw M2= %.4f\n", out_Te_M2.extrema.pp_raw);
fprintf("P-P robust M1= %.4f\n", out_Te_M1.ripple_pp_robust);
fprintf("P-P robust M2= %.4f\n", out_Te_M2.ripple_pp_robust);



fprintf("\n---------------------------------------------Análisis We ---------------------------------------------\n")
T_We = table_fund_and_3harm( ...
    {out_We_M1, out_We_M2}, ...
    {'We M1','We M2'});

disp(T_We)

fprintf("mean We M1 = %.2f \n", out_We_M1.mean);
fprintf("mean We M2 = %.2f \n", out_We_M2.mean);
fprintf("ripple P-P We M1 = %.2f \n", out_We_M1.ripple_pp);
fprintf("ripple P-P We M2 = %.2f \n", out_We_M2.ripple_pp);
fprintf("ripple rms We M1 = %.2f \n", out_We_M1.ripple_rms_ac);
fprintf("ripple rms We M2 = %.2f \n", out_We_M2.ripple_rms_ac);
fprintf("Max M1= %.4f en t = %.6f s\n", out_We_M1.extrema.max, out_We_M1.extrema.t_max);
fprintf("Max M2= %.4f en t = %.6f s\n", out_We_M2.extrema.max, out_We_M2.extrema.t_max);
fprintf("Min M1= %.4f en t = %.6f s\n", out_We_M1.extrema.min, out_We_M1.extrema.t_min);
fprintf("Min M2= %.4f en t = %.6f s\n", out_We_M2.extrema.min, out_We_M2.extrema.t_min);
fprintf("P-P raw M1= %.4f\n", out_We_M1.extrema.pp_raw);
fprintf("P-P raw M2= %.4f\n", out_We_M2.extrema.pp_raw);
fprintf("P-P robust M1= %.4f\n", out_We_M1.ripple_pp_robust);
fprintf("P-P robust M2= %.4f\n", out_We_M2.ripple_pp_robust);

fprintf("\n---------------------------------------------Análisis Wm ---------------------------------------------\n")
T_Wm = table_fund_and_3harm( ...
    {out_Wm_M1, out_Wm_M2}, ...
    {'Wm M1','Wm M2'});

disp(T_Wm)

fprintf("mean Wm M1 = %.2f \n", out_Wm_M1.mean);
fprintf("mean Wm M2 = %.2f \n", out_Wm_M2.mean);
fprintf("ripple P-P Wm M1 = %.2f \n", out_Wm_M1.ripple_pp);
fprintf("ripple P-P Wm M2 = %.2f \n", out_Wm_M2.ripple_pp);
fprintf("ripple rms Wm M1 = %.2f \n", out_Wm_M1.ripple_rms_ac);
fprintf("ripple rms Wm M2 = %.2f \n", out_Wm_M2.ripple_rms_ac);
fprintf("Max M1= %.4f en t = %.6f s\n", out_Wm_M1.extrema.max, out_Wm_M1.extrema.t_max);
fprintf("Max M2= %.4f en t = %.6f s\n", out_Wm_M2.extrema.max, out_Wm_M2.extrema.t_max);
fprintf("Min M1= %.4f en t = %.6f s\n", out_Wm_M1.extrema.min, out_Wm_M1.extrema.t_min);
fprintf("Min M2= %.4f en t = %.6f s\n", out_Wm_M2.extrema.min, out_Wm_M2.extrema.t_min);
fprintf("P-P raw M1= %.4f\n", out_Wm_M1.extrema.pp_raw);
fprintf("P-P raw M2= %.4f\n", out_Wm_M2.extrema.pp_raw);
fprintf("P-P robust M1= %.4f\n", out_Wm_M1.ripple_pp_robust);
fprintf("P-P robust M2= %.4f\n", out_Wm_M2.ripple_pp_robust);