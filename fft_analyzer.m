function out = fft_analyzer(ts, varargin)
% out = fft_analyzer(ts, 'Name',"Ia", 'TimeWindow',[t1 t2], 'Fmax',5000, ...)
%
% Ventanas:
%   'TimeWindow' : [t1 t2] en segundos (mismo sistema que ts.Time)
%                  Si se entrega, recorta exactamente a ese intervalo.
%   'Tsteady'    : analiza solo los últimos Tsteady segundos (se usa solo si TimeWindow está vacío)
%
% Otras:
%   'Fmax'       : limita espectro hasta Fmax Hz (default: Nyquist)
%   'Fref'       : frecuencia referencia para buscar fundamental (Hz)
%   'BwRel'      : banda relativa +/-BwRel alrededor de Fref (default: 0.02)
%   'BwMinHz'    : banda mínima absoluta (Hz) (default: 2)
%   'Kmax'       : armónicos 2..Kmax para THD (default: 8)
%
% Extra (ripple robusto):
%   'RobustPPPrct' : percentiles [low high] para p-p robusto (default: [0.5 99.5])
%
% Salidas importantes:
%   out.mean
%   out.ripple_pp, out.ripple_rms_ac
%   out.extrema.max, out.extrema.min, out.extrema.t_max, out.extrema.t_min
%   out.ripple_pp_robust (percentil-based)

p = inputParser;
p.addParameter('Name',"signal");
p.addParameter('MinPeakFrac', 0.05);
p.addParameter('MaxPeaks', 20);
p.addParameter('DoPlot', false);

% ---- ventana manual ----
p.addParameter('TimeWindow', [], @(v) isempty(v) || (isnumeric(v) && numel(v)==2 && v(2)>v(1)));

% ---- steady-state ----
p.addParameter('Tsteady', [], @(v) isempty(v) || (isscalar(v) && v>0));

% ---- espectro ----
p.addParameter('Fmax', [], @(v) isempty(v) || (isscalar(v) && v>0));

% ---- fundamental/armónicos ----
p.addParameter('Fref', [], @(v) isempty(v) || (isscalar(v) && v>0));
p.addParameter('BwRel', 0.02, @(v) isscalar(v) && v>0);
p.addParameter('BwMinHz', 2, @(v) isscalar(v) && v>0);
p.addParameter('Kmax', 8, @(v) isscalar(v) && v>=2);

% ---- ripple robusto ----
p.addParameter('RobustPPPrct', [0.5 99.5], @(v) isnumeric(v) && numel(v)==2 && v(2)>v(1) && v(1)>=0 && v(2)<=100);

p.parse(varargin{:});

name        = string(p.Results.Name);
minPeakFrac = p.Results.MinPeakFrac;
maxPeaks    = p.Results.MaxPeaks;
doPlot      = logical(p.Results.DoPlot);

TimeWindow  = p.Results.TimeWindow;
Tsteady     = p.Results.Tsteady;

Fmax        = p.Results.Fmax;

Fref        = p.Results.Fref;
BwRel       = p.Results.BwRel;
BwMinHz     = p.Results.BwMinHz;
Kmax        = p.Results.Kmax;

robPrct     = p.Results.RobustPPPrct;

% --- leer timeseries ---
if ~isa(ts,'timeseries')
    error('Entrada debe ser timeseries. Pasa la variable timeseries (no .Data/.Time).');
end

t = ts.Time(:);
x = ts.Data;

% si x es NxM (multicanal), toma canal 1
if ~isvector(x)
    x = x(:,1);
end
x = x(:);

% --- seguridad ---
if numel(t) ~= numel(x)
    error('Time y Data no coinciden en longitud.');
end
if numel(x) < 64
    error('Muy pocas muestras para FFT confiable (N < 64).');
end

% =======================
%   RECORTE POR VENTANA
% =======================
if ~isempty(TimeWindow)
    t1 = TimeWindow(1);
    t2 = TimeWindow(2);

    % clamp a rango disponible
    t1 = max(t1, t(1));
    t2 = min(t2, t(end));

    if t2 <= t1
        error('TimeWindow inválida luego de ajustar a datos: [%.6g %.6g]', t1, t2);
    end

    idx = (t >= t1) & (t <= t2);
    if nnz(idx) < 64
        error('TimeWindow deja muy pocas muestras (N < 64). Elige una ventana más larga.');
    end

    t = t(idx);
    x = x(idx);

elseif ~isempty(Tsteady)
    % tramo estacionario: últimos Tsteady segundos
    tEnd = t(end);
    idx = t >= (tEnd - Tsteady);

    if nnz(idx) < 64
        error('Tsteady deja muy pocas muestras (N < 64). Usa un Tsteady mayor.');
    end

    t = t(idx);
    x = x(idx);
end

% --- muestreo robusto ---
dt = diff(t);
Ts = median(dt);
Fs = 1/Ts;

nonuniform = max(abs(dt - Ts)) > 1e-6*max(1,Ts);

% --- métricas de ripple en tiempo ---
x_mean = mean(x);
x_ac   = x - x_mean;

% extremos (para diagnosticar p-p inflado)
[xmax, imax] = max(x);
[xmin, imin] = min(x);
tmax = t(imax);
tmin = t(imin);

% ripple robusto por percentiles
x_hi = prctile(x, robPrct(2));
x_lo = prctile(x, robPrct(1));
pp_rob = x_hi - x_lo;

out = struct();
out.name = name;

out.t = t;
out.x = x;

out.mean = x_mean;

out.ripple_pp     = xmax - xmin;     % pico-a-pico raw
out.ripple_rms_ac = rms(x_ac);       % RMS sin DC
out.ripple_pp_robust = pp_rob;       % p-p robusto (percentiles)

out.nonuniform_time = nonuniform;

out.extrema = struct();
out.extrema.max   = xmax;
out.extrema.min   = xmin;
out.extrema.t_max = tmax;
out.extrema.t_min = tmin;
out.extrema.pp_raw = out.ripple_pp;
out.extrema.robust_percentiles = robPrct;
out.extrema.p_lo = x_lo;
out.extrema.p_hi = x_hi;

out.window = struct();
out.window.TimeWindow = TimeWindow;
out.window.Tsteady    = Tsteady;
out.window.t1_used    = t(1);
out.window.t2_used    = t(end);
out.window.N          = numel(x);

% --- FFT ---
N  = numel(x_ac);
Nh = floor(N/2);

% ventana Hann
n = (0:N-1)';
w = 0.5*(1 - cos(2*pi*n/(N-1)));
w_gain = mean(w);

xw = x_ac .* w;

Y  = fft(xw);
P2 = abs(Y/N) / w_gain;
P1 = P2(1:Nh+1);
if Nh > 1
    P1(2:end-1) = 2*P1(2:end-1);
end
f = Fs*(0:Nh)/N;

% limitar Fmax
if ~isempty(Fmax)
    idxF = f <= min(Fmax, Fs/2);
    f  = f(idxF);
    P1 = P1(idxF);
end

out.Ts = Ts;
out.Fs = Fs;
out.f  = f;
out.P1 = P1;

% --- picos (umbral relativo al máximo, excluyendo DC) ---
P1_noDC = P1;
if ~isempty(P1_noDC)
    P1_noDC(1) = 0;
end

if isempty(P1_noDC)
    threshold = 0;
else
    threshold = max(P1_noDC) * minPeakFrac;
end

pkf = [];
pkm = [];
for k = 2:length(P1)-1
    if P1(k) > P1(k-1) && P1(k) > P1(k+1) && P1(k) > threshold
        pkf(end+1,1) = f(k); %#ok<AGROW>
        pkm(end+1,1) = P1(k); %#ok<AGROW>
    end
end

if ~isempty(pkm)
    [pkm, idx] = sort(pkm, 'descend');
    pkf = pkf(idx);

    nkeep = min(maxPeaks, numel(pkm));
    pkf = pkf(1:nkeep);
    pkm = pkm(1:nkeep);
end

out.peaks = struct();
out.peaks.freq  = pkf;
out.peaks.mag   = pkm;
out.peaks.table = table(pkf, pkm, 'VariableNames', {'Frecuencia_Hz','Magnitud'});

out.f_base = ternary(~isempty(pkf), pkf(1), NaN);

% --- fundamental cerca de Fref + armónicos + THD ---
out.Fref = Fref;
out.fund = struct('freq_Hz',NaN,'mag',NaN);
out.harmonics = table();
out.THD = NaN;

if ~isempty(Fref) && isfinite(Fref) && Fref > 0
    bw = max(BwMinHz, BwRel*Fref);
    inBand = (f >= (Fref-bw)) & (f <= (Fref+bw));

    if any(inBand)
        [fundMag, imax] = max(P1(inBand));
        fBand = f(inBand);
        fundHz = fBand(imax);
        out.fund.freq_Hz = fundHz;
        out.fund.mag = fundMag;
    end

    harmK  = (2:Kmax).';
    harmHz = nan(size(harmK));
    harmMag= nan(size(harmK));

    for ii = 1:numel(harmK)
        kH = harmK(ii);
        fk = kH*Fref;
        bwk = max(BwMinHz, BwRel*fk);

        inBandK = (f >= (fk-bwk)) & (f <= (fk+bwk));
        if any(inBandK)
            [a, imax2] = max(P1(inBandK));
            fBandK = f(inBandK);
            harmHz(ii)  = fBandK(imax2);
            harmMag(ii) = a;
        end
    end

    out.harmonics = table(harmK, harmHz, harmMag, ...
        'VariableNames', {'k','Frecuencia_Hz','Magnitud'});

    if isfinite(out.fund.mag) && out.fund.mag > 0
        valid = ~isnan(harmMag);
        out.THD = sqrt(sum(harmMag(valid).^2)) / out.fund.mag;
    end
end

% --- plot opcional ---
if doPlot
    figure('Name', "FFT - " + name);
    plot(f, P1); grid on;
    xlabel('Frecuencia (Hz)'); ylabel('Magnitud');
    title("FFT de " + name + " | [" + out.window.t1_used + ", " + out.window.t2_used + "] s");

   
end

end

% helper simple (para MATLAB sin inline if)
function y = ternary(cond, a, b)
if cond, y = a; else, y = b; end
end
