function T = table_fund_and_3harm(outList, labels)

n = numel(outList);

Caso = strings(n,1);
Fund_Hz = nan(n,1);
Fund_Amp = nan(n,1);

H2_Amp = nan(n,1);
H3_Amp = nan(n,1);
H4_Amp = nan(n,1);

H2_pct = nan(n,1);
H3_pct = nan(n,1);
H4_pct = nan(n,1);

for i = 1:n
    o = outList{i};
    Caso(i) = string(labels{i});

    if isfield(o,'fund') && isfinite(o.fund.mag)
        Fund_Hz(i)  = o.fund.freq_Hz;
        Fund_Amp(i) = o.fund.mag;
    end

    if isfield(o,'harmonics') && ~isempty(o.harmonics)
        H = o.harmonics;

        % buscar k=2,3,4
        for k = 2:4
            idx = find(H.k == k);
            if ~isempty(idx)
                amp = H.Magnitud(idx);
                switch k
                    case 2
                        H2_Amp(i) = amp;
                        H2_pct(i) = 100*amp/Fund_Amp(i);
                    case 3
                        H3_Amp(i) = amp;
                        H3_pct(i) = 100*amp/Fund_Amp(i);
                    case 4
                        H4_Amp(i) = amp;
                        H4_pct(i) = 100*amp/Fund_Amp(i);
                end
            end
        end
    end
end

T = table(Caso, ...
          Fund_Hz, Fund_Amp, ...
          H2_Amp, H2_pct, ...
          H3_Amp, H3_pct, ...
          H4_Amp, H4_pct);

end
function [outputArg1,outputArg2] = untitled(inputArg1,inputArg2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
outputArg1 = inputArg1;
outputArg2 = inputArg2;
end