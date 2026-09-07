function [phi_deg, mag_ratio] = fase_frecuencia(t, x, y, f0)

    % x = referencia
    % y = señal que sigue a x

    t = t(:);
    x = x(:);
    y = y(:);

    % quitar componente DC
    x = x - mean(x);
    y = y - mean(y);

    % frecuencia angular
    w0 = 2*pi*f0;

    % proyección compleja en f0
    Ex = exp(-1j*w0*t);

    X = sum(x .* Ex);
    Y = sum(y .* Ex);

    % y respecto de x
    phi = angle(Y/X);

    phi_deg = rad2deg(phi);

    % llevar a [-180,180]
    phi_deg = mod(phi_deg + 180,360) - 180;

    mag_ratio = abs(Y)/abs(X);

end