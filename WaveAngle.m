function [theta_i, theta_t, tir] = WaveAngle(beta, eps_r1, eps_r2, plane)
% TransmissionAngle  Incidence and transmission angles from propagation vector.

%   beta   : 3x1 real propagation vector (from exponent in e^{-j beta·r})
%   eps_r1 : relative permittivity of medium 1 (incident side)
%   eps_r2 : relative permittivity of medium 2 (transmitted side)
%   plane  : 'xy', 'yz', or 'xz' interface plane

%   theta_i : incidence angle (deg) w.r.t. normal
%   theta_t : transmission angle (deg) w.r.t. normal (NaN if TIR)
%   tir     : logical, true if total internal reflection

    beta = beta(:);

    switch lower(plane)
        case 'xz', n_hat = [0; 1; 0];  % y = 0 plane
        case 'xy', n_hat = [0; 0; 1];  % z = 0 plane
        case 'yz', n_hat = [1; 0; 0];  % x = 0 plane
        otherwise
            error('Unknown plane. Use ''xy'', ''yz'', or ''xz''.');
    end

    % normalize
    k_hat = beta / norm(beta);
    n_hat = n_hat / norm(n_hat);

    % incidence angle (always 0..90 deg)
    cos_theta_i = abs(dot(k_hat, n_hat));
    cos_theta_i = max(min(cos_theta_i, 1), -1);
    theta_i = acosd(cos_theta_i);
    sin_theta_i = sind(theta_i);

    % refractive indices (μr = 1 assumed)
    n1 = sqrt(eps_r1);
    n2 = sqrt(eps_r2);

    % Snell's law
    arg = (n1/n2) * sin_theta_i;

    tir = false;
    if arg > 1
        theta_t = NaN;
        tir = true;
        return;
    end

    theta_t = asind(arg);
end
