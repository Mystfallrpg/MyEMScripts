function result = PowerFraction(kind, pol, Ephasor, k_hat, n_hat, ...
                                eps_r1, mu_r1, eps_r2, mu_r2)

% Compute fraction of incident power that is REFLECTED or
% TRANSMITTED at a planar interface between two lossless media.
%
% kind    : 'reflected' or 'transmitted'
% pol     : 'TE', 'TM', or 'mixed'
% Ephasor : 3x1 complex column vector (incident E field in medium 1)
% k_hat   : 3x1 real unit vector (incident propagation direction)
% n_hat   : 3x1 real unit vector (surface normal from medium 1 -> 2)
% eps_r1, mu_r1 : relative permittivity/permeability of medium 1
% eps_r2, mu_r2 : relative permittivity/permeability of medium 2
%
% result fields:
%   .R_total, .T_total   : total power reflection / transmission coeff.
%   .value               : chosen fraction (R or T, depending on 'kind')
%   .percent             : same in percent
%   .R_TE, .R_TM         : TE/TM power reflection coeffs
%   .theta_i, .theta_t   : angles (deg)
%   .f_TE, .f_TM         : TE/TM incident power fractions (for 'mixed')
%   .tir                 : true if total internal reflection


    % normalize inputs 
    kind    = lower(kind);
    pol     = upper(pol);
    Ephasor = Ephasor(:);
    k_hat   = k_hat(:) / norm(k_hat);
    n_hat   = n_hat(:) / norm(n_hat);

    % --- constants (eta0 cancels in R/T for lossless, kept for clarity) ---
    eta0 = 376.730313668;

    % intrinsic impedances
    eta1 = eta0 * sqrt(mu_r1/eps_r1);
    eta2 = eta0 * sqrt(mu_r2/eps_r2);

    % refractive indices (for Snell)
    n1 = sqrt(eps_r1*mu_r1);
    n2 = sqrt(eps_r2*mu_r2);

    % incidence angle (0..pi/2) 
    cos_theta_i = abs(dot(k_hat, n_hat));
    cos_theta_i = max(min(cos_theta_i,1),-1);  % clip numerically
    theta_i = acos(cos_theta_i);
    sin_theta_i = sin(theta_i);

    % transmission angle (Snell)
    sin_theta_t = (n1/n2) * sin_theta_i;

    % total internal reflection check
    if abs(sin_theta_t) > 1
        % TIR: |Gamma| = 1 for both pol, no transmitted power
        R_TE = 1;
        R_TM = 1;
        R_total = 1;
        T_total = 0;

        % pick reflected / transmitted
        switch kind
            case 'reflected'
                value = R_total;
            case 'transmitted'
                value = T_total;
            otherwise
                error('Unknown kind "%s": use reflected or transmitted.', kind);
        end

        result.R_total = R_total;
        result.T_total = T_total;
        result.value   = value;
        result.percent = 100*value;
        result.R_TE    = R_TE;
        result.R_TM    = R_TM;
        result.theta_i = theta_i*180/pi;
        result.theta_t = NaN;
        result.f_TE    = NaN;
        result.f_TM    = NaN;
        result.tir     = true;
        return;
    end

    theta_t = asin(sin_theta_t);
    cos_theta_t = cos(theta_t);

    % Fresnel amplitude coefficients
    Gamma_TE = (eta2*cos_theta_i - eta1*cos_theta_t) / ...
               (eta2*cos_theta_i + eta1*cos_theta_t);

    Gamma_TM = (eta2*cos_theta_t - eta1*cos_theta_i) / ...
               (eta2*cos_theta_t + eta1*cos_theta_i);

    R_TE = abs(Gamma_TE)^2;
    R_TM = abs(Gamma_TM)^2;

    % TE/TM split of the actual field (needed for 'mixed')
    % TE direction = perpendicular to plane of incidence
    e_TE = cross(k_hat, n_hat);
    if norm(e_TE) < 1e-12
        % normal incidence → TE/TM not uniquely defined
        e_TE = [0;0;0];
    else
        e_TE = e_TE / norm(e_TE);
    end

    % projection of E on TE direction
    E_TE_vec = dot(Ephasor, e_TE) * e_TE;
    E_TM_vec = Ephasor - E_TE_vec;

    P_TE_inc = norm(E_TE_vec)^2;   % ∝ incident TE power
    P_TM_inc = norm(E_TM_vec)^2;   % ∝ incident TM power
    P_sum    = P_TE_inc + P_TM_inc;

    if P_sum < 1e-18
        f_TE = 0;
        f_TM = 0;
    else
        f_TE = P_TE_inc / P_sum;
        f_TM = P_TM_inc / P_sum;
    end

    % choose how to combine TE/TM, based on 'pol'
    switch pol
        case 'TE'
            R_total = R_TE;
        case 'TM'
            R_total = R_TM;
        case 'MIXED'
            if P_sum < 1e-18
                R_total = 0;
            else
                % weighted by incident TE/TM powers
                R_total = (P_TE_inc*R_TE + P_TM_inc*R_TM) / P_sum;
            end
        otherwise
            error('Unknown pol "%s": use TE, TM, or mixed.', pol);
    end

    T_total = 1 - R_total;

    % pick reflected / transmitted
    switch kind
        case 'reflected'
            value = R_total;
        case 'transmitted'
            value = T_total;
        otherwise
            error('Unknown kind "%s": use reflected or transmitted.', kind);
    end

    % pack result
    result.R_total = R_total;
    result.T_total = T_total;
    result.value   = value;
    result.percent = 100*value;
    result.R_TE    = R_TE;
    result.R_TM    = R_TM;
    result.theta_i = theta_i*180/pi;
    result.theta_t = theta_t*180/pi;
    result.f_TE    = f_TE;
    result.f_TM    = f_TM;
    result.tir     = false;
end
