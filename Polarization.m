function pol = Polarization(E0, beta)
%PW_POLARIZATION  Intrinsic polarization of a plane wave.

%   Inputs:
%     E0   : 3x1 complex vector, electric-field phasor
%     beta : 3x1 real vector, propagation vector β (rad/m)
%            (only its direction is used)

%   Output struct 'pol' with fields:
%     pol.type  : 'linear', 'circular', or 'elliptic'
%     pol.hand  : 'LH', 'RH', or ''  (for linear)
%     pol.AR    : axial ratio ≥ 1  (major/minor; Inf for linear)
%     pol.psi   : orientation angle (rad) of major axis
%     pol.chi   : ellipticity angle (rad)

    E0   = E0(:);
    beta = beta(:);
    k_hat = beta / norm(beta);   % direction of propagation

    % 1. Remove any longitudinal component (safety)
    Et = E0 - (dot(k_hat, E0))*k_hat;   % transverse part

    % 2. Build an orthonormal basis {e1, e2} ⟂ k_hat 
    if abs(k_hat(1)) < 0.9
        tmp = [1;0;0];
    else
        tmp = [0;1;0];
    end
    e1 = cross(k_hat, tmp);  e1 = e1 / norm(e1);
    e2 = cross(k_hat, e1);   % already ⟂ both

    % Components of E in this basis
    E1 = e1' * Et;    % complex scalar
    E2 = e2' * Et;

    % 3. Handle degenerate cases (purely linear along e1 or e2)
    eps = 1e-9;
    if abs(E1) < eps && abs(E2) < eps
        pol.type = 'linear';
        pol.hand = '';
        pol.AR   = Inf;
        pol.psi  = 0;
        pol.chi  = 0;
        return;
    elseif abs(E1) < eps
        % purely along e2
        pol.type = 'linear';
        pol.hand = '';
        pol.AR   = Inf;
        pol.psi  = pi/2;   % along e2
        pol.chi  = 0;
        return;
    elseif abs(E2) < eps
        % purely along e1
        pol.type = 'linear';
        pol.hand = '';
        pol.AR   = Inf;
        pol.psi  = 0;      % along e1
        pol.chi  = 0;
        return;
    end

    % 4. Amplitude ratio and phase difference
    r  = abs(E2) / abs(E1);                    % |E2| / |E1|
    d  = angle(E2) - angle(E1);                % δ

    % Normalize δ to [-pi, pi] for robustness
    d = atan2(sin(d), cos(d));

    % 5. Orientation ψ and ellipticity χ (Born & Wolf / Ulaby)
    % tan(2ψ) = 2 r cosδ / (1 - r^2)
    psi = 0.5 * atan2(2*r*cos(d), 1 - r^2);

    % sin(2χ) = 2 r sinδ / (1 + r^2)
    s2chi = 2*r*sin(d) / (1 + r^2);
    s2chi = max(min(s2chi,1),-1);   % clip for numerical safety
    chi = 0.5 * asin(s2chi);

    pol.psi = psi;
    pol.chi = chi;

    % 6. Axial ratio
    if abs(chi) < 1e-3          % ≈ 0 → linear
        pol.type = 'linear';
        pol.AR   = Inf;
        pol.hand = '';
        return;
    end

    AR = 1/abs(tan(chi));       % major/minor ≥ 1
    if AR < 1
        AR = 1/AR;
    end
    pol.AR = AR;

    % 7. Type: circular vs elliptic
    if abs(abs(chi) - pi/4) < 1e-3
        pol.type = 'circular';
    else
        pol.type = 'elliptic';
    end

    % 8. Handedness (Ulaby convention)
    if sin(d) > 0
        pol.hand = 'LH';
    elseif sin(d) < 0
        pol.hand = 'RH';
    else
        pol.hand = '';   % linear limit
    end
end
