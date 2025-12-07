function power = WavePower(E0, beta, medium, n_hat, A)
%WavePower  Time-average power quantities for a plane wave.

%   power = WavePower(E0, beta, medium)
%   power = WavePower(E0, beta, medium, n_hat)
%   power = WavePower(E0, beta, medium, n_hat, A)

%   Inputs:
%     E0     : 3x1 complex phasor of E (V/m)
%     beta   : 3x1 real propagation vector from exp(-j beta·r)
%     medium : what your E2H() accepts (eta or struct)
%     n_hat  : (optional,[]) 3x1 real surface normal
%     A      : (optional,[]) surface area (m^2)

%   Output: struct 'power' with fields:
%     power.H0         : 3x1 complex H phasor (A/m)
%     power.Savg       : 3x1 real time-average Poynting vector (W/m^2)
%     power.Savg_mag   : scalar |<S>| (W/m^2)
%     power.S_normal   : (if n_hat given) <S> · n_hat (W/m^2)
%     power.P          : (if n_hat & A given) total power through surface (W)

    % Ensure column vectors
    E0   = E0(:);
    beta = beta(:);

    % 1) H from E
    H0 = E2H(E0, beta, medium);

    % 2) Time-average Poynting vector
    Savg = 0.5 * real(cross(E0, conj(H0)));
    Savg_mag = norm(Savg);

    % Fill basic outputs
    power.H0       = H0;
    power.Savg     = Savg;
    power.Savg_mag = Savg_mag;

    % 3) If a surface normal is provided
    if nargin >= 4 && ~isempty(n_hat)
        n_hat = n_hat(:);
        n_hat = n_hat / norm(n_hat);

        S_normal = dot(Savg, n_hat);
        power.S_normal = S_normal;  % W/m^2

        % 4) If an area is also provided
        if nargin >= 5 && ~isempty(A)
            power.P = S_normal * A; % W
        end
    end
end
