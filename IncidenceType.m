function result = IncidenceType(E, k, surface)
% IncidenceType  Classify incidence as normal / TE / TM / mixed

%   E       : 3x1 complex phasor (electric field amplitude)
%   k       : 3x1 real vector (propagation vector or direction)
%   surface : 'xy', 'yz', or 'xz' for interface plane

%   result.type = "Normal" | "TE" | "TM" | "Mixed"
%   result.description = string description

    % Ensure column vectors
    E = E(:);
    k = k(:);

    % Surface normal
    switch surface
        case 'xz'     % y = 0
            n = [0; 1; 0];
        case 'xy'     % z = 0
            n = [0; 0; 1];
        case 'yz'     % x = 0
            n = [1; 0; 0];
        otherwise
            error('surface must be ''xy'', ''yz'', or ''xz''.');
    end

    % Normalize directions
    if norm(k) == 0
        error('k must be non-zero.');
    end
    k_hat = k / norm(k);
    n_hat = n / norm(n);

    % Tolerances (relative)
    tol = 1e-9;

    % 1) Normal incidence test: k ‖ n
    if norm(cross(k_hat, n_hat)) <= tol
        result.type        = "Normal";
        result.description = "Normal incidence";
        return;
    end

    % 2) Oblique incidence: build normal to plane of incidence
    Ni = cross(k_hat, n_hat);      % normal to plane of incidence
    Ni_hat = Ni / norm(Ni);

    % For complex E, we use its phasor direction
    E_mag = norm(E);
    if E_mag == 0
        result.type        = "Mixed";
        result.description = "Zero field (undefined incidence type)";
        return;
    end

    % 3) TE test: E ∥ Ni  →  E × Ni ≈ 0
    cross_ENi = norm(cross(E, Ni_hat));
    is_TE = cross_ENi <= tol * E_mag;

    % 4) TM test: E ⟂ Ni  →  E · Ni ≈ 0
    dot_ENi = abs(dot(E, Ni_hat));
    is_TM = dot_ENi <= tol * E_mag;

    if is_TE && ~is_TM
        result.type        = "TE";
        result.description = "Oblique transverse electric incidence (E ⟂ plane of incidence)";
    elseif is_TM && ~is_TE
        result.type        = "TM";
        result.description = "Oblique transverse magnetic incidence (E in plane of incidence)";
    else
        % Numerically could be close to both or neither
        result.type        = "Mixed";
        result.description = "Oblique incidence with both TE and TM components";
    end
end
