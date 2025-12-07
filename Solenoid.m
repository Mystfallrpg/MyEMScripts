function sol = Solenoid(mu_r, N, A, ell, varargin)

%   sol = Solenoid(mu_r, N, A, ell, 'I', I)
%   sol = Solenoid(mu_r, N, A, ell, 'B', B)

%   Inputs:
%     mu_r : relative permeability of core
%     N    : number of turns
%     A    : cross-section area (m^2)
%     ell  : length of solenoid (m)

%   Optional (name–value, exactly one of them must be given):
%     'I'   : current (A)
%     'B'   : magnetic flux density inside solenoid (T)

%   Output struct 'sol' fields:
%     sol.mu      : absolute permeability (H/m)
%     sol.n       : turns per meter (1/m)
%     sol.I       : current (A)
%     sol.B       : B inside solenoid (T)
%     sol.L       : inductance (H)
%     sol.Phi     : flux per turn (Wb)
%     sol.lambda  : flux linkage N*Phi (Wb·turn)
%     sol.W       : stored magnetic energy (J)
%     sol.u       : energy density B^2/(2*mu) (J/m^3)

    % parse optional arguments
    I_given = [];
    B_given = [];

    if ~isempty(varargin)
        for k = 1:2:numel(varargin)
            name = lower(varargin{k});
            val  = varargin{k+1};
            switch name
                case 'i'
                    I_given = val;
                case 'b'
                    B_given = val;
                otherwise
                    error('Unknown option "%s". Use ''I'' or ''B''.', varargin{k});
            end
        end
    end

    if ~isempty(I_given) && ~isempty(B_given)
        error('Specify either I or B, not both.');
    elseif isempty(I_given) && isempty(B_given)
        error('You must specify either ''I'',I or ''B'',B.');
    end

    % constants and basic geometry
    mu0 = 4*pi*1e-7;
    mu  = mu0 * mu_r;
    n   = N / ell;       % turns per meter

    % find I and B consistently
    if ~isempty(I_given)
        I = I_given;
        B = mu * n * I;          % B = mu * n * I
    else
        B = B_given;
        I = B / (mu * n);        % I = B / (mu * n)
    end

    % inductance
    L = mu * N^2 * A / ell;

    % flux per turn, flux linkage
    Phi    = B * A;
    lambda = N * Phi;

    % stored energy and energy density
    W = 0.5 * L * I^2;
    u = B^2 / (2*mu);            % energy density

    % pack results
    sol.mu     = mu;
    sol.n      = n;
    sol.I      = I;
    sol.B      = B;
    sol.L      = L;
    sol.Phi    = Phi;
    sol.lambda = lambda;
    sol.W      = W;
    sol.u      = u;
end
