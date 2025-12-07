function cap = Capacitor(kind, eps_r, varargin)
% Capacitor  Basic capacitance + energy for common geometries.

%   cap = Capacitor('parallel', eps_r, A, d, 'V', V)
%   cap = Capacitor('coax',     eps_r, a, b, L, 'Q', Q)
%   cap = Capacitor('twowire',  eps_r, r, D, L, 'V', V)

%   kind   : 'parallel', 'coax', or 'twowire'
%   eps_r  : relative permittivity of dielectric

%   Geometry arguments (after eps_r):
%     'parallel' : A (m^2), d (m)
%     'coax'     : a (m) inner radius, b (m) outer radius, l (m)
%     'twowire'  : r (m) wire radius, D (m) center spacing, l (m)

%   Optional name–value (exactly 0 or 1 of them):   
%     'V', V : voltage across capacitor (V)
%     'Q', Q : total charge magnitude (C)

%   Output struct cap:
%     cap.C       : capacitance (F)
%     cap.C_per   : capacitance per unit length (F/m) where relevant
%     cap.V       : voltage (V)  (if Q or V given)
%     cap.Q       : charge (C)   (if Q or V given)
%     cap.W       : stored energy (J), if V or Q given
%     cap.kind    : geometry string

%   Note: Two-wire uses C' = pi*eps/acosh(D/radius_ratio),
%         i.e. C' = pi*eps / acosh(D/(2*r)), valid for typical D >> r.

    kind = lower(kind);
    eps0 = 8.854187817e-12;
    eps  = eps0*eps_r;

    % parse geometry
    switch kind
        case 'parallel'
            if numel(varargin) < 2
                error('parallel: need A, d, then optional ''V''/''Q''.');
            end
            A = varargin{1};
            d = varargin{2};
            idx = 3;

            C = eps*A/d;
            C_per = NaN;  % not per-length in this geometry

        case 'coax'
            if numel(varargin) < 3
                error('coax: need a, b, L, then optional ''V''/''Q''.');
            end
            a = varargin{1};
            b = varargin{2};
            L = varargin{3};
            idx = 4;

            C   = 2*pi*eps*L / log(b/a);       % total C
            C_per = 2*pi*eps / log(b/a);       % per unit length

        case 'twowire'
            if numel(varargin) < 3
                error('twowire: need r, D, L, then optional ''V''/''Q''.');
            end
            r = varargin{1};
            D = varargin{2};
            L = varargin{3};
            idx = 4;

            % Capacitance per unit length for two-wire line
            % General form: C' = pi*eps / acosh(D/radius_ratio).
            % Here radius = r, center spacing = D → acosh(D/r).
            C_per = pi*eps / acosh(D/(2*r));
            C     = C_per * L;                 % total C

        otherwise
            error('Unknown kind "%s". Use ''parallel'', ''coax'', or ''twowire''.', kind);
    end

    % parse optional V/Q
    V_given = [];
    Q_given = [];

    while idx <= numel(varargin)
        name = lower(varargin{idx});
        val  = varargin{idx+1};
        switch name
            case 'v'
                V_given = val;
            case 'q'
                Q_given = val;
            otherwise
                error('Unknown option "%s". Use ''V'' or ''Q''.', varargin{idx});
        end
        idx = idx + 2;
    end

    if ~isempty(V_given) && ~isempty(Q_given)
        error('Specify either V or Q or neither, not both.');
    end

    % compute Q, V, W if possible
    if ~isempty(V_given)
        V = V_given;
        Q = C*V;
        W = 0.5*C*V^2;
    elseif ~isempty(Q_given)
        Q = Q_given;
        V = Q/C;
        W = 0.5*Q^2/C;
    else
        V = NaN; Q = NaN; W = NaN;
    end

    % pack output
    cap.kind  = kind;
    cap.C     = C;
    cap.C_per = C_per;
    cap.V     = V;
    cap.Q     = Q;
    cap.W     = W;
end
