function H = E2H(E, beta, medium)
% E2H  Magnetic field phasor H from electric field phasor E for a uniform plane wave.

%   Inputs:
%     E     : 3x1 complex vector, electric-field phasor (V/m)
%     beta  : 3x1 real vector, propagation vector β (rad/m) from ωt - β·r
%              (a) medium.eta  : complex intrinsic impedance (Ohm)
%                  [all other fields optional]
%              (b) medium.f    : frequency (Hz)
%                  medium.eps_r: relative permittivity
%                  medium.mu_r : relative permeability
%                  medium.sigma: conductivity (S/m)
%   Output:
%     H     : 3x1 complex vector, magnetic-field phasor (A/m)

    E    = E(:);
    beta = beta(:);

    % Direction of propagation
    k_hat = beta / norm(beta);

    % Intrinsic impedance
    if isfield(medium,'eta')
        eta = medium.eta;
    else
        eps0 = 8.854187817e-12;
        mu0  = 4*pi*1e-7;
        eps_r = medium.eps_r;
        mu_r  = medium.mu_r;
        sigma = medium.sigma;
        f     = medium.f;

        omega = 2*pi*f;
        eps   = eps0*eps_r;
        mu    = mu0*mu_r;
        eps_c = eps - 1j*sigma/omega;
        eta   = sqrt(mu./eps_c);
    end

    % Core relation
    H = (1/eta)*cross(k_hat, E);
end
