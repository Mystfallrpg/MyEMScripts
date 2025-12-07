function E0 = E2Phasor(Ex_cos, Ey_cos, Ez_cos, Ex_sin, Ey_sin, Ez_sin)
%E2PHASOR Build phasor E0 from cos/sin amplitudes at same phase.

%   E(t) = [Ex_cos*cos(θ) + Ex_sin*sin(θ);
%           Ey_cos*cos(θ) + Ey_sin*sin(θ);
%           Ez_cos*cos(θ) + Ez_sin*sin(θ)]

    Ex0 = Ex_cos - 1j*Ex_sin;
    Ey0 = Ey_cos - 1j*Ey_sin;
    Ez0 = Ez_cos - 1j*Ez_sin;
    E0  = [Ex0; Ey0; Ez0];
    
end
