function [L, D] = Shortstub(ZL,Z0_load,Z0_stub)
syms l d real

lam = 1;
beta = 2*pi/lam;

Zi = Z0_load*(ZL + 1i*Z0_load*tan(beta*d)) / ...
     (Z0_load + 1i*ZL*tan(beta*d));
Yi = 1/Zi;

Ystub = -1i/Z0_stub * 1/tan(beta*l);

eq1 = imag(Yi + Ystub) == 0;
eq2 = real(Yi) == 1/Z0_load;

sol = vpasolve([eq1,eq2],[l,d],[0.1,0.2]);

L = double(sol.l);
D = double(sol.d);

fprintf("l_stub = %.4f λ\n", mod(L,0.5));
fprintf("l_tl   = %.4f λ\n", mod(D,0.5));
end

