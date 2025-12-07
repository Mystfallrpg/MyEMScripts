function [Gamma, VSWR] = Reflection(plane,Z0,Z,l,Gamma)

syms lambda
beta = 2*pi/lambda;

switch plane
    case 'A'
        if ~isempty(Gamma)

            GammaA = Gamma .* exp(-1j * 2 * beta * l);
            Gamma = vpa(GammaA,8);
        else
            GammaA = (Z - Z0)/(Z + Z0);
            Gamma = vpa(GammaA,8);
        end


    case 'B'
        if ~isempty(Gamma)

            GammaB = Gamma .* exp(1j * 2 * beta * l);
            Gamma = vpa(GammaB,8);
        else
            GammaB = (Z - Z0)/(Z + Z0);
            Gamma = vpa(GammaB,8);
        end

end

VSWR = (1+abs(Gamma))/(1-abs(Gamma));

end