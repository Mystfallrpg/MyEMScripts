function Z = Impedance(plane, Z0, Z, l, Gamma, Gamma_plane)
% plane = 'A'  find ZA
% plane = 'B'  find ZB 
% Z = ZL hvis plane ='A'
% Z = ZA hvis plane ='B'
% Gamma = reflectionskoefficient
% Gamma_plane = GammaA || GammaB
% Lambda er givet som up/f

syms lambda
beta = 2*pi/lambda;

switch plane

    case 'A'
        if ~isempty(Gamma)

            if Gamma_plane == 'A'
            ZA = Z0 * (1 + Gamma) / (1 - Gamma);
            Z = vpa(ZA,8);
          
            elseif Gamma_plane == 'B'
            GammaA = Gamma .* exp(-1j * 2 * beta * l);
            ZA = Z0 * (1 + GammaA) / (1 - GammaA);
            Z = vpa(ZA,8);
            end

        else 
            if simplify(l/lambda) == 1/4
            % ZA from ZL
            ZA = Z0^2 / Z;
            Z = vpa(ZA,8);
            return
            end

            ZA = Z0 * (Z + 1j*Z0*tan(beta*l)) / (Z0 + 1j*Z*tan(beta*l));
            Z = vpa(ZA,8);

        end

    case 'B'

        if ~isempty(Gamma)

            if Gamma_plane == 'B'
            ZL = Z0 * (1 + Gamma) / (1 - Gamma);
            Z = vpa(ZL,8);

            elseif Gamma_plane == 'A'
            GammaL = Gamma .* exp(1j * 2 * beta * l);
            ZL = Z0 * (1 + GammaL) / (1 - GammaL);
            Z = vpa(ZL,8);
            end

        else
            if simplify(l/lambda) == 1/4
            % ZL from ZA
            ZL = Z0^2 / Z;
            Z = vpa(ZL,8);
            return
            end
            
            ZL = Z0 * (Z - 1j*Z0*tan(beta*l)) / (Z0 - 1j*Z*tan(beta*l));
            Z = vpa(ZL,8);
        end

end


end

