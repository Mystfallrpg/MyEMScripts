function isPW = PlaneWaveCheck(E0,H0,gamma)

% Find beta automatisk
beta = imag(gamma);

    if all(beta==0)
        beta = real(gamma);
    end

beta_hat = beta / norm(beta);

% Plane-wave tests
c1 = dot(beta,E0);
c2 = dot(beta,H0);

cp = cross(beta_hat,E0);

% Beregn effektiv eta
eta = norm(cp)/norm(H0);
c3 = norm(cp - eta*H0);

% Plane wave hvis alle ca 0
isPW = abs(c1)<1e-9 && abs(c2)<1e-9 && c3<1e-9;

    if isPW
        disp("Is a plane wave");
    else
        disp("Not a plane wave");
    end
end
