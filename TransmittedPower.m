function T = TransmittedPower(mu1, eps1, mu2, eps2, theta_i_deg, theta_t_deg, incidence)

ti = deg2rad(theta_i_deg);
tt = deg2rad(theta_t_deg);

eta1 = sqrt(mu1/eps1);
eta2 = sqrt(mu2/eps2);

switch upper(incidence)

    case 'TM'
        GammaTM = (eta2*cos(tt) - eta1*cos(ti)) / ...
                  (eta2*cos(tt) + eta1*cos(ti));

        T = 1 - abs(GammaTM)^2;

    case 'TE'
        GammaTE = (eta2*cos(ti) - eta1*cos(tt)) / ...
                  (eta2*cos(ti) + eta1*cos(tt));

        T = 1 - abs(GammaTE)^2;

    otherwise
        error('Choose TE or TM');
end
end
