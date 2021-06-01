function rho = density(p,rho0,c0,gamma)
%DENSITY Computes density from pressure via EoS

rho = rho0.*(gamma*p./(rho0.*c0.^2)+1).^(1/gamma);
end

