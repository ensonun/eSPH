function rho = density(p,rho0,c0,gamma)
%DENSITY Summary of this function goes here
%   Detailed explanation goes here
rho = rho0.*(gamma*p./(rho0.*c0.^2)+1).^(1/gamma);
end

