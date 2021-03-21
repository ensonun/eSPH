function p = pressure(rho,rho0,c0,gamma)
%PRESSURE Summary of this function goes here
%   Detailed explanation goes here
p = rho0.*c0.^2/gamma.*((rho./rho0).^gamma - 1);
end

