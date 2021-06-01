function p = pressure(rho,rho0,c0,gamma)
%PRESSURE Computes pressure from pressure via EoS

p = rho0.*c0.^2/gamma.*((rho./rho0).^gamma - 1);
end

