function [F_L,F_R] = recon(rho,u,p,scheme)
%RECON Summary of this function goes here
%   Detailed explanation goes here

if scheme == 1 %MUSCL: 4 points
    %b = 2; phi = @(r) max(0, min(b,r)); %Osher/minmod(b=1)
    %phi = @(r) (r+abs(r))./(1+abs(r)); %van Leer
    b = 2; phi = @(r) max([0*r; min(1,b*r); min(b,r)]); %Sweby/superbee(b=2)
    
    dF_1 = [rho(2), u(2), p(2)] - [rho(1), u(1), p(1)];
    dF_2 = [rho(3), u(3), p(3)] - [rho(2), u(2), p(2)];
    dF_3 = [rho(4), u(4), p(4)] - [rho(3), u(3), p(3)];
    
    F_L = [rho(2), u(2), p(2)] + 0.5*phi(dF_2./(dF_1+0.001*(abs(dF_1)==0))).*dF_1;
    F_R = [rho(3), u(3), p(3)] - 0.5*phi(dF_2./(dF_3+0.001*(abs(dF_3)==0))).*dF_3;
    
elseif scheme == 2 %modified WENO (Zhang 2019): 4 points
    
else %linear/1st order: 2 points
    F_L = [rho(2), u(2), p(2)];
    F_R = [rho(3), u(3), p(3)];
end

end