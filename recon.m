function [F_L,F_R] = recon(rho,u,p,scheme)
%RECON Summary of this function goes here
%   Detailed explanation goes here

if scheme == 1 %MUSCL/2nd order: 4 points
    %b = 1; phi = @(r) max(0, min(b,r)); %Osher/minmod(b=1)
    %phi = @(r) (r+abs(r))./(1+abs(r)); %van Leer
    b = 2; phi = @(r) max([0*r; min(1,b*r); min(b,r)]); %Sweby/superbee(b=2)
    
    dF_1 = [rho(2), u(2), p(2)] - [rho(1), u(1), p(1)];
    dF_2 = [rho(3), u(3), p(3)] - [rho(2), u(2), p(2)];
    dF_3 = [rho(4), u(4), p(4)] - [rho(3), u(3), p(3)];
    
    F_L = [rho(2), u(2), p(2)] + 0.5*phi(dF_2./(dF_1+0.001*(abs(dF_1)==0))).*dF_1;
    F_R = [rho(3), u(3), p(3)] - 0.5*phi(dF_2./(dF_3+0.001*(abs(dF_3)==0))).*dF_3;
        
elseif scheme == 2 %Parabolic/3rd order: 4 points
    phi = @(r) 2*r./(1+r.^2); %Albada
    
    dF_1 = [rho(2), u(2), p(2)] - [rho(1), u(1), p(1)];
    dF_2 = [rho(3), u(3), p(3)] - [rho(2), u(2), p(2)];
    dF_3 = [rho(4), u(4), p(4)] - [rho(3), u(3), p(3)];
    
    F_L = [rho(2), u(2), p(2)] + 0.25*phi(dF_2./(dF_1+0.001*(abs(dF_1)==0))).*(2/3*dF_1 + 4/3*dF_2);
    F_R = [rho(3), u(3), p(3)] - 0.25*phi(dF_2./(dF_3+0.001*(abs(dF_3)==0))).*(2/3*dF_3 + 4/3*dF_2);
    
elseif scheme == 6 %modified 6th order WENO (Zhang 2019): 4 points
    
    
else %piecewise constant/1st order: 2 points
    F_L = [rho(2), u(2), p(2)];
    F_R = [rho(3), u(3), p(3)];
end

end