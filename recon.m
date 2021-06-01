function [F_L,F_R] = recon(rho,u,p,scheme)
%RECON Reconstruct left and right states of the local Riemann problem
% rho, u, p - arrays of the corresponding quantity at [h,i,j,k]

if scheme == 1 %MUSCL linear/2nd order: 4 points
    %b = 2; phi = @(r) max(0, min(b,r)); %Osher/minmod(b=1)
    %phi = @(r) (r+abs(r))./(1+abs(r)); %van Leer
    b = 2; phi = @(r) max([0*r; min(1,b*r); min(b,r)]); %Sweby/superbee(b=2)/minmod(b=1)
      
    dF_1 = [rho(2), u(2), p(2)] - [rho(1), u(1), p(1)];
    dF_2 = [rho(3), u(3), p(3)] - [rho(2), u(2), p(2)];
    dF_3 = [rho(4), u(4), p(4)] - [rho(3), u(3), p(3)];
    
    F_L = [rho(2), u(2), p(2)] + 0.5*phi(dF_2./(dF_1+0.001*(abs(dF_1)==0))).*dF_1;
    F_R = [rho(3), u(3), p(3)] - 0.5*phi(dF_2./(dF_3+0.001*(abs(dF_3)==0))).*dF_3;
    
elseif scheme == 2 %MUSCL parabolic/3rd order: 4 points
    phi = @(r) 2*r./(1+r.^2); %Albada
    %b = 2; phi = @(r) max([0*r; min(1,b*r); min(b,r)]); %Sweby/superbee(b=2)
    
    dF_1 = [rho(2), u(2), p(2)] - [rho(1), u(1), p(1)];
    dF_2 = [rho(3), u(3), p(3)] - [rho(2), u(2), p(2)];
    dF_3 = [rho(4), u(4), p(4)] - [rho(3), u(3), p(3)];
    
    F_L = [rho(2), u(2), p(2)] + 0.25*phi(dF_2./(dF_1+0.001*(abs(dF_1)==0))).*(2/3*dF_1 + 4/3*dF_2);
    F_R = [rho(3), u(3), p(3)] - 0.25*phi(dF_2./(dF_3+0.001*(abs(dF_3)==0))).*(2/3*dF_3 + 4/3*dF_2);

elseif scheme == 3 %3th order WENO-Z (Borges et al. 2007): 4 points
    %Left
    q0 = [rho(1), u(1), p(1)];
    q1 = [rho(2), u(2), p(2)];
    q2 = [rho(3), u(3), p(3)];
    
    IS0 = (q1-q0).^2;
    IS1 = (q2-q1).^2;
    tau = abs(IS0-IS1);
    a0 = (1+tau./(IS0+1e-6))/3;
    a1 = 2*(1+tau./(IS1+1e-6))/3;
    w0 = a0./(a0+a1);
    w1 = a1./(a0+a1);
    
    F_L = w0.*(-0.5*q0+1.5*q1) + w1.*(0.5*q1+0.5*q2);
    
    %Right
    q0 = [rho(4), u(4), p(4)];
    q1 = [rho(3), u(3), p(3)];
    q2 = [rho(2), u(2), p(2)];
   
    IS0 = (q1-q0).^2;
    IS1 = (q2-q1).^2;
    tau = abs(IS0-IS1);
    a0 = (1+tau./(IS0+1e-6))/3;
    a1 = 2*(1+tau./(IS1+1e-6))/3;
    w0 = a0./(a0+a1);
    w1 = a1./(a0+a1);
    
    F_R = w0.*(-0.5*q0+1.5*q1) + w1.*(0.5*q1+0.5*q2);
            
else %piecewise constant/1st order: 2 points
    F_L = [rho(2), u(2), p(2)];
    F_R = [rho(3), u(3), p(3)];
end

end