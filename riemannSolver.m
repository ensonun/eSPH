function [us, ps] = riemannSolver(F_L,F_R,c0,gh,type,eta)
%RIEMANNSOLVER Summary of this function goes here
%   Detailed explanation goes here

if type == 0
    us = 0.5*(F_L(2)+F_R(2));   
    ps = (F_R(1)*F_L(3)+F_L(1)*F_R(3))/(F_L(1)+F_R(1));
    
elseif type == 1 %Roe solver with dissipation limters
    alpha = (abs(F_L(3)-F_R(3)) > gh);  %no diss when dp <= g*h 
    us = (F_L(1)*c0(1)*F_L(2) + F_R(1)*c0(2)*F_R(2) + alpha*(F_L(3)-F_R(3)))/(F_L(1)*c0(1) + F_R(1)*c0(2));
    
    c_bar = (F_L(1)*c0(1)+F_R(1)*c0(2))/(F_L(1)+F_R(1));   
    beta = min( 1, eta*max( (F_L(2)-F_R(2))/c_bar, 0 ) ); %no diss across rarefraction
    ps = (F_L(1)*c0(1)*F_R(3) + F_R(1)*c0(2)*F_L(3) + F_L(1)*c0(1)*F_R(1)*c0(2)*beta*(F_L(2)-F_R(2)))/(F_L(1)*c0(1) + F_R(1)*c0(2));

    % rho_bar = 0.5*(F_L(1)+F_R(1));
    % beta = min( c0, eta*max(F_L(2)-F_R(2),0) );
    % us = 0.5*(F_L(2)+F_R(2) + (F_L(3)-F_R(3))/rho_bar/c0);
    % ps = 0.5*(F_L(3)+F_R(3) + beta*rho_bar*(F_L(2)-F_R(2)));
end

end

