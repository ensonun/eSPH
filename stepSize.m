function dt = stepSize(c0,v,h,f,nu,maxdt)
%STEPSIZE Summary of this function goes here
%   Detailed explanation goes here

% CFL condition
dt(1) = 0.25*h/(c0+v);

% Viscous condition
dt(2) = 0.25*h^2/nu;

% Body force condition
dt(3) = 0.25*sqrt(h/f);

dt(4) = maxdt;

% global timestep
dt = min(dt);
end

