function dt = stepSize(c0,v,h,f,nu,C)
%STEPSIZE Compute timestep size

% CFL condition
dt(1) = 0.25*C*h/(c0+v);

% Viscous condition
dt(2) = 0.25*h^2/nu;

% Body force condition
dt(3) = 0.25*sqrt(h/f);

% global timestep
dt = min(dt);
end

