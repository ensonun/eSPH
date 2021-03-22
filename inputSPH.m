function [settings, f, dt_save, file_name] = inputSPH()

settings(1) = 2; %kernel function type
settings(2) = 0.15; %h
settings(3) = 7; %gamma for equation of state (usually =7)
settings(4) = 0; %reconstruction scheme (0 for linear)
settings(5) = 1; %Riemann solver type (0 for traditional, 1 for linear)
settings(6) = 3; %Riemann solver limiter parameter, eta
settings(7) = 2; %time integration order
settings(8) = 6; %simulation end time (non dimensionalised by sqrt(g/H))
settings(9) = 0.1; %max dt

td = 1;
%f = @(x,y,t) [0;-1]*( 0.5*(sin((t/td-0.5)*pi) + 1).*(t<=td) + (t>td) ) + 0*x';
%f = @(x,y,t) [0;0]*x';
f = @(x,y,t) [-x';-y'];

dt_save = 0.1;
file_name = 0; %use 0 for default
end