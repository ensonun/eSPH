% TEST CASE for oscillating droplet in a conservative force field
clear
clc

% INPUTS
dx = 0.02; %particles spacing
A0 = 1.5; %velocity scale of the case

settings(1) = 5; %kernel function type
settings(2) = 3*dx; %compact support radius, kh
settings(3) = 7; %gamma for equation of state (usually =7)
settings(4) = 3; %reconstruction scheme
settings(5) = 1; %Riemann solver type (0 for traditional, 1 for Roe)
settings(6) = 1; %second derivative, 0=off 1=on
settings(7) = 8;%simulation end time (non dimensionalised by sqrt(g/H))
settings(8) = 1; %cfl#

dt_save = 0.1;
dir_name = "od_"+num2str(settings(4))+"_"+num2str(dx*10^2)+"_"+num2str(settings(2)/dx);

%% GENERATE INPUT FILE
f = @(x,y,t) A0^2*[-x';-y'];

rho = 1; %density [unit/area]
mu = 0;
v_max = 1.5*A0;

% Generate fluid particles
F = [-1 1 -1 1]; %fluid domain =[xlim, ylim]
[xf,yf] = meshgrid(F(1)+dx/2:dx:F(2)-dx/2,F(3)+dx/2:dx:F(4)-dx/2);
in = sqrt(xf.^2+yf.^2)<=1;
xf = reshape(xf(in),[],1); 
yf = reshape(yf(in),[],1);

N = length(xf) %number of particles
m = rho*dx*dx; %pi*rho/N

fluid(:,1:2) = [xf,yf];
fluid(:,3) = rho*ones(N,1);
fluid(:,4) = m*ones(N,1);
fluid(:,5) = 0.5*rho*A0^2*(1-xf.^2-yf.^2); %pressure
fluid(:,6) = A0*xf; %x velocity
fluid(:,7) = -A0*yf; %y velocity
fluid(:,9) = rho*ones(N,1); %rho0
fluid(:,10) = 10*v_max*ones(N,1); %c0
fluid(:,11) = mu*ones(N,1);

figure(1)
hold on
scatter(xf,yf,600,fluid(:,5),'.')
quiver(xf,yf,fluid(:,6),fluid(:,7))
axis equal
colorbar

% Generate wall particles
wall = ones(0,4);

save(dir_name+".mat",'wall','fluid','f','settings','dt_save','dir_name')

% Calculate analytical sol
odefun = @(t,x) [(x(1)^2+A0^2)*(1-x(2)^4)/(1+x(2)^4); x(1)*x(2)];
tspan = [0 settings(7)];
x0 = [A0; 1];
[ta,xa] = ode23(odefun,tspan,x0);

figure(1) % semi-axis vs t
plot(ta, xa(:,2))
hold on
ylabel('a(t)')

figure(2) %ke vs t
plot(ta, pi/8*xa(:,1).^2.*(xa(:,2).^2 + xa(:,2).^(-2))*rho)
hold on
ylabel('KE(t)')

%% RUN eSPH
eSPH(dir_name+".mat")

%% LOAD DATA
clearvars tn ke pe ie maxrho minrho meanrho x y u v p phase
for n = 1:80
    load(dir_name+"/SPHout_"+num2str(n),'fluid','wall','t');
    
    tn(n) = t;
    ke(n) = 0.5*sum( fluid(:,4).*(fluid(:,6).^2 + fluid(:,7).^2) );
    pe(n) = sum( fluid(:,4).*fluid(:,2) );
    ie(n) = sum( fluid(:,4).*fluid(:,10).^2/6.*((fluid(:,3)./fluid(:,9)).^6 ...
        + (fluid(:,9)./fluid(:,3)).*6 - 7) );
    
    x(:,n) = fluid(:,1);
    y(:,n) = fluid(:,2);
    u(:,n) = fluid(:,6);
    v(:,n) = fluid(:,7);
    p(:,n) = fluid(:,5);
end


figure(1)
hold on
plot(tn, 0.25*(max(x)-min(x))+1./(max(y)-min(y)),'o')

figure(2) %ke, pe vs t
yyaxis('left')
hold on
plot(tn, ke)

yyaxis('right')
hold on
plot(tn, pe)

%% ANIMATION
theta = linspace(-pi,pi);
for tt = 1:length(tn)
    figure(3)
    scatter(x(:,tt),y(:,tt),200,p(:,tt),'.')
    hold on    
    plot(cos(theta),sin(theta),'k')
    plot(1.9314*cos(theta),0.5179*sin(theta),'k--')
    plot(0.5179*cos(theta),1.9314*sin(theta),'k--')
    
    colorbar
    caxis([0 1.125])
    title(['t = ',num2str(tn(tt))])
    axis equal
    hold off    
    pause(0.1)
end