% TEST CASE for hydrostatic tank
clear
clc

% INPUTS
F = [0 1 0 0.5]; %fluid domain =[xlim, ylim]
L = [0 1 0 1]; %wall domain =[xlim, ylim]
dx = 0.02; %particles spacing
td = 1; %daming time
lid = 0; % 1=with lid, 2=no lid

settings(1) = 5; %kernel function type
settings(2) = 2*dx; %compact support radius, kh
settings(3) = 7; %gamma for equation of state (usually =7)
settings(4) = 0; %reconstruction scheme
settings(5) = 1; %Riemann solver type (0 for traditional, 1 for Roe)
settings(6) = 1; %second derivative, 0=off 1=on
settings(7) = 10;%simulation end time (non dimensionalised by sqrt(g/H))
settings(8) = 1; %cfl#

dt_save = 0.1;
dir_name = "hs_"+num2str(settings(4))+"_"+num2str(dx*10^2)+"_"+num2str(settings(2)/dx);

%% GENERATE INPUT FILE
f = @(x,y,t) [0;-1]*( 0.5*(sin((t/td-0.5)*pi) + 1).*(t<=td) + (t>td) ) + 0*x';
dw = settings(2); %wall thickness (=kh)

mu = 0;
rho = 1; %density [unit/area]
v_lid = 0;
v_max = sqrt(F(4));

% Generate fluid particles
[xf,yf] = meshgrid(F(1)+dx/2:dx:F(2)-dx/2, F(3)+dx/2:dx:F(4)-dx/2);
xf = reshape(xf,[],1); 
yf = reshape(yf,[],1);

N = length(xf) %number of particles
m = rho*dx*dx; %mass

fluid(:,1:2) = [xf,yf];
fluid(:,3) = rho*ones(N,1);
fluid(:,4) = m*ones(N,1);
fluid(:,5) = zeros(N,1); %pressure
fluid(:,6) = zeros(N,1); %x velocity
fluid(:,7) = zeros(N,1); %y velocity
fluid(:,9) = rho*ones(N,1); %rho0
fluid(:,10) = 10*v_max*ones(N,1); %c0
fluid(:,11) = mu*ones(N,1);

% Generate wall particles
% Side walls
[xside,yside] = meshgrid([L(1)-dw+dx/2:dx:L(1)-dx/2,L(2)+dx/2:dx:L(2)+dw-dx/2],L(3)-dw+dx/2:dx:L(4)+dw-dx/2);
xside = reshape(xside,[],1); yside = reshape(yside,[],1);
% Horizontal walls
[xhor,yhor] = meshgrid(L(1)+dx/2:dx:L(2)-dx/2,L(3)-dw+dx/2:dx:L(3)-dx/2);
xhor = reshape(xhor,[],1); yhor = reshape(yhor,[],1);
xwall = [xside;xhor];
ywall = [yside;yhor];
N = length(xwall);

if lid == 1
    [xlid,ylid] = meshgrid(L(1)+dx/2:dx:L(2)-dx/2,L(4)+dx/2:dx:L(4)+dw-dx/2);
    xlid = reshape(xlid,[],1); ylid = reshape(ylid,[],1);
    xwall = [xwall;xlid];
    ywall = [ywall;ylid];
    N_lid = length(xlid);
else
    N_lid = 0;
end

wall(:,1:2) = [xwall,ywall];
wall(:,3) = v_lid*(ywall>L(4)&xwall<L(2)); %x velocity
wall(:,4) = zeros(N+N_lid,1); %y velocity

% figure(1)
% scatter(xf,yf,600,fluid(:,5),'.')
% hold on
% quiver(xf,yf,fluid(:,6),fluid(:,7))
% plot(xwall,ywall,'kx')
% quiver(xwall,ywall,wall(:,3),wall(:,4))
% axis equal
% colorbar

save(dir_name+".mat",'wall','fluid','f','settings','dt_save','dir_name')
  
%% RUN eSPH
eSPH(dir_name+".mat")

%% LOAD DATA
clearvars tn ke pe ie maxrho minrho meanrho x y u v p phase
for n = 1:100
    load(dir_name+"/SPHout_"+num2str(n),'fluid','wall','t');
    
    tn(n) = t;
    ke(n) = 0.5*sum( fluid(:,4).*(fluid(:,6).^2 + fluid(:,7).^2) );
    pe(n) = sum( fluid(:,4).*fluid(:,2) );
    ie(n) = sum( fluid(:,4).*fluid(:,10).^2/6.*((fluid(:,3)./fluid(:,9)).^6 ...
        + (fluid(:,9)./fluid(:,3)).*6 - 7) );
    
    maxrho(n) = max(fluid(:,3));
    minrho(n) = min(fluid(:,3));
    meanrho(n) = mean(fluid(:,3));    
    
    x(:,n) = fluid(:,1);
    y(:,n) = fluid(:,2);
    u(:,n) = fluid(:,6);
    v(:,n) = fluid(:,7);
    p(:,n) = fluid(:,5);    
end

figure(2)
yyaxis('left')
plot(tn,ke+pe+ie)

yyaxis('right')
hold on
plot(tn,meanrho)
plot(tn,maxrho)

%% ANIMATION
for tt = 1:length(tn)
    figure(3)
    scatter(x(:,tt),y(:,tt),200,p(:,tt),'.')
    %scatter(x(:,tt),y(:,tt),200,p(:,tt)-0.5+y(:,tt),'.')
    
    hold on    
    scatter(wall(:,1),wall(:,2),'x')    
    colorbar
    
    title(['t = ',num2str(tn(tt))])
    axis equal
    hold off    
    pause(0.1)
end

%% PLOT PRESSURE vs Y
H = 0.5; % F(4);
for tt = 1:length(tn)
    figure(3)
    plot([0 H],[H 0])
    hold on    
    plot(y(:,tt),p(:,tt),'.')
    
    title(['t = ',num2str(tn(tt))])
    axis equal
    hold off    
    pause(0.1)
end
