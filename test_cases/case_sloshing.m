% TEST CASE for linear sloshing
clear
clc

% INPUTS
F = [0 1 0 0.4]; %fluid domain =[xlim, ylim]
L = [0 1 0 0.5]; %wall domain =[xlim, ylim]
dx = 0.02; %particles spacing
lid = 0; % 1=with lid, 2=no lid

A = 0.25; %non-dimensional amplitude
g = 1;
freq = 0.3*sqrt(pi/F(2)*tanh(pi*F(4)/F(2)))/(2*pi) %excitation frequency
T = 3/freq %simulation end time = 3 period

settings(1) = 5; %kernel function type
settings(2) = 3*dx; %compact support radius, kh
settings(3) = 7; %gamma for equation of state (usually =7)
settings(4) = 0; %reconstruction scheme
settings(5) = 1; %Riemann solver type (0 for traditional, 1 for Roe)
settings(6) = 1; %second derivative, 0=off 1=on
settings(7) = T;%simulation end time (non dimensionalised by sqrt(g/H))
settings(8) = 1; %cfl#

dt_save = 0.1;
dir_name = "sl_"+num2str(settings(4))+"_"+num2str(dx*10^2)+"_"+num2str(settings(2)/dx);

%%
dw = ceil(settings(2)/dx)*dx; %wall thickness (=kh)
v_lid = 0;

f = @(x,y,t) [-A*(2*pi*freq)^2*sin(2*pi*freq*t);-g] + 0*x';

rho = 1; %density [unit/area]
v_max = sqrt(F(4)*g);

% Generate fluid particles
[xf,yf] = meshgrid(F(1)+dx/2:dx:F(2)-dx/2,F(3)+dx/2:dx:F(4)-dx/2);
xf = reshape(xf,[],1);
yf = reshape(yf,[],1);

N = length(xf) %number of particles
m = rho*dx*dx; %mass

fluid(:,1:2) = [xf,yf];
fluid(:,3) = rho*ones(N,1);
fluid(:,4) = m*ones(N,1);
fluid(:,5) = (F(4)-yf)*g; %pressure
fluid(:,6) = zeros(N,1); %x velocity
fluid(:,7) = zeros(N,1); %y velocity
fluid(:,9) = rho*ones(N,1); %rho_0
fluid(:,10) = 10*v_max*ones(N,1); %c_0
fluid(:,11) = zeros(N,1);

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
wall(:,3) = v_lid*(ywall>L(4)); %x velocity
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

% Calculate analytical solution
tz = 0:dt_save:T+dt_save;
xz = linspace(0,L(2))';
ele = zeros(length(xz),length(tz)); %interp1(t,y(:,1),tz);

a = @(t) -A*(2*pi*freq)^2*sin(2*pi*freq*t);
for m = 1:2:15
    w_m = sqrt(pi*m/L(2)*g*tanh(pi*m/L(2)*F(4)));
    dydt = @(t,y) [y(2); -w_m^2*y(1) - 2*a(t)*(2/pi/m)*tanh(pi*m/L(2)*F(4))];
    
    [t,fn] = ode45(dydt,[0 T+dt_save],[0 0]);
    
    ele = ele + cos(m*pi/L(2)*xz)*interp1(t,fn(:,1),tz);
end

figure(1)
subplot(1,2,1)
surf(xz.*ones(length(xz),length(tz)),tz.*ones(length(xz),length(tz)),ele,'EdgeColor','none')
xlabel('x')
ylabel('t')
zlabel('H-H_0')

LT = length(tz);
Ya = fft(ele(1,:));
P1 = abs(Ya/LT);
E_ana = P1(1:LT/2+1);
E_ana(2:end-1) = 2*E_ana(2:end-1);
fplot = 1/dt_save*(0:(LT/2))/LT;

subplot(1,2,2)
plot(fplot,E_ana)
hold on
plot([freq freq],[0 max(E_ana)])
xlabel('freq [s^{-1}]')
ylabel('Fourier transform')

%% RUN eSPH
eSPH(dir_name+".mat")

%% LOAD DATA
clearvars tn ke pe ie maxrho minrho meanrho x y u v p phase V
for n = 1:256
    load(['sl_0_2_3/SPHout_',num2str(n)],'fluid','t');
    %load(dir_name+"/SPHout_"+num2str(n),'fluid','wall','t');
    
    tn(n) = t;
    ke(n) = 0.5*sum( fluid(:,4).*(fluid(:,6).^2 + fluid(:,7).^2) );
    pe(n) = sum( fluid(:,4).*fluid(:,2) );
    ie(n) = sum( fluid(:,4).*fluid(:,10).^2/7.*((fluid(:,3)./fluid(:,9)).^6/6 ...
        + fluid(:,9)./fluid(:,3) -1-1/6) );
    
    %maxrho(n) = max(fluid(:,3));
    %minrho(n) = min(fluid(:,3));
    %meanrho(n) = mean(fluid(:,3));
    
    x(:,n) = fluid(:,1);
    y(:,n) = fluid(:,2);
    p(:,n) = fluid(:,5);    
end

%% ELEVATION vs TIME
dx = 0.02;
z0 = 0*tn;
z1 = 0*tn;
for tt = 1:length(tn)
    leftsurf = find(y(:,tt)==max(y(x(:,tt)<1*dx,tt)));
    z0(tt) = max(y(x(:,tt)<1.5*dx,tt)) + dx/2;
    
    rightsurf = find(y(:,tt)==max(y(x(:,tt)>(1-1*dx),tt)));
    z1(tt) = max(y(x(:,tt)>(1-1.5*dx),tt)) + dx/2;    
end

figure(2)
subplot(2,1,1)
plot(tz,ele(1,:),'--',tz,ele(end,:),'--')
hold on
plot(tn,z0-F(4))
plot(tn,z1-F(4))

err0 = abs(z0-F(4)-interp1(tz,ele(1,:),tn));
err1 = abs(z1-F(4)-interp1(tz,ele(end,:),tn));
L2 = sqrt(sum([err0,err1].^2)/length([err0,err1]))
subplot(2,1,2)
hold on
plot(tn,[err0;err1])

%% FOURIER SPECTRUM
Yn = fft(z0-F(4));

L = length(tn);
fs = 1/0.1; %1/dt_save

P2 = abs(Yn/L);
E_num = P2(1:L/2+1);
E_num(2:end-1) = 2*E_num(2:end-1);

fplot2 = fs*(0:(L/2))/L;

figure(3)
plot(fplot,E_ana)
hold on
plot(fplot2,E_num,'--')

%% ANIMATION
for tt = 1:length(tn)
    figure(4)
    scatter(x(:,tt),y(:,tt),[],p(:,tt),'filled')
    hold on
    scatter(wall(:,1),wall(:,2),'x')
    colorbar
    caxis([0 F(4)*g])
    
    plot(xz,ele(:,tt+1)+F(4),'r','LineWidth',2)
    title(['t = ',num2str(tn(tt))])
    axis equal
    hold off
    pause(0.02)
end