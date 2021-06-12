% TEST CASE for dam-breaking
clear
clc

% INPUTS
dx = 0.02; %particles spacing
lid = 0; % 1=with lid, 2=no lid

settings(1) = 5; %kernel function type
settings(2) = 3*dx; %compact support radius, kh
settings(3) = 7; %gamma for equation of state (usually =7)
settings(4) = 0; %reconstruction scheme
settings(5) = 1; %Riemann solver type (0 for traditional, 1 for Roe)
settings(6) = 1; %second derivative, 0=off 1=on
settings(7) = 10;%simulation end time (non dimensionalised by sqrt(g/H))
settings(8) = 1; %cfl#

dt_save = 0.1;
dir_name = "db_"+num2str(settings(4))+"_"+num2str(dx*10^2)+"_"+num2str(settings(2)/dx);

%% GENERATE INPUT FILE
F = [0 2 0 1]; %fluid domain =[xlim, ylim]
f = @(x,y,t) [0;-1] + 0*x';

L = [0 floor(5.366/dx)*dx 0 3]; %wall domain =[xlim, ylim]
dw = settings(2); %wall thickness (=h)
rho = 1; %density [unit/area]
mu = 0;
v_max = 2*sqrt(F(4));

% Generate fluid particles
[xf,yf] = meshgrid(F(1)+dx/2:dx:F(2)-dx/2,F(3)+dx/2:dx:F(4)-dx/2);
xf = reshape(xf,[],1);
yf = reshape(yf,[],1);

N = length(xf) %number of particles
m = rho*dx*dx; %mass

fluid(:,1:2) = [xf,yf];
fluid(:,3) = rho*ones(N,1);
fluid(:,4) = m*ones(N,1);
fluid(:,5) = F(4)-yf; %pressure
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
wall(:,3) = zeros(N+N_lid,1); %x velocity
wall(:,4) = zeros(N+N_lid,1); %y velocity

% figure(1)
% scatter(xf,yf,100,fluid(:,5),'.')
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
clearvars tn ke pe ie maxrho minrho meanrho rho m x y u v p phase
for n = 1:100
    load(dir_name+"/SPHout_"+num2str(n),'fluid','t');
    
    tn(n) = t;
    x(:,n) = fluid(:,1);
    y(:,n) = fluid(:,2);
    p(:,n) = fluid(:,5);
    
    indomain = (y(:,n)>0)&(x(:,n)<5.366);
    ke(n) = 0.5*sum( fluid(indomain,4).*(fluid(indomain,6).^2 + fluid(indomain,7).^2) );
    pe(n) = sum( fluid(indomain,4).*fluid(indomain,2) );
    ie(n) = sum( fluid(:,4).*fluid(:,10).^2/6.*((fluid(:,3)./fluid(:,9)).^6 ...
        + (fluid(:,9)./fluid(:,3)).*6 - 7) );
    
    %maxrho(n) = max(fluid(:,3));
    %minrho(n) = min(fluid(:,3));
    %meanrho(n) = mean(fluid(:,3));
    
    rho(:,n) = fluid(:,3);
    m(:,n) = fluid(:,4);
    u(:,n) = fluid(:,6);
    v(:,n) = fluid(:,7);
end

figure(2)
hold on
plot(tn,ke+pe+ie)

% Plot wavefront
dx = 0.02;
figure(3)
plot(tn(1:25),max(x(:,1:25))+dx/2)
hold on
plot([1 1.5 2.5],[0.5 1.5 3.5]+2,'k--')

%% ANIMATION
for tt = 1:length(tn)
    figure(6)
    scatter(x(:,tt),y(:,tt),50,p(:,tt),'.')
    hold on
    scatter(wall(:,1),wall(:,2),'x')
    colorbar
    caxis([0 1])
    
    title(['t = ',num2str(tn(tt))])
    axis equal
    hold off
    pause(0.1)
end

%% PRESSURE PROBE
ploc = [5.36 0.2]; %location of pressure probe
T_len = 80;

pp = zeros(1,T_len);
exp = [1.875594978	0.014385121;    2.094612445	0.00812887;    2.253029476	0.011256996;    2.287838428	-0.025642117;
    2.528215066	0.015023514;    2.616975983	0.362841274;    2.649530568	0.598067797;    2.687740175	0.65122465;
    2.757358079	0.691273168;    2.855098253	0.662502926;    2.879819869	0.598684911;    2.998881004	0.598684911;
    3.128067686	0.569276275;    3.263291485	0.546060051;    3.532860262	0.536675675;    3.672134279	0.536675675;
    3.935016376	0.552316302;    4.205731441	0.555444428;    4.365218341	0.555444428;    4.533722707	0.54980529;
    4.694355895	0.525418679;    4.901031659	0.517906922;    5.181872271	0.512267785;    5.353700873	0.53354755;
    5.532330786	0.58608729;    5.678367904	0.614240419;    5.746877729	0.713084928;    5.817641921	0.773774818;
    5.849088428	0.840082565;    5.873810044	0.891388079;    5.959170306	0.726852936;    6.028826419	0.58421467;
    6.117549127	0.556061541;    6.236610262	0.561062286;    6.339967249	0.63677569;    6.399497817	0.537931181;
    6.491621179	0.474730279;    6.59944869	0.489753793;    6.725272926	0.51539591;    6.791528384	0.55793416;
    6.827483624	0.477241291;    6.894885371	0.410933543;    7.024072052	0.430957802;    7.113902838	0.501649182;
    7.218367904	0.44470453;    7.376784934	0.442831911];
kh = settings(2);
w = @(r) (1+4*r/kh).*(2-2*r/kh).^4*(7/16/pi/kh^2).*(r<=kh);

for tt = 1:T_len
    idj = rangesearch([x(:,tt),y(:,tt)],ploc,kh);
    sumW = 0;
    for i = idj{1}
        r = ploc - [x(i,tt),y(i,tt)];
        pp(tt) = pp(tt) + 2*m(i,tt)*p(i,tt)*w(norm(r))/rho(i,tt);
    end
end

figure(5)
plot(tn(1:T_len),pp,'LineWidth',1.5)
hold on
plot(exp(:,1),exp(:,2),'k^')