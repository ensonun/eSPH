function [drhodt, dvxdt, dvydt] = timeDer(fluid,wall,settings,f,n_w)
%TIMEDER Summary of this function goes here
%   Detailed explanation goes here

N = size(fluid,1); %number of fluid particles (= number of d?/dt)
x = fluid(:,1);
y = fluid(:,2);
rho = fluid(:,3);
m = fluid(:,4);
p = fluid(:,5);
vx = fluid(:,6);
vy = fluid(:,7);
%phase = fluid(:,8);
rho0 = fluid(:,9);
c0 = fluid(:,10);
nu = fluid(:,11);

w_type = settings(1);
h = settings(2);
gamma = settings(3);

% Range search for fluid particles
j_id = rangesearch([x,y],[x,y],h);

if size(wall,2) > 0
    xw = wall(:,1);
    yw = wall(:,2);
    %mw = wall(:,3);
    vw_x = wall(:,4);
    vw_y = wall(:,5);
    
    % Range search for wall particles
    k_id = rangesearch([xw,yw],[x,y],h);
else
    xw = 0; %must be here for parfor
    yw = 0;
    vw_x = 0;
    vw_y = 0;
    
    k_id = cell(N,1); %empty cell
end

% Initialise variables
drhodtf = zeros(N,1);
dvxdtf = zeros(N,1);
dvydtf = zeros(N,1);

drhodtw = zeros(N,1);
dvxdtw = zeros(N,1);
dvydtw = zeros(N,1);

% Prepare for extrapolation
gradrho = zeros(N,2);
gradvx  = zeros(N,2);
gradvy  = zeros(N,2);
gradp   = zeros(N,2);
if settings(4) ~= 0       
    %par
    for i = 1:N        
        for j = j_id{i}(2:end)
            rij_x = x(i) - x(j);
            rij_y = y(i) - y(j);            
            norm_rij = sqrt(rij_x^2 + rij_y^2);        
            eij = [rij_x, rij_y]/max( norm_rij, 0.001*h );
            dwdr = wDer(norm_rij,h,w_type);
            
            gradrho(i,:) = gradrho(i,:) + m(j)/rho(j)*(rho(j)-rho(i))*dwdr*eij;
            gradvx(i,:)  = gradvx(i,:)  + m(j)/rho(j)*(vx(j)-vx(i))*dwdr*eij;
            gradvy(i,:)  = gradvy(i,:)  + m(j)/rho(j)*(vy(j)-vy(i))*dwdr*eij;
            gradp(i,:)   = gradp(i,:)   + m(j)/rho(j)*(p(j)-p(i))*dwdr*eij;
        end
    end    
end


%par
for i = 1:N %parallel loop through each fluid particle
    
    for j = j_id{i}(2:end) %loop through neighbour fluid particles
        % Compute distance vectors
        rij_x = x(i) - x(j);
        rij_y = y(i) - y(j);
        norm_rij = sqrt(rij_x^2 + rij_y^2);
        eij_x = rij_x/max( norm_rij, 0.001*h );
        eij_y = rij_y/max( norm_rij, 0.001*h );
        dwdr = wDer(norm_rij,h,w_type);
        
        % Construct Riemann problem: 4-point stencil h--i--j--k
        rhoh = rho(i) - gradrho(i,:)*[rij_x;rij_y];
        vxh = vx(i) - gradvx(i,:)*[rij_x;rij_y];
        vyh = vy(i) - gradvy(i,:)*[rij_x;rij_y];
        ph = p(i) - gradp(i,:)*[rij_x;rij_y];
        
        rhok = rho(j) + gradrho(j,:)*[rij_x;rij_y];
        vxk = vx(j) + gradvx(j,:)*[rij_x;rij_y];
        vyk = vy(j) + gradvy(j,:)*[rij_x;rij_y];
        pk = p(j) + gradp(j,:)*[rij_x;rij_y];
        
        u = -[vxh,vx(i),vx(j),vxk]*eij_x -[vyh,vy(i),vy(j),vyk]*eij_y;
        [F_L,F_R] = recon([rhoh,rho(i),rho(j),rhok],u,[ph,p(i),p(j),pk],settings(4));
        
        % Solve Riemann problem
        [us, ps] = riemannSolver(F_L,F_R,[c0(i) c0(j)],norm(f(:,i))*h,settings(5),settings(6));
        vsx = -(us - 0.5*(F_L(2)+F_R(2)))*eij_x + 0.5*(vx(i)+vx(j));
        vsy = -(us - 0.5*(F_L(2)+F_R(2)))*eij_y + 0.5*(vy(i)+vy(j));
        
        % Viscosity
        if nu(i)==0 || nu(j)==0
            visc = 0;
        else
            visc = -2*(nu(i)+nu(j))*((vx(i)-vx(j))*eij_x+(vy(i)-vy(j))*eij_x)/norm_rij/rho(i)/rho(j);
        end
        
        % Evaluate fluxes
        drhodtf(i) = drhodtf(i) + 2*rho(i)*m(j)/rho(j)*((vx(i)-vsx)*eij_x+(vy(i)-vsy)*eij_y)*dwdr;
        dvxdtf(i) = dvxdtf(i) - 2*m(j)/rho(i)/rho(j)*(ps+visc)*dwdr*eij_x;
        dvydtf(i) = dvydtf(i) - 2*m(j)/rho(i)/rho(j)*(ps+visc)*dwdr*eij_y;
    end
    
    for k = k_id{i} %loop through neighbour wall particles
        % Compute distance vectors
        rij_x = x(i) - xw(k);
        rij_y = y(i) - yw(k);
        norm_rij = sqrt(rij_x^2 + rij_y^2);
        eij_x = rij_x/max( norm_rij, 0.001*h );
        eij_y = rij_y/max( norm_rij, 0.001*h );
        dwdr = wDer(norm_rij,h,w_type);
        
        % Velocity
        nw_x = n_w(k,1);
        nw_y = n_w(k,2);
        u = [0, -vx(i)*nw_x-vy(i)*nw_y, ...
              vx(i)*nw_x+vy(i)*nw_y +2*sqrt(vw_x(k)^2+vw_y(k)^2), 0];
        % Pressure
        pw = p(i) + rho(i)*max( 0, -(f(1,i)*rij_x+f(2,i)*rij_y) );        
        % Wall density
        rhow = density(pw,rho0(i),c0(i),gamma);
        
        % Reconstruction
        [F_L,F_R] = recon([0,rho(i),rhow,0],u,[0,p(i),pw,0],0);
        
        % Solve Riemann problem
        [us, ps] = riemannSolver(F_L,F_R,[c0(i) c0(i)],norm(f(:,i))*h,settings(5),settings(6));
        vsx = -(us - 0.5*(F_L(2)+F_R(2)))*eij_x + 0.5*(vx(i)+2*vw_x(k)-vx(i));
        vsy = -(us - 0.5*(F_L(2)+F_R(2)))*eij_y + 0.5*(vy(i)+2*vw_y(k)-vy(i));
        
        % Viscosity
        if nu(i)==0
            visc = 0;
        else
            visc = -4*nu(i)*((vx(i)-vw_x(k))*eij_x+(vy(i)-vw_y(k))*eij_x)/norm_rij/rho(i)/rhow;
        end
        
        % Evaluate fluxes
        drhodtw(i) = drhodtw(i) + 2*rho(i)*m(i)/rhow*((vx(i)-vsx)*eij_x+(vy(i)-vsy)*eij_y)*dwdr;
        dvxdtw(i) = dvxdtw(i) - 2*m(i)/rho(i)/rhow*(ps+visc)*dwdr*eij_x;
        dvydtw(i) = dvydtw(i) - 2*m(i)/rho(i)/rhow*(ps+visc)*dwdr*eij_y;
    end
end

drhodt = drhodtf + drhodtw;
dvxdt = dvxdtf + dvxdtw + f(1,:)';
dvydt = dvydtf + dvydtw + f(2,:)';

end