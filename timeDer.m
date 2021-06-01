function [drhodt, dvxdt, dvydt] = timeDer(fluid,wall,settings,f,n_w)
%TIMEDER Computes time derivatives [drhodt, dvxdt, dvydt]

N = size(fluid,1); %number of fluid particles (= number of d?/dt)
x = fluid(:,1);
y = fluid(:,2);
rho = fluid(:,3);
m = fluid(:,4);
p = fluid(:,5);
vx = fluid(:,6);
vy = fluid(:,7);
rho0 = fluid(:,9);
c0 = fluid(:,10);
nu = fluid(:,11);

w_type = settings(1);
kh = settings(2);
gamma = settings(3);
recon_order = settings(4);

% Range search for fluid particles
j_id = rangesearch([x,y],[x,y],kh);

if size(wall,2) > 0
    xw = wall(:,1);
    yw = wall(:,2);
    vw_x = wall(:,3);
    vw_y = wall(:,4);
    
    % Range search for wall particles
    k_id = rangesearch([xw,yw],[x,y],kh);
else
    %must be here for parfor
    xw = 0;
    yw = 0;
    vw_x = 0;
    vw_y = 0;
    
    k_id = cell(N,1); %empty cells
end

% Initialise variables
drhodt = zeros(N,1);
dvxdt = zeros(N,1);
dvydt = zeros(N,1);

% Prepare for extrapolation
gradrho = zeros(N,2); %[drho/dx,drho/dy]
gradvx  = zeros(N,2);
gradvy  = zeros(N,2);
gradp   = zeros(N,2);
if settings(4) > 0
    parfor i = 1:N
        B = [0 0;0 0];
        for j = j_id{i}(2:end)
            rij_x = x(i) - x(j);
            rij_y = y(i) - y(j);
            norm_rij = sqrt(rij_x^2 + rij_y^2);
            eij = [rij_x, rij_y]/max( norm_rij, 0.001*kh );
            dwdr = wDer(norm_rij,kh,w_type);
            
            gradrho(i,:) = gradrho(i,:) + m(j)/rho(j)*(rho(j)-rho(i))*dwdr*eij;
            gradvx(i,:)  = gradvx(i,:)  + m(j)/rho(j)*(vx(j)-vx(i))*dwdr*eij;
            gradvy(i,:)  = gradvy(i,:)  + m(j)/rho(j)*(vy(j)-vy(i))*dwdr*eij;
            gradp(i,:)   = gradp(i,:)   + m(j)/rho(j)*(p(j)-p(i))*dwdr*eij;
            
            B = B - m(j)/rho(j)*dwdr*[rij_x;rij_y]*eij;
        end
        if abs(det(B)) > 1e-5 %prevent singular B
            invBt = (inv(B))';
            
            gradrho(i,:) = gradrho(i,:)*invBt;
            gradvx(i,:)  = gradvx(i,:)*invBt;
            gradvy(i,:)  = gradvy(i,:)*invBt;
            gradp(i,:)   = gradp(i,:)*invBt;
        end
    end
end

hessrho = zeros(N,4); %[rho_xx, rho_yx, rho_xy, rho_yy]
hessvx  = zeros(N,4);
hessvy  = zeros(N,4);
hessp   = zeros(N,4);
if settings(4) > 1
    parfor i = 1:N
        B = [0 0;0 0];
        for j = j_id{i}(2:end)
            rij_x = x(i) - x(j);
            rij_y = y(i) - y(j);
            norm_rij = sqrt(rij_x^2 + rij_y^2);
            eij = [rij_x, rij_y]/max( norm_rij, 0.001*kh );
            dwdr = wDer(norm_rij,kh,w_type);
            
            hessrho(i,:) = hessrho(i,:) + m(j)/rho(j)*[(gradrho(j,1)-gradrho(i,1))*eij,(gradrho(j,2)-gradrho(i,2))*eij]*dwdr;
            hessvx(i,:)  = hessvx(i,:) + m(j)/rho(j)*[(gradvx(j,1)-gradvx(i,1))*eij,(gradvx(j,2)-gradvx(i,2))*eij]*dwdr;
            hessvy(i,:)  = hessvy(i,:) + m(j)/rho(j)*[(gradvy(j,1)-gradvy(i,1))*eij,(gradvy(j,2)-gradvy(i,2))*eij]*dwdr;
            hessp(i,:)   = hessp(i,:) + m(j)/rho(j)*[(gradp(j,1)-gradp(i,1))*eij,(gradp(j,2)-gradp(i,2))*eij]*dwdr;
            
            B = B - m(j)/rho(j)*dwdr*[rij_x;rij_y]*eij;
        end
        if abs(det(B)) > 1e-5 %prevent singular B
            invB = inv(B);
            
            hessrho(i,:) = reshape( invB*reshape(hessrho(i,:),2,2) ,1,4);
            hessvx(i,:) = reshape( invB*reshape(hessvx(i,:),2,2) ,1,4);
            hessvy(i,:) = reshape( invB*reshape(hessvy(i,:),2,2) ,1,4);
            hessp(i,:) = reshape( invB*reshape(hessp(i,:),2,2) ,1,4);
        end
    end
end

parfor i = 1:N %parallel loop through each fluid particle
    
    for j = j_id{i}(2:end) %loop through neighbour fluid particles
        % Compute distance vectors
        rij_x = x(i) - x(j);
        rij_y = y(i) - y(j);
        norm_rij = sqrt(rij_x^2 + rij_y^2);
        eij_x = rij_x/max( norm_rij, 0.001*kh );
        eij_y = rij_y/max( norm_rij, 0.001*kh );
        dwdr = wDer(norm_rij,kh,w_type);
        
        % Viscosity
        if nu(i)==0 || nu(j)==0
            visc = 0;
        else
            visc = -2*(nu(i)+nu(j))*((vx(i)-vx(j))*eij_x+(vy(i)-vy(j))*eij_y)...
                /norm_rij/rho(i)/rho(j);
        end
        
        % Construct Riemann problem: 4-point stencil h--i--j--k
        dx2 = [rij_x^2;rij_x*rij_y;rij_x*rij_y;rij_y^2];
        rhoh = rho(i) + gradrho(i,:)*[rij_x;rij_y] + 0.5*hessrho(i,:)*dx2;
        vxh = vx(i) + gradvx(i,:)*[rij_x;rij_y] + 0.5*hessvx(i,:)*dx2;
        vyh = vy(i) + gradvy(i,:)*[rij_x;rij_y] + 0.5*hessvy(i,:)*dx2;
        ph = p(i) + gradp(i,:)*[rij_x;rij_y] + 0.5*hessp(i,:)*dx2;
        
        rhok = rho(j) - gradrho(j,:)*[rij_x;rij_y] + 0.5*hessrho(j,:)*dx2;
        vxk = vx(j) - gradvx(j,:)*[rij_x;rij_y] + 0.5*hessvx(j,:)*dx2;
        vyk = vy(j) - gradvy(j,:)*[rij_x;rij_y] + 0.5*hessvy(j,:)*dx2;
        pk = p(j) - gradp(j,:)*[rij_x;rij_y] + 0.5*hessp(j,:)*dx2;
        
        u = -[vxh,vx(i),vx(j),vxk]*eij_x -[vyh,vy(i),vy(j),vyk]*eij_y;
        [F_L,F_R] = recon([rhoh,rho(i),rho(j),rhok],u,[ph,p(i),p(j),pk],recon_order);
        
        % Solve Riemann problem for ij
        [us, ps] = riemannSolver(F_L,F_R,[c0(i) c0(j)],norm(f(:,i))*kh,settings(5));
        vsx = -(us - 0.5*(F_L(2)+F_R(2)))*eij_x + 0.5*(vx(i)+vx(j));
        vsy = -(us - 0.5*(F_L(2)+F_R(2)))*eij_y + 0.5*(vy(i)+vy(j));
        
        % Evaluate fluxes for i
        drhodt(i) = drhodt(i) + 2*rho(i)*m(j)/rho(j)*((vx(i)-vsx)*eij_x+(vy(i)-vsy)*eij_y)*dwdr;
        dvxdt(i) = dvxdt(i) - 2*m(j)/rho(i)/rho(j)*(ps+visc)*dwdr*eij_x;
        dvydt(i) = dvydt(i) - 2*m(j)/rho(i)/rho(j)*(ps+visc)*dwdr*eij_y;
    end
    
    for k = k_id{i} %loop through neighbour wall particles
        % Compute distance vectors
        rij_x = x(i) - xw(k);
        rij_y = y(i) - yw(k);
        norm_rij = sqrt(rij_x^2 + rij_y^2);
        eij_x = rij_x/max( norm_rij, 0.001*kh );
        eij_y = rij_y/max( norm_rij, 0.001*kh );
        dwdr = wDer(norm_rij,kh,w_type);
        
        % Velocity
        nw_x = n_w(k,1);
        nw_y = n_w(k,2);
        ui = -vx(i)*nw_x-vy(i)*nw_y;
        uw = vx(i)*nw_x+vy(i)*nw_y - 2*(vw_x(k)*nw_x+vw_y(k)*nw_y);
        
        % Pressure
        pw = p(i) + rho(i)*max( 0, -(f(1,i)*rij_x+f(2,i)*rij_y) );
        % Wall density
        rhow = density(pw,rho0(i),c0(i),gamma);
        
        % Reconstruction
        [F_L,F_R] = recon([0,rho(i),rhow,0],[0,ui,uw,0],[0,p(i),pw,0],0);
        
        % Solve Riemann problem
        [us, ps] = riemannSolver(F_L,F_R,[c0(i) c0(i)],norm(f(:,i))*kh,settings(5));
        vsx = -(us - 0.5*(F_L(2)+F_R(2)))*eij_x + 0.5*(vx(i)+2*vw_x(k)-vx(i));
        vsy = -(us - 0.5*(F_L(2)+F_R(2)))*eij_y + 0.5*(vy(i)+2*vw_y(k)-vy(i));
        
        % Viscosity
        if nu(i)==0
            visc = 0;
        else
            visc = -4*nu(i)*((vx(i)-vw_x(k))*eij_x+(vy(i)-vw_y(k))*eij_y)...
                /norm_rij/rho(i)/rhow;
        end
        
        % Evaluate fluxes
        drhodt(i) = drhodt(i) + 2*rho(i)*m(i)/rhow*((vx(i)-vsx)*eij_x+(vy(i)-vsy)*eij_y)*dwdr;
        dvxdt(i) = dvxdt(i) - 2*m(i)/rho(i)/rhow*(ps+visc)*dwdr*eij_x;
        dvydt(i) = dvydt(i) - 2*m(i)/rho(i)/rhow*(ps+visc)*dwdr*eij_y;
    end
end

% Body force
dvxdt = dvxdt + f(1,:)';
dvydt = dvydt + f(2,:)';

end