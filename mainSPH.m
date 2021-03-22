% A 2D Riemann solver based weakly compressible SPH code
% MEng FYP, Enson Tin Hang Un, 2021
% Department of Aeronautics, Imperial College London
%
% Inputs:
% fluid - particles(:,1) = x coordinate
%         particles(:,2) = y coordinate
%         particles(:,3) = density
%         particles(:,4) = mass (constant throughout simulation)
%         particles(:,5) = pressure (only used when output)
%         particles(:,6) = x velocity
%         particles(:,7) = y velocity
%         particles(:,8) = /
%         particles(:,9) = total density of fluid, rho0
%         particles(:,10) = artifical speed of sound, c0
%         particles(:,11) = kinematic viscosity, nu
% wall -- particles(:,1) = x coordinate
%         particles(:,2) = y coordinate
%         particles(:,4) = mass (constant throughout simulation)
%         particles(:,5) = wall x velocity
%         particles(:,6) = wall y velocity
% settings -- settings(1) = kernel function type
%             settings(2) = h
%             settings(3) = gamma for equation of state (usually =7)
%             settings(4) = reconstruction scheme
%             settings(5) = Riemann solver type
%             settings(6) = dissipation limiter parameter (usually = 3)
%             settings(7) = time integration order
%             settings(8) = simulation end time
%             settings(9) = max dt
% ref_values - ref_values(:,1) = 
%              ref_values(:,2) = artifical speed of sound, c0
%              ref_values(:,3) = kinematic viscosity, nu
% f(t) ------- forcing per unit mass = [fx;fy]
% dt_save --- time between each output file (may not be exact)
% file_name - output file prefix
%
% Functions:
% dwdr(r,type,h) -------------- parameter : type = kernel function
%                                           h = kernel compact support length
%                               input     : r = norm(distance)
%                               output    : dwdr = derivative of the kernel function
% pressure(rho,rho0,c0,gamma) - parameters: rho0 = initial density
%                                           c0 = speed of sound (=10*v_max)
%                                           gamma = tunnable parameter (usually =7)
%                               inputs    : rho = density (=particles(:,3))
%                               output    : p = pressure (=particles(:,5))
% rho(pressure,rho0,c0,gamma) - parameters: rho0 = initial density
%                                           c0 = speed of sound (=10*v_max)
%                                           gamma = tunnable parameter (usually =7)
%                               inputs    : p = pressure (=particles(:,5))
%                               output    : rho = density (=particles(:,3))
% recon(particles,type) ------- parameter : type = reconstruction scheme
%                               inputs    : ?
%                               outputs   : F_L, F_R for riemannSolver(..)
% boundaryRecon(particles) ---- dependency: dwdr(r), rho(pressure)
%                               inputs    : ?
%                               outputs   : F_L, F_R for riemannSolver(..)
% riemannSolver(F_L,F_R,type) - parameter : type = which Riemann solver
%                                           eta = dissipation limiter parameter (~3)
%                               inputs    : F_L, F_R from reconstruction
%                               outputs   : p_star, u_star for drhodt and dvdt
% timeDer(particles,nu) ------- dependency: dwdr, recon, boundaryRecon,
%                                           riemannSolver, rangesearch
%                               inputs    : m,rho,v,r,u_star,p_star,nu
%                               outputs   : drhodt, dvxdt, dvydt
% stepSize(c0,norm(v),h,f,nu,maxdt) - parameter : maxdt = max step size
%                                     inputs    : c0,norm(v),h,f,nu
%                                     outputs   : dt
% timeMarch(particles,type) --- dependency: timeDer, stepsize
%                               parameter : order = order of scheme
%                                           T = end time
%                               inputs    : particles
%                               outputs   : fluid particles v, r, rho
clear


%% Initialisation
% Inputs
load('input.mat');

% Pre-compute wall normal vector (static wall, i.e. wall particle x,y fixed)
if size(wall,1) > 0
    n_w = wallNormal(wall,settings(1),settings(2));
else
    n_w = nan;
end
%quiver(wall(:,1),wall(:,2),n_w(:,1),n_w(:,2))

% Initialise saving settings
if file_name == 0
    file_name = 'SPHout';
end
if exist(file_name,'dir') == 7
    warning(['"',file_name,'" already exist.'])
    if input('Overwrite existing data? (YES - 1) \n') ~= 1
        error('Terminated.')
    end
else
    mkdir(file_name) %create output directory
end
% Save settings to a txt file
txtID = fopen([file_name,'/detail.txt'],'w');
fprintf(txtID,'%s\n',datetime);
fprintf(txtID,'No. of fluid particles: %d\n',size(fluid,1));
fprintf(txtID,'No. of wall particles : %d\n',size(wall,1));
fprintf(txtID,'%7s = %.3f\n',[["kernel","h","gamma","recon","riemann","beta"];settings(1:6)]);
% vtuSave(particles, 0, file_name); %save boundary particles

% Initialise time loop
t = 0;
n_dt = 0;
n_save = 1;
order = settings(7);

%% Time loop
tic
%for nn = 1:10
while t < settings(8)
    % Find dt
    dt = stepSize(max(fluid(:,10)), ... max(c0)
                  max(sqrt(fluid(:,6).^2 + fluid(:,7).^2)), ... norm(v)
                  settings(2), ... h
                  max( sqrt(sum(f(fluid(:,1),fluid(:,2),t).^2)) ), ... norm(f)
                  max(fluid(:,11)), ... max(nu)
                  settings(9)); % max_dt              
    if dt <= 1e-10
        save([file_name,'/',file_name,'_end.mat'],'fluid','wall','t');
        error('dt is too small')
    end
    
    % Time integration
    if order == 2 %2nd order: Stormer-Verlet
        if n_dt == 0
            [~, dvxdt, dvydt] = timeDer(fluid,wall,settings,f(fluid(:,1),fluid(:,2),t),n_w);
        end
        fluid(:,6) = fluid(:,6) + 0.5*dvxdt*dt;
        fluid(:,7) = fluid(:,7) + 0.5*dvydt*dt;
        
        fluid(:,1) = fluid(:,1) + 0.5*fluid(:,6)*dt;
        fluid(:,2) = fluid(:,2) + 0.5*fluid(:,7)*dt;
        
        t = t + 0.5*dt;
        [drhodt, ~, ~] = timeDer(fluid,wall,settings,f(fluid(:,1),fluid(:,2),t),n_w);
        fluid(:,3) = fluid(:,3) + drhodt*dt;
        fluid(:,5) = pressure(fluid(:,3),fluid(:,9),fluid(:,10),settings(3));
        
        fluid(:,1) = fluid(:,1) + 0.5*fluid(:,6)*dt;
        fluid(:,2) = fluid(:,2) + 0.5*fluid(:,7)*dt;
        
        t = t + 0.5*dt;
        [~, dvxdt, dvydt] = timeDer(fluid,wall,settings,f(fluid(:,1),fluid(:,2),t),n_w);
        fluid(:,6) = fluid(:,6) + 0.5*dvxdt*dt;
        fluid(:,7) = fluid(:,7) + 0.5*dvydt*dt;
        
    end
    n_dt = n_dt + 1;
    
    % Output particles to file
    if t >= n_save*dt_save
        % Save all particle as .mat-file
        save([file_name,'/SPHout_',num2str(n_save),'.mat'],'fluid','wall','t');
        
        % Save fluid particle matrix as .vtu-file
        %vtuSave(particles, n_save, file_name);
        
        n_save = n_save + 1;
    end
    
    % Display progress 
    fprintf('t = %.4f, dt = %.6f, step %d\n',t,dt,n_dt); 
    
    % User defined error message
    if exist('catch_err','var') == 1
        for err_n = 1:size(catch_err,2)
            if catch_err(1,err_n)
                error(catch_err(2,err_n))
            end
        end
    end
%     if any(fluid(:,2)<-2)
%         error('Fluid leaks!')
%     end
    if max(abs(fluid(:,3)./fluid(:,9)-1)) > 0.01 %check 1% density flucation holds
        warning(['Density flucation = ',num2str(max(abs(fluid(:,3)./fluid(:,9)-1))*100,'%.2f'),'%'])
    end
end

% Save final workspace

runtime = toc;
fprintf('Simulation completed after %.2fs \n',runtime);

fprintf(txtID,'runtime: %.4fs\n',runtime);
fclose(txtID);