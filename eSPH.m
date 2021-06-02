function eSPH(inputf)
% A 2D Riemann solver based weakly-compressible SPH code
% written by Enson Tin Hang Un for MEng FYP, 2021
% Department of Aeronautics, Imperial College London
%
% Dependencies: density.m, pressure.m, recon.m, riemannSolver.m,
%               stepSize.m, timeDer.m, wallNormal.m, wDer.m (total 8)

%% Initialisation
% Inputs
load(inputf,'fluid','wall','f','settings','dt_save','dir_name');

% Pre-compute wall normal vector (static wall, i.e. wall particle x,y fixed)
if size(wall,1) > 0
    n_w = wallNormal(wall,settings(1),settings(2));
else
    n_w = nan;
end

% Initialise saving settings
if exist(dir_name,'dir') == 7
    warning(['"',dir_name,'" already exist.'])
    usrinput = input('Overwrite existing data - 1 \nContinue simulation     - 2 \nTerminate - something else:\n');
    if (usrinput ~= 1)&&(usrinput ~= 2)
        error('Terminated.')
    end
else
    mkdir(dir_name) %create output directory
end

% Initialise time loop
n_dt = 0;
n_save = 1;
if (exist('usrinput','var') == 0) 
    %initialise directly if no existing data or overwrite existing data
    t = 0;    
elseif usrinput == 1
    t = 0;
else    
    %start from the last existing data
    while exist(dir_name+"/SPHout_"+num2str(n_save)+".mat","file") == 2
        n_save = n_save + 1;
    end
    load(dir_name+"/SPHout_"+num2str(n_save-1)+".mat","t","fluid","wall");
end

%% Time loop
tic
while t < settings(7)
    % Find dt
    dt = stepSize(max(fluid(:,10)), ... max(c0)
                  max(sqrt(fluid(:,6).^2 + fluid(:,7).^2)), ... norm(v)
                  settings(2), ... h
                  max( sqrt(sum(f(fluid(:,1),fluid(:,2),t).^2)) ), ... norm(f)
                  max(fluid(:,11)),... max(nu)  
                  settings(8)); %CFL#
    if dt <= 1e-10
        save([dir_name,'/',dir_name,'_end.mat'],'fluid','wall','t');
        error('dt is too small')
    end
    
    % Time integration: 2nd order Stormer-Verlet
    if n_dt == 0
        [~, dvxdt, dvydt] = timeDer(fluid,wall,settings(1:6),f(fluid(:,1),fluid(:,2),t),n_w);
    end
    fluid(:,6) = fluid(:,6) + 0.5*dvxdt*dt;
    fluid(:,7) = fluid(:,7) + 0.5*dvydt*dt;
    
    fluid(:,1) = fluid(:,1) + 0.5*fluid(:,6)*dt;
    fluid(:,2) = fluid(:,2) + 0.5*fluid(:,7)*dt;
    
    t = t + 0.5*dt;
    [drhodt, ~, ~] = timeDer(fluid,wall,settings(1:6),f(fluid(:,1),fluid(:,2),t),n_w);
    fluid(:,3) = fluid(:,3) + drhodt*dt;
    fluid(:,5) = pressure(fluid(:,3),fluid(:,9),fluid(:,10),settings(3));
    
    fluid(:,1) = fluid(:,1) + 0.5*fluid(:,6)*dt;
    fluid(:,2) = fluid(:,2) + 0.5*fluid(:,7)*dt;
    
    t = t + 0.5*dt;
    [~, dvxdt, dvydt] = timeDer(fluid,wall,settings(1:6),f(fluid(:,1),fluid(:,2),t),n_w);
    fluid(:,6) = fluid(:,6) + 0.5*dvxdt*dt;
    fluid(:,7) = fluid(:,7) + 0.5*dvydt*dt;
    
    n_dt = n_dt + 1; %timestep++
    
    % Write particles to output file
    if t >= n_save*dt_save
        % Save all particle as .mat-file
        save([dir_name,'/SPHout_',num2str(n_save),'.mat'],'fluid','wall','t');
        
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
    % Warning if density variation > 1%
    if max(abs(fluid(:,3)./fluid(:,9)-1)) > 0.01
        warning(['Density flucation = ',num2str(max(abs(fluid(:,3)./fluid(:,9)-1))*100,'%.2f'),'%'])
    end
end

runtime = toc;
fprintf('Simulation completed after %.2fs \n',runtime);

%% Save settings to a txt file
txtID = fopen([dir_name,'/detail.txt'],'w');
fprintf(txtID,'%s\n',datetime);
fprintf(txtID,'No. of fluid particles: %d\n',size(fluid,1));
fprintf(txtID,'No. of wall particles : %d\n',size(wall,1));
fprintf(txtID,'%7s = %.3f\n',[["kernel","h","gamma","recon","riemann","beta","cfl#"];settings([1:6,8])]);
fprintf(txtID,'runtime: %.4fs\n',runtime);
fclose(txtID);