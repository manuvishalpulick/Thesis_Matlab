function main

clear all
close all
clc

k_dom_lsa = 0.6586;     % LSA prediction of dom wave number  
ho=1;
error_index = 1;  % Assign 1 if you are performing error analysis
%% initializing plot properties
omega = @(k) (-ho.^3.*k.^4) + (k.^2.*(1/ho)) - (0.4/3.*k.^2.*(1/ho.^2));  %dispersion relation
marker = ['*','o','+','d','.'];
colour = ['r','g','b','c','b'];
%% Simulation parameters

    L_flat_array = [20, 40, 60, 80, 100]; 
    deltaX = [0.025,0.05,0.1,0.125];
    %deltaX = 0.125;
    %legend_string = [num2str(deltaX(1)),num2str(deltaX(2)),num2str(deltaX(3)),num2str(deltaX(4))];
    c = [2.5: 0.25 : 3.25];     %DeltaT parameter
    L_flat = L_flat_array(1);   %length of domain
    Tmp=0.000; % Dimensionless Noise 
%% Heterogeneity initialization    
e=0.0;
P_het=0;   
if e == 0.0
            strhet='homogeneous';
            het=0;
else
            strhet='heterogeneous';
            het=1;
end
%%    
    
    for iter_p = 1:4 % for deltaX change
%         delta_X(iter_p) = L_flat/(2.^p(iter_p));
        for iter_c = 1:4 % for deltaT change
            [t_ruptavg(iter_p,iter_c), k_dom_sim(iter_p,iter_c),omega_max_sim(iter_p,iter_c),t_calc_avg(iter_p,iter_c)]=thin_films(L_flat,deltaX(iter_p),c(iter_c),P_het,e,Tmp);
            deltaT(iter_p,iter_c) = deltaX(iter_p)^(c(iter_c));
        end
        %% To ignore all zero entries
%         deltaT(:,all(deltaT==0)) = [];
%         deltaT(all(deltaT==0,2),:) = [];
%         t_ruptavg(:,all(t_ruptavg==0)) = [];
%         t_ruptavg(all(t_ruptavg==0,2),:) = [];
%         k_dom_sim(:,all(k_dom_sim==0)) = [];
%         k_dom_sim(all(k_dom_sim==0,2),:) = [];
%         omega_max_sim(:,all(omega_max_sim==0)) = [];
%         omega_max_sim(all(omega_max_sim==0,2),:) = [];
       
        if error_index == 1        
            truptfig = figure(1);
            plot(deltaT(iter_p,:),t_ruptavg(iter_p,:),'Color',colour(iter_p),'Marker',marker(iter_p))
            xlabel('deltaT')
            ylabel('time of rupture')
            title('Time of rupture change with deltaT')
            hold on
            %legend(strcat('deltaX=',num2str(deltaX(1))),strcat('deltaX=',num2str(deltaX(2))),strcat('deltaX=',num2str(deltaX(3))),strcat('deltaX=',num2str(deltaX(4))),'Location','southeast');
            legend(strcat('deltaX=',num2str(deltaX(iter_p))),'Location','southeast');
            err_k(iter_p,:) = k_dom_sim(iter_p,:) - k_dom_lsa;
            err_omega_norm(iter_p,:) = sqrt((omega(k_dom_lsa) - omega_max_sim(iter_p,:))).^2;
            err_kfig = figure(2);
            plot(deltaT(iter_p,:),err_k(iter_p,:),'Color',colour(iter_p),'Marker',marker(iter_p))
            xlabel('deltaT')
            ylabel('error in dominant wave number')
            title('Wave number error change with deltaT')
            %legend(strcat('deltaX=',num2str(deltaX(1))),strcat('deltaX=',num2str(deltaX(2))),strcat('deltaX=',num2str(deltaX(3))),strcat('deltaX=',num2str(deltaX(4))),'Location','southeast');
            legend(strcat('deltaX=',num2str(deltaX(iter_p))),'Location','southeast');
            hold on
            err_omegafig = figure(3);
            plot( deltaT(iter_p,:),err_omega_norm(iter_p,:),'Color',colour(iter_p),'Marker',marker(iter_p))
            xlabel('deltaT')
            ylabel('error in growth rate')
            title('Growth rate error change with deltaT')
           %legend(strcat('deltaX=',num2str(deltaX(1))),strcat('deltaX=',num2str(deltaX(2))),strcat('deltaX=',num2str(deltaX(3))),strcat('deltaX=',num2str(deltaX(4))),'Location','southeast');
            legend(strcat('deltaX=',num2str(deltaX(iter_p))),'Location','southeast');
            hold on
            %[err_min,min_index] = min(err_omega_norm(iter_p,:));
            %c_min(iter_p) = c(min_index);
            t_calc_fig = figure(4);
            plot( deltaT(iter_p,:),t_calc_avg(iter_p,:),'Color',colour(iter_p),'Marker',marker(iter_p))
            xlabel('deltaT')
            ylabel('Time of simulation')
            title('Simulation time with deltaT')
           %legend(strcat('deltaX=',num2str(deltaX(1))),strcat('deltaX=',num2str(deltaX(2))),strcat('deltaX=',num2str(deltaX(3))),strcat('deltaX=',num2str(deltaX(4))),'Location','southeast');
            legend(strcat('deltaX=',num2str(deltaX(iter_p))),'Location','northeast');
        end
    end
    
    if error_index ==1
        save(strcat('Errors',num2str(L_flat),'.mat'),'err_k','err_omega_norm','t_ruptavg','k_dom_sim')
        savefig(strcat('truptfig',num2str(L_flat)))
        savefig(strcat('err_kfig',num2str(L_flat)))
        savefig(strcat('err_omegafig',num2str(L_flat)))
        savefig(strcat('t_calcfig',num2str(L_flat)))
        mk2 = strcat(strhet,'_Lf_',num2str(L_flat), '_Tmp_', num2str(Tmp),'_P_het_', num2str(P_het), '_e_', num2str(e));
        movefile('Errors*',mk2)
        movefile('err_*',mk2)
        movefile('truptfig*',mk2)
        movefile('t_calcfig*',mk2)
    end
    
end

function [t_ruptavg ,t_calcavg] = thin_films(L_flat,deltaX,c,P_het,e,Tmp)

kappa = 0.0;                                                                % dimensionless curvature (= 0 for flat films)

%deltaX = 0.05                                                              % grid/mesh size
if kappa == 0  
    L_curv = 0;      % conditions for flat films
    x = 0:deltaX:L_flat;        % domain for flat films
else
    L_curv = 300;       % length of the curved portion of the film, for kappa > 1, one needs a smaller value of of L_curv
    x1 = (-L_curv ):deltaX:-deltaX;    % domain for the curved portion of the film without one ghost point that gets added later on. This is used now to generate noise
    x2 = 0:deltaX:(L_flat);            % domain for the flat portion of the film without the two ghost points. Same as above
    x = [x1 x2];                       % full domain
end
L = L_curv + L_flat;   % total length of the film (curved+flat)  Still L_flat for flat films
N = round(L/deltaX);   % adjusted number of grid points -- different from earlier value of N only for curved films 
deltaT = deltaX^c;          % time step size
deltaX
%% initialization of heterogeneous parameters
wave_dom_lsa=9.54;   % Prediction from theory
Pc = wave_dom_lsa./sqrt(2);     % Critical wavelength from theory
%P_het= 0;          %Periodicity of heterogeneity
ratio_het = P_het/Pc;          % Ratio of Phet:Pc
%e=0.0;  %Amplitude of wettability
%P_het= ratio_het.*Pc;          %Periodicity of heterogeneity

N_nodes_het = round(N/(L_flat/P_het));  % No of nodes in one heterogeneous stripe


if e == 0.0
            strhet='homogeneous';
            het=0;
else
            strhet='heterogeneous';
            het=1;
end
%% Noise intro
%Tmp = 0.000;             % dimensionless noise strength (= 0, for deterministic)       
if Tmp == 0.000
    N_Reals = 1;           % number of realizations for deterministic
else
    N_Reals = 25;           % number of realizations for stochastic
end
gx = gx_generator(N,L,x);  % generates a matrix that is going to be used when we finally implement noise
%% Processing parameters
post_pro=1;                   % if only postprocessing has to be done, set to 1
animationSkip = 800;        % To fix after how many time steps should the animation take the next plot values
continue_index = 0;         % If we want to continue a previous simulation; assign 1 else 0
continue_index_post = 0;   % If the data files have already been read and saved to a .mat file once, 
%set to 1 to avoid reading them again 
% endTime = 45;          % end time of a realization
if Tmp == 0   
    if e < 0.1 && e > 0
        endTime = 40;          % experimentation
    elseif e ==0
        endTime = 60;
    else
        endTime = 20;
    end
else
    if e == 0
        endTime = 60;          % experimentation
    else
        endTime = 30;
    
    end
end

seN = 5000;              % save every these many time steps    
realization = 0;                % counter for the number of realizations
t_rupt = zeros(N_Reals,1);      % preallocate rupture times vector
tt = seN*deltaT;                % time between saving two files

    if post_pro == 0
        %% Allocation of space for making A and B matrices for solving
                h_adjusted=N+5;         % 2 ghost points on each side and extra 1 for the additional grid point (for x=0)
                % preallocate sparse matrix below: the right most entry is for
                % identifying how many non-zero elements do we have in spars mat
                A=spalloc( h_adjusted , h_adjusted ,8 + (h_adjusted-4)*5) ;        
                % preallocate sparse matrix: 5 nonzeros points for (N+1) or (h_adjusted - 4) grid points, and 2 for each points (2*4--> 8) on the boundary conditions

         %Generate a 'p' vector that is going to be used to fill up the band in the sparse matrix

                % in C++ or in fortran, we dont have to do it in this complicated fashion, just a simple for loop would be effective
                p1 = linspace(3, h_adjusted-2,h_adjusted-4)*(1 + h_adjusted);    
                p  = [];
                k  = linspace(1,h_adjusted,h_adjusted);
                for i =-3:1     % -3: for leftmost/lowest band, -1: for the diagonal and 1: for the rightmost/upper band 
                    p=[p p1+i*h_adjusted ] ;    
                end

                %% now fill up the boundary conditions
                p_ones= [1 h_adjusted+2 h_adjusted^2 h_adjusted*(h_adjusted - 1)-1];
                A(p_ones)=ones (1,4);
                p_minusOnes=[h_adjusted*2+h_adjusted-1 h_adjusted*3+h_adjusted (h_adjusted-4)*h_adjusted+1 (h_adjusted-3)*h_adjusted+2];
                A(p_minusOnes)=ones (1,4) *-1;
        %%
        if kappa == 0                     %For flat films
            for m = 1:N_Reals           % perform 'NReals' number of realizations
         %% call the solver 
                tFlatfilms_solver = tic;
                % Here we get the rupture time of the realization
                if strcmp(strhet,'homogeneous')== 1
                    
                   %[t_rupt] = flatFilms_repulsion_gen(L,N,deltaX,c,Tmp,gx,h_adjusted,A,p,endTime,seN,N_Reals,realization);      % for general simulations
                   %toc
                   [t_rupt(m), k_dom_sim(m)] =flatFilms_homo(L,N,deltaX,c,Tmp,gx,h_adjusted,A,p,endTime,seN,N_Reals,strhet,m,animationSkip,continue_index);    % For Homogeneous with reuplsion
                    
                else
        %            t_rupt(m) = flatFilms_het_repulsion_2(L,N,c,Tmp,gx,h_adjusted,A,p,endTime,seN,N_nodes_het,P_het,e);  % For Heterogeneous with reuplsion
        %           t_rupt(m) = flatFilms_het(L,N,c,Tmp,gx,h_adjusted,A,p,endTime,seN);           %For heterogeneous without repulsion    
                    [t_rupt(m)] = flatFilms_het(L,N,deltaX,c,Tmp,gx,h_adjusted,A,p,endTime,seN,N_nodes_het,P_het,e,continue_index, N_Reals,strhet,m,animationSkip);
                end
                t_calc(m) = toc(tFlatfilms_solver);
                fprintf('Time taken by the flatFilms solver: %d min %f s\n',floor(t_calc(m)/60),mod(t_calc(m),60))
%                 reali_series(m) = m;                    % to keep a log of the realization (not needed, but I keep it)
%                 realization = realization + 1;          % counter
            
            end
            %% following is for non-flat simulations (quite similar to the flat ones, except for the boundary conditions
        else
            for m = 1:N_Reals      
                h_adjusted= N + 4; % extra 1+2 for the ghost point and extra 1 for the additional grid point
                A=spalloc( h_adjusted , h_adjusted ,10 + h_adjusted *5-10) ;% preallocate sparse matrix
                %% generate a 'p' vector that is going to be used to fill up band in the sparse matrix
                p1=linspace (3 , h_adjusted-2,h_adjusted-4) *(1 + h_adjusted ) ;
                p = [] ;
                k = linspace(1,h_adjusted,h_adjusted);
                for i =-3:1
                    p=[p p1+ i*h_adjusted ] ;
                end
                %% boundary conditions : dh/dx = 0 and d3h/dx3 = 0 for right far field || and h = 1 + kappa*x^2 and d2h/dx2 = 2 kappa for the far left (curved) film
                p_ones=[1 2*h_adjusted+1 h_adjusted+2 h_adjusted*(h_adjusted-3)+(h_adjusted-1) h_adjusted^2];
                A(p_ones )=ones(1,5);
                p_minusOnes=[h_adjusted*(h_adjusted-4) (h_adjusted-1)*(h_adjusted)+(h_adjusted-1)];
                A(p_minusOnes)=ones (1,2) *-1;
                A(h_adjusted+1) = -2;
                A(h_adjusted*(h_adjusted-3)) = 3;
                A(h_adjusted*(h_adjusted-1)) = -3;
                %% call the solver

                t_rupt(m) = nonFlatFilms(L,L_flat,L_curv,N,c,Tmp,gx,h_adjusted,A,p,endTime,seN,deltaX,kappa);

                reali_series(m) = m;
                realization = realization + 1;
                %post_processor_saved(animationSkip, x, tt, L_flat, deltaX, c, deltaT, N, endTime, t_rupt(m), het, P_het, wave_dom_lsa, e, Tmp, N_Reals,strhet);  

            end
        end
        tpp = tic;
         t_ruptavg= mean(t_rupt);
         fprintf('Rupture time  for domain %d, P_het %d, e %d is %d\n',L_flat,P_het,e, t_ruptavg)
         S_t_rupt = std(t_rupt);;
         fprintf('Std dev. of rupture time  for domain %d, P_het %d, e %d is %d\n',L_flat,P_het,e, S_t_rupt)
         t_calcavg = mean(t_calc);
         fprintf('Simulation time  for domain %d, P_het %d, e %d is %d\n',L_flat,P_het,e, t_calcavg)
         S_t_calc = std(t_calc);
         fprintf('Std dev. of simulation time for domain %d, P_het %d, e %d is %d\n',L_flat,P_het,e, S_t_calc)

         if het == 0
             k_dom_sim_avg = mean(k_dom_sim);
             fprintf('dom. wave number for domain %d, P_het %d, e %d is %d\n',L_flat,P_het,e, k_dom_simavg)
             S_k_dom_sim = std(k_dom_sim);
             fprintf('Std. dav. of dom. wave number for domain %d, P_het %d, e %d is %d\n',L_flat,P_het,e, S_k_dom_sim)
%              omega_max_sim_avg = mean(omega_max_sim);
%              fprintf('dom_growth rate for domain %d, P_het %d, e %d is %d\n',L_flat,P_het,e, omega_max_simavg)
%              S_omega_max_sim = std(omega_max_sim);
%              fprintf('Std. dav. of dom growth rate for domain %d, P_het %d, e %d is %d\n',L_flat,P_het,e, S_omega_max_sim)
         else 
             k_dom_sim_avg = 0;
             S_k_dom_sim = 0;
             omega_max_sim_avg = 0;
             S_omega_max_sim = 0;
         end
         
        %% store the data (but more importantly the rupture times) into a .mat file, so that there is no further post processing required if we are just looking for T_r
        filename = ['simulationdata_','kappa_',num2str(kappa),'_Lf_',num2str(L),'_N_',num2str(N), '_Tmp_', num2str(Tmp),'.mat'];
        save(filename,'t_ruptavg','k_dom_sim_avg','S_t_rupt','S_k_dom_sim','t_calcavg','S_t_calc')   % saves the .mat file including all the variable values in a filename of the specified format
        
        mk = strcat(strhet,'_Lf_',num2str(L_flat),'_deltaX_',num2str(deltaX),'_c_',num2str(c), '_Tmp_', num2str(Tmp),'_P_het_', num2str(P_het), '_e_', num2str(e));  % name your realization folder
        mk2 = mkdir(mk);                                  % make its directory
        movefile('*.mat',mk)                              % move all the data files to that directory
        a=0;
        R_f = 65e-6;                     % radius of flat surface
        h0_init = 150e-9;                % initial height =150nm
        A_vw = 2.026e-20;
        gam = 0.034;                     %for length scale calculation
        Rc = 1.8e-3;
        visc = 0.00089;
        %% derived quantities
        % t_scale = 12*pi^2*visc*gam*h0_init^5/A_vw^2;
        % l_scale = h0_init^2*sqrt(2*pi*gam/A_vw); 
        %continue_index_post =0;
        
        move_results(mk)
        
        tElapsed_pp = toc(tpp);
        fprintf('Time taken for post processing: %d min %f s\n',floor(tElapsed_pp/60),mod(tElapsed_pp,60))
    else
        for realization = 1:N_Reals  
            %t_ruptavg=28.0164; % If simulation data file was not created
            file = strcat('*rzn',realization)
            mk = strcat(strhet,'_Lf_',num2str(L_flat),'_deltaX_',num2str(deltaX),'_c_',num2str(c), '_Tmp_', num2str(Tmp),'_P_het_', num2str(P_het), '_e_', num2str(e));
            str2 = strcat('./',mk,'/',file);
            load(str2,t_rupt(realization)); % reading rupture time to pass to post processor

            %load(filename)
            post_processor_saved(animationSkip, x, tt, L_flat, deltaX, c, deltaT, N, endTime, t_rupt(realization), het, P_het, wave_dom_lsa, e, Tmp, N_Reals,strhet);
            move_results(mk)
        end
    end

end