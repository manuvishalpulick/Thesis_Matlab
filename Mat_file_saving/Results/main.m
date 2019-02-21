function main
clc
close all
clearvars -global
clearvars

T_StartMain = tic;
T_StartMatAssembly = tic;

L_FlatArray = [60];
k1 = 0.2;
k2 = 0;

DeltaX = 0.05859375;                                                        % grid/mesh size

C = 2.75;                                                                   % exponent used in deciding deltaT = deltaX^c
Tmp = 0;                                                                    % dimensionless noise strength (= 0, for deterministic)
Kappa = 0;                                                                  % dimensionless curvature (= 0 for flat films)

EndTime = 5000;                                                             % end time of a realization
seN = 20;                                                                   % save every these many time steps
N_Reals = 1;                                                                % number of realizations. 1 in case of deterministic simulation
SaveThreshold = 5;
N_Sim = size(L_FlatArray,2);

RunFolders = dir(fullfile('Results','Run*'));
if isempty(RunFolders)
    RunNo = 0;
else
    for FolderNo = 1:max(size(RunFolders))
        f_name = RunFolders(FolderNo).name;
        if contains(f_name,'run','IgnoreCase',true)
            Segments = split(f_name,' ');
            S_Index = contains(Segments,'run','IgnoreCase',true);
            RunNo(FolderNo) = str2double(Segments(find(S_Index,1)+1));
        else
            RunNo(FolderNo) = 0;
        end
    end
end

if ~exist('Run','var')
    Run = max(RunNo)+1;
else
    if isempty(Run)
        Run = max(RunNo)+1;
    else
        if max(ismember(RunNo,Run))
            error(['Given input Run number is %d\nBut a folder with '...
                'the same Run number already exists in the Results ' ...
                'folder! \nPlease try again with a different Run '...
                'number.'],Run);
        end
    end
end

[~, ~] = mkdir(fullfile('Results',strcat("Run ",num2str(Run)),'Height profiles'));
MatFiles = dir(fullfile('Results',strcat("Run ",num2str(Run)),'Height profiles','*.mat'));
if isempty(MatFiles)
    SimNo = 0;
else
    for FileNo = 1:max(size(MatFiles))
        f_name = MatFiles(FileNo).name;
        if contains(f_name,"sim",'IgnoreCase',true)
            Segments = split(f_name,'_');
            S_Index = contains(Segments,"sim",'IgnoreCase',true);
            SimNo(FileNo) = str2double(Segments(find(S_Index,1)+1));
        else
            SimNo(FileNo) = 0;
        end
    end
end

for SimLocal = 1:N_Sim
    SimGlobal = SimLocal + max(SimNo);
    L_flat = L_FlatArray(SimLocal);

    if Kappa == 0  
        L_curv = 0;                         % conditions for flat films
        x = 0:DeltaX:L_flat;                % domain for flat films
    else
        L_curv = 300;                       % length of the curved portion of the film, for kappa > 1, one needs a smaller value of of L_curv
        x1 = (-L_curv ):DeltaX:-DeltaX;     % domain for the curved portion of the film without one ghost point that gets added later on. This is used now to generate noise
        x2 = 0:DeltaX:(L_flat);             % domain for the flat portion of the film without the two ghost points. Same as above
        x = [x1 x2];                        % full domain
    end
    
    L = L_curv + L_flat;                    % total length of the film (curved+flat)
    N = round(L/DeltaX);                    % adjusted number of grid points -- different from earlier value of N only for curved films 
    N_Adjusted = N+5;                       % 2 ghost points on each side and extra 1 for the additional grid point (for x=0)
    
    gx = gx_generator(N,L,x);               % generates a matrix that is going to be used when we finally implement noise
    h0 = initial_condition(k1,k2,x);              % initial film
%%
    Realization = 0;                        % counter for the number of realizations
%%
    if Kappa == 0   
        for Reals = 1:N_Reals                   % perform 'NReals' number of realizations
            
            % preallocate sparse matrix below: the right most entry is for
            % identifying how many non-zero elements do we have in sparse matrix
            A=spalloc( N_Adjusted , N_Adjusted ,8 + (N_Adjusted-4)*5) ;        % preallocate sparse matrix: 5 nonzeros points for (N+1) or (h_adjusted - 4) grid points, and 2 for each points (2*4--> 8) on the boundary conditions
            %% generate a 'p' vector that is going to be used to fill up the band in the sparse matrix
            % in C++ or in fortran, we dont have to do it in this complicated fashion, just a simple for loop would be effective
            p1 = linspace(3, N_Adjusted-2,N_Adjusted-4)*(1 + N_Adjusted);    
            p  = [];
            k  = linspace(1,N_Adjusted,N_Adjusted);
            for i =-3:1     % -3: for leftmost/lowest band, -1: for the diagonal and 1: for the rightmost/upper band 
                p=[p p1+i*N_Adjusted ] ;    
            end
        
%% now fill up the boundary conditions
            p_ones= [1 N_Adjusted+2 N_Adjusted^2 N_Adjusted*(N_Adjusted - 1)-1];
            A(p_ones)=ones (1,4);
            p_minusOnes=[N_Adjusted*2+N_Adjusted-1 N_Adjusted*3+N_Adjusted (N_Adjusted-4)*N_Adjusted+1 (N_Adjusted-3)*N_Adjusted+2];
            A(p_minusOnes)=ones (1,4) *-1;
        
            tElapsed_matAssembly = toc(T_StartMatAssembly);
            fprintf('Time taken to assemble the empty matrix: %d min %f s\n',floor(tElapsed_matAssembly/60),mod(tElapsed_matAssembly,60));
            tStart_flatFilms = tic;
%% call the solver
            [T_rupt, H_Profile, T_Steps, Sk, f] = flatFilms(L,N,C,Tmp,gx,N_Adjusted,A,p,EndTime,seN,h0,SaveThreshold);
            
            tElapsed_flatFilms = toc(tStart_flatFilms);
            fprintf('Time taken by the flatFilms solver: %d min %f s\n',floor(tElapsed_flatFilms/60),mod(tElapsed_flatFilms,60));
            tStart_create_mat = tic;
        
%% store the data (but more importantly the rupture times) into a mat file, so that there is no further post processing required if we are just looking for T_r
            filename = ['Workspace_','Run_',num2str(Run),'_Simulation_',num2str(SimGlobal),...
                '_Realisation_',num2str(Reals),'_Lf_',num2str(L),...
                '_Tmp_', num2str(Tmp),'.mat'];
            save(filename)
            movefile(filename,fullfile('Results',strcat("Run ",num2str(Run)),'Height profiles'));
            
            tElapsed_create_mat = toc(tStart_create_mat);
            fprintf('Time taken to create the .mat file: %d min %f s\n',floor(tElapsed_create_mat/60),mod(tElapsed_create_mat,60));
            tStart_makeAnimation = tic;
            %PostProcessOptions = {[animationSkip threshold res], 0, 0, 0};
            PostProcessOptions = {[20 5 0], 0, 0, 0};
            %PostProcessing(Run,SimGlobal,Reals,PostProcessOptions)
            
            tElapsed_makeAnimation = toc(tStart_makeAnimation);
            fprintf('Time taken to make the animation: %d min %f s\n',floor(tElapsed_makeAnimation/60),mod(tElapsed_makeAnimation,60));

    end
    %% following is for non-flat simulations (quite similar to the flat ones, except for the boundary conditions
    else
        for Reals = 1:N_Reals      
            N_Adjusted= N + 4; % extra 1+2 for the ghost point and extra 1 for the additional grid point
            A=spalloc( N_Adjusted , N_Adjusted ,10 + N_Adjusted *5-10) ;% preallocate sparse matrix
            %% generate a 'p' vector that is going to be used to fill up band in the sparse matrix
            p1=linspace (3 , N_Adjusted-2,N_Adjusted-4) *(1 + N_Adjusted ) ;
            p = [] ;
            k = linspace(1,N_Adjusted,N_Adjusted);
            for i =-3:1
                p=[p p1+ i*N_Adjusted ] ;
            end
            %% boundary conditions : dh/dx = 0 and d3h/dx3 = 0 for right far field || and h = 1 + kappa*x^2 and d2h/dx2 = 2 kappa for the far left (curved) film
            p_ones=[1 2*N_Adjusted+1 N_Adjusted+2 N_Adjusted*(N_Adjusted-3)+(N_Adjusted-1) N_Adjusted^2];
            A(p_ones )=ones(1,5);
            p_minusOnes=[N_Adjusted*(N_Adjusted-4) (N_Adjusted-1)*(N_Adjusted)+(N_Adjusted-1)];
            A(p_minusOnes)=ones (1,2) *-1;
            A(N_Adjusted+1) = -2;
            A(N_Adjusted*(N_Adjusted-3)) = 3;
            A(N_Adjusted*(N_Adjusted-1)) = -3;
            %% call the solver

            T_rupt(Reals) = nonFlatFilms(L,L_flat,L_curv,N,C,Tmp,gx,N_Adjusted,A,p,EndTime,seN,DeltaX,Kappa);
        
            reali_series(Reals) = Reals;
            Realization = Realization + 1;
            mk = strcat('realization',num2str(Realization));  % name your realization folder
            mk2 = mkdir(mk);                                  % make its directory
            movefile('Data*',mk)                              % move all the data files to that directory
        end
    end
    
    tElapsed_main = toc(T_StartMain);
    fprintf('Total time taken by the run: %d min %f s\n',floor(tElapsed_main/60),mod(tElapsed_main,60));
end
end
