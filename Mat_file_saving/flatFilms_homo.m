
function [t_rupt,k_dom_sim,t_calc] = flatFilms_repulsion_save(L,N,deltaX,c,B,Tmp,gx,h_adjusted,A,p,endTime,seN,N_Reals,strhet,realization,animationSkip,continue_index);


format long g
P_het= 0;
e = 0;



%% equation
%                                                                            ___
%  d H      d  /   3 d^3 H     1 dH       B   dH    \       __    d  /    /           \
%  --- = - ---|   H  ----- +   ----- -  ------ --    |  + _/2T   ---|    / H^3 N(X,T)  |
%  d T     d X \     d X^3     H dX       H^2  dX   /           d X \ \/              /


%% boundary conditions

% periodic

%% discretization
%           j      j 
%       2. H^2  * H^2                      _            _   _           _
%          i-1     i                    |   1       1  | |   j    j    |
%h1_l = _------------ _     h2_l = 0.5* |  ---  +  --- | |  H  - H     | 
%      |   j      j    |                |   j       j  | |_  i    i-1 _|
%      |  H    + H     |                |  H       H   |
%      |_  i-1    i   _|                |_  i-1     i _|  

%                 _                _   _           _
%                |   h_m     h_m    | |   j    j    |
%  h3_l =B* 0.5* |  ---  +  ---     | |  H  - H     | 
%                |     j         j  | |_  i    i-1 _|
%                |  H^2       H^2   |
%                |_    i-1       i _|  

%                 _                _   _           _
%                |   1       1      | |   j    j    |
%  h3_r = B*0.5* |  ---  +  ---     | |  H  - H     | 
%                |     j         j  | |_  i+1  i   _|
%                |  H^2       H^2   |
%                |_    i        i+1_|  

%           j      j 
%       2. H^2  * H^2                      _            _   _           _
%          i+1     i                    |   1       1  | |   j      j  |
%h1_r = _------------ _     h2_r = 0.5* |  ---  +  --- | |  H    - H   |
%      |   j      j    |                |   j       j  | |_  i+1    i _|
%      |  H    + H     |                |  H       H   |
%      |_  i+1    i   _|                |_  i+1     i _|

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now here is the final discretization
%  _              
% |         j+1                             j+1                              j+1  
% | p2*h2_l* H        -p2*(3*h1_l + h1_r)* H        ( 3*p2*(h1_1+h1_r) + 1 ) H     ... 
% |_         i-2                            i-1                               i
%                                          _
%                     j+1             j+1   |                
%   -p2*(3*h1 + h2)* H        p2*h1* H      |    =  LHS --> A matrix
%                     i+1             i+2  _|                       


% _                                                                               _          
%| 0                                                                               | 
%| 0                                                                               | 
%|                                                                                 | 
%|previous        disj pressure                   noise term                       |     
%| soln             vdW only                                                       |
%|  j                                                                              |
%| H       -     p1 (h2_r - h2_l)  +    p3*(sqrt(h1_r).*noi1 - sqrt(h1_l).*noi2)
%|  i                                                                              |        
%|                                                                                 |  
%|0                                                                                |  
%|_0                                                                              _| 

tInside_solver = tic;
%deltaX = L/N;                 % grid size
if continue_index == 0
    x = 0:deltaX:L;               % domain of the simulation
    deltaT = deltaX^c            % time step size
    flag = 1;                     % counter for storing the data files
    p1 = deltaT/(deltaX^2);       % p1 - used in the explicit part
    p2 = deltaT/(deltaX^4);       % p2 - used in the implicit part
    p3 = 1/(deltaX)*sqrt(2*deltaT*Tmp);    % p3 - used for the noise term: Important to realize that factor 3 is not present in the avg mobility, if you do the proper scaling, so no need to have it in the sqrt here
    %% Solver
    rupture=0;              % index for rupture event
    t_rupt=0;
    iter =1;  % index for iteration

    k_dom_sim=1; %initialization
    k_dom_index =1 ;
    lin_index = 0;
    h_flat = ones(max(size(x)),1);        % Flat film initialization for linear regime determination
    f =(2*pi.*(linspace(0,(N)/2,(N+2)/2)./L))';
    % endTime = 2;
    stable = 1;                     % Stability index
    het=0;
    rng shuffle                      % change the seed for every realization

    tElapsed_strt = toc(tInside_solver);
    fprintf('Time taken to start the solver flatFilms solver: %d min %f s\n',floor(tElapsed_strt/60),mod(tElapsed_strt,60));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%----------- Setting up the initial conditions ------------%%%
    %P=(N/6*    deltaX);
    %h = ones(size(x));                        %  flat initial film
    %%
    %f =(2*pi.*(linspace(0,(N/2),(N+2)/2)./L))';     %Frequency distribution




    %sinsum = zeros(size(x));
    % for ini = 1:N/2+1
    %     s = sin(x.*f(ini));
    %     sinsum(:) = sinsum(:)+s' ;
    %  end
    %h(1,:) = ones(size(x))+0.001.*(sinsum);%+ 0.01.*(sinsum(:,5)+sinsum(:,10)+sinsum(:,2)+sinsum(:,N-2))';
    %h=ones(size(x))+0.001*sin(6.*(x-7).^2);  %for L=15
    max_iter = round(endTime/deltaT) +1;

    %h_save=ones(max_iter,max(size(x)))';  % initializing saved h for speed

    % introducing an initial perturbation whichexcites multiple frequencies
    %h = ones(size(x))+0.001*sin(6.*(x-(L/2)).^2); 
    h = ones(size(x))+0.001.*sin(0.2.*x);
    h = [h(end-1), h(end), h, h(1), h(2)]';   % add two ghost points each side, and ensure periodicity employing periodic boundary conditions
    h_save(:,1) = h;
    % Energy of initial film
        E_st(1) = deltaX.*sum(1.5.*((h(4:h_adjusted-1,1)-h(2:h_adjusted-3,1))./(2.*deltaX)).^2);
        E_pi(1) = deltaX.*sum((1./(2.*(h(3:h_adjusted-2,1)).^2)) - (0.1./(3.*(h(3:h_adjusted-2,1)).^2)));
        E(1) = E_st(1) + E_pi(1);% Y_diff=h(3:h_adjusted-2)- mean(h(3:h_adjusted-2));

    %     Y_diff = h(3:end-2)-mean(h(3:end-2));
    %     hk=fft(Y_diff);
    %     P2 = abs(hk/(L)); % Taking only the real values of the fourier transform and normalizing
    %     %which can still give us the required info
    %     %P2=fftshift(hk2);
    %     P1 = P2(1:(N/2)+1); %since the middle value would be N/2 +1
    %     P1(2:end-1) = 2*P1(2:end-1);   %Have to double the magnitude as we 
    %     %ignored one side of the domain in the previous step
    %     %why is it 2:end-1 and not 1:end --> found in net dont know why
    % %    k= ((N/L).*(linspace(0,(N)/2,(N)/2+1))/N)';
    %     f =(2*pi.*(linspace(0,(N/2),(N+2)/2)./L))';     %Frequency distribution
    % %    f =(2*pi.*(linspace(0,(N+1),(N+1))./L))'
    %     S=P1.^2;
    %     plot(f(2:end-1),S(2:end-1),'*','linewidth',2);
    %     xlabel('k [-]','Fontsize',16)
    %     ylabel('(S_k) [-]','Fontsize',16)
    %     %ylim([0 8e-6]);
    %     %xlim([0 100])
    %  figure
    %      plot(x,h(3:h_adjusted-2))
    %  %ylim([0.998 1.002]);
    %  xlim([0 L])
    %  xlabel('x')
    %  ylabel('h')
    %  title('initial height profile')
    % figure
    %%
    %h = [h(end-1), h(end), h, h(1), h(2)]';   % add two ghost points each side, and ensure periodicity employing periodic boundary conditions
    %InitArea = trapz(x,h(3:end-2));           % initial area of the film (mass) Just for confirmation

    %%  preallocation

    h1_r = zeros(size(h));
    h1_l = zeros(size(h));
    h2_r = zeros(size(h));
    h2_l = zeros(size(h));
    h3_r = zeros(size(h));
    h3_l = zeros(size(h));
    R = zeros(size(x));
    noi1 = zeros(size(x));
    noi2 = zeros(size(x));
    lin_time = 0;
    B=0.1;  % h_min--> minimum height (pre-defined precursor film height) constant of Part of disjoining pressure
    %B=0; % For No repulsion case
    k = linspace(1,h_adjusted,h_adjusted);   % index used in vectorizing the code; see the solver
    t_start = deltaT;
else
    
    load('temp.mat')
    h = h_save;
    x = 0:deltaX:L;               % domain of the simulation
    deltaT = deltaX^c            % time step size
    flag = 1;                     % counter for storing the data files
    p1 = deltaT/(deltaX^2);       % p1 - used in the explicit part
    p2 = deltaT/(deltaX^4);       % p2 - used in the implicit part
    p3 = 1/(deltaX)*sqrt(2*deltaT*Tmp);    % p3 - used for the noise term: Important to realize that factor 3 is not present in the avg mobility, if you do the proper scaling, so no need to have it in the sqrt here
    %% Solver
    rupture=0;              % index for rupture event
    t_rupt=0;
    iter =1;  % index for iteration

    h_flat = ones(max(size(x)),1);        % Flat film initialization for linear regime determination
    f =(2*pi.*(linspace(0,(N)/2,(N+2)/2)./L))';
    % endTime = 2;
    stable = 1;                     % Stability index
    het=0;
    rng shuffle                      % change the seed for every realization
    max_iter = round(endTime/deltaT) +1;

    h1_r = zeros(size(h));
    h1_l = zeros(size(h));
    h2_r = zeros(size(h));
    h2_l = zeros(size(h));
    h3_r = zeros(size(h));
    h3_l = zeros(size(h));
    R = zeros(size(x));
    noi1 = zeros(size(x));
    noi2 = zeros(size(x));
    B=0.1;  % h_min--> minimum height (pre-defined precursor film height) constant of Part of disjoining pressure
    %B=0; % For No repulsion case
    k = linspace(1,h_adjusted,h_adjusted);   % index used in vectorizing the code; see the solver
    t_start = t + deltaT;
    tElapsed_strt = toc(tInside_solver);
    fprintf('Time taken to start the solver flatFilms solver: %d min %f s\n',floor(tElapsed_strt/60),mod(tElapsed_strt,60));
end

t_loop =tic;
t_iter = tic;
for t = t_start:deltaT:endTime    % time marching 
%     clear b
   
    if Tmp == 0
        noi1 = 0;      % for the main domain
        noi2 = 0;      % for the main domain 
    else
        R = randn(N+1,1);            % random numbers for each grid point (changes every time)
        gx_f = gx.*R;                % use the gx matrix as per the model
        g = sum(gx_f);               % Q-Wiener process definition as outlined in Grun et al, Diez et al, Lord et al
        
        g = [g(end-1); g(end); g'; g(1); g(2)];         % periodicity also in noise
        noi1 = (g(3:h_adjusted-2) + g(4:h_adjusted-1))./2;      % for the main domain
        noi2 = (g(2:h_adjusted-3) + g(3:h_adjusted-2))./2;      % for the main domain
    end
    
    %% generate the pentagonal band in the sparse matrix
    %tic
    % schemes outlined in Diez et al (2000), Grun et al (2004)[in appendix]
    % to discretize h^3 term
    h1_r = 2*h(k(3:h_adjusted-2)).^2.*h(k(4:h_adjusted-1)).^2./(h(k(3:h_adjusted-2))+h(k(4:h_adjusted-1)));          % see the discretization

    h1_l = 2*h(k(3:h_adjusted-2)).^2.*h(k(2:h_adjusted-3)).^2./(h(k(3:h_adjusted-2))+h(k(2:h_adjusted-3)));          % see the discretization

    h2_r =  0.5.*(1./h(k(3:h_adjusted-2)) + 1./h(k(4:h_adjusted-1))).*(h(k(4:h_adjusted-1)) - h(k(3:h_adjusted-2))); % see the discretization

    h2_l =  0.5.*(1./h(k(3:h_adjusted-2)) + 1./h(k(2:h_adjusted-3))).*(h(k(3:h_adjusted-2)) - h(k(2:h_adjusted-3))); % see the discretization

    h3_r =  -0.5*(4/3)*B.*(1./(h(k(3:h_adjusted-2))).^2 + 1./(h(k(4:h_adjusted-1))).^2).*(h(k(4:h_adjusted-1)) - h(k(3:h_adjusted-2))); % see the discretization

    h3_l =  -0.5*(4/3)*B.*(1./(h(k(3:h_adjusted-2))).^2 + 1./(h(k(2:h_adjusted-3))).^2).*(h(k(3:h_adjusted-2)) - h(k(2:h_adjusted-3))); % see the discretization

    
%      h3_r=(1./h(k(3:h_adjusted-2)) + 1./h(k(2:h_adjusted-3)));
%      h4_r=(h(k(3:h_adjusted-2)) - h(k(2:h_adjusted-3)));
%      h2_2=h3_r.*h4_r.*h3_r;

    A(p) =[h1_l*p2; -(h1_r+3*h1_l)*p2; 3*(h1_r+h1_l)*p2+1; -(3*h1_r+h1_l)*p2; h1_r*p2];   % see the discretization - Arranging the 5 element space in the middle 
    %of the matrix through vectorization- p was already defined, each element in the RHS is of size=N+1 and p is of size N+1*5 (5--> -3,-2,-1,0,1) 
    
    %% generate the b-vector in Ax = b
    
    b = [0; 0; h(3:h_adjusted-2)-p1*(h2_r - h2_l+ h3_r- h3_l) + p3.*(sqrt(h1_r).*noi1 - sqrt(h1_l).*noi2); 0; 0];                  % see the discretization
    % First and last 2 elements will be 0 as its is for the periodic
    % boundary conditions (eg: h(end)-h(1)=0 and so on)
%    b = [0; 0; h(3:h_adjusted-2)-p1*(h2_r - h2_l) + p3.*(sqrt(h1_r).*noi1 - sqrt(h1_l).*noi2); 0; 0]; 
    iter = iter+1;
    h = A\b;                                                                                                           % solve for h
    %toc
    h_save(:,iter) = h;
    
   h_min=min(h);
   h_max = max(h);
   
     if h_min <= 0.11 && rupture==0    % if the film height goes below a certain height stop the realization; 
        %it is insensitive to a value below 0.1. Rapid thinning
        t_rupt = t         % rupture time obtained from this simulation
        rupture=1;          % indicating rupture event has happened
     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Energy determination %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%     [E,E_pi,E_st]=energy_calc(x, L_flat, deltaX, c, deltaT, N, endTime, P_het, e, Tmp, realization);
%     E_st(iter) = -deltaX.*sum(1.5.*((h(4:h_adjusted-1)-h(2:h_adjusted-3))./(2.*deltaX)).^2);
%     E_pi(iter) = -deltaX.*sum((1./(2.*(h(3:h_adjusted-2)).^2)) - (0.1./(3.*(h(3:h_adjusted-2)).^3)));
%     E(iter) = E_st(iter) + E_pi(iter);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
     if min(h(:)) <= 0.0001    % if the film height goes below a certain height stop the realization; it is insensitive to a value below 0.1. Rapid thinning
         stable = 0;
         break; % here we break out if the simulation and call it the rupture event and the time at this moment is the rupture time
        
     end
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

% toc
%% to determine the dominant wave number
    if (h_min <= 0.995 || h_max >= 1.005) && lin_index == 0
           lin_index =iter;
            lin_time = t
            [S,f] = CalculateSk_f(h(3:h_adjusted-2),L,N) ;                      
           [S_max,S_max_index(1,iter)] = max(S);  % Assign the indexes of max values
           k_dom_index = S_max_index(1,iter);
    end
   %toc  
   
%====================================================================== 
%               Save data after every seN time steps
%======================================================================
    flag = flag+1;  % count for saving data 
    % save data eevery seN time steps
    if flag == seN
        %tElapsed_iter = toc(tInside_solver);
        %fprintf('Time taken for seN iterations: %d min %f s\n',floor(tElapsed_iter/60),mod(tElapsed_iter,60));
        flag = 1;                     % counter for storing the data files
        %t                            % Display t after every seN time steps if required to know the extent of simulation  
        %tic
        save('temp.mat','-v7.3');
        %toc
        %t_iter = tic;
    end
    %toc
    
end
tElapsed_loop = toc(tInside_solver);
fprintf('\nTime taken to complete %d iterations: %d min %f s\n', iter, floor(tElapsed_loop/60),mod(tElapsed_loop,60));

k_dom_sim = f(k_dom_index)   % Identifying the dominant wave number in the defined linear regime

% If omega has to be calculated
%omega_max_sim = k_dom_sim_calc_save(animationSkip,iter,h_save(3:end-2,:),x,deltaX,c,L,time,het,P_het,wave_dom_lsa,N,realization,lin_index,k_dom_sim,e,Tmp)

filename = [strhet,'T_rupt',num2str(t_rupt),'_Lf_',num2str(L),'_deltaX_',num2str(deltaX),'_c_',num2str(c), '_Tmp_', num2str(Tmp),'_P_het_', num2str(P_het), '_e_', num2str(e),'rzn',num2str(realization),'.mat']; 

save(filename,'-v7.3');

end

