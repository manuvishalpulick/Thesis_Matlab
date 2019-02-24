function [t_rupt, H_Profile, T_Steps, Sk, f] = flatFilms(L,N,c,Tmp,gx,h_adjusted,A,p,endTime,seN,h0,save_threshold)
digits(32)
format longG

%% equation
%                                                           ___
%  d H      d  /   3 d^3 H     1 dH   \      __    d  /    /           \
%  --- = - ---|   H  ----- +   -----  |  + _/2T   ---|    / H^3 N(X,T)  |
%  d T     d X \     d X^3     H dX   /           d X \ \/             / 


%% boundary conditions

% periodic

%% discretization


%           j      j 
%        H^2  * H^2                      _            _   _           _
%          i-1     i                    |   1       1  | |   j    j    |
%h1_l = _------------ _     h2_l = 0.5* |  ---  +  --- | |  H  - H     |
%      |   j      j    |                |   j       j  | |_  i    i-1 _|
%      |  H    + H     |                |  H       H   |
%      |_  i-1    i   _|                |_  i-1     i _|  

%           j      j 
%        H^2  * H^2                      _            _   _           _
%          i+1     i                    |   1       1  | |   j      j  |
%h1_r = _------------ _     h2_r = 0.5* |  ---  +  --- | |  H    - H   |
%      |   j      j    |                |   j       j  | |_  i+1    i _|
%      |  H    + H     |                |  H       H   |
%      |_  i+1    i   _|                |_  i+1     i _|

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now here is the final discretization
%  _              
% |           j+1                           j+1                               j+1  
% | p2*h1_l* H        -p2*(3*h1_l + h1_r)* H        ( 3*p2*(h1_l+h1_r) + 1 ) H     ... 
% |_          i-2                           i-1                               i
%                                                _
%                         j+1               j+1   |                
%   -p2*(3*h1_r + h1_l)* H        p2*h1_r* H      |    =  LHS --> A matrix
%                         i+1               i+2  _|                       


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


DeltaX = L/N;                 % grid size
x = 0:DeltaX:L;               % domain of the simulation
DeltaT = DeltaX^c;            % time step size
Flag = 0;                     % counter for storing the data files
SaveCounter = 1;
p1 = DeltaT/(DeltaX^2);       % p1 - used in the explicit part
p2 = DeltaT/(DeltaX^4);       % p2 - used in the implicit part
p3 = 1/(DeltaX)*sqrt(2*DeltaT*Tmp);    % p3 - used for the noise term: Important to realize that factor 3 is not present in the avg mobility, if you do the proper scaling, so no need to have it in the sqrt here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%----------- Setting up the initial conditions ------------%%%

h = h0;                                   % initial film
h = [h(end-1), h(end), h, h(1), h(2)]';   % add two ghost points each side, and ensure periodicity
InitArea = trapz(x,h(3:end-2));           % initial area of the film (mass)

%%  preallocation

h1_r = zeros(size(h));
h1_l = zeros(size(h));
h2_r = zeros(size(h));
h2_l = zeros(size(h));
R = zeros(size(x));
noi1 = zeros(size(x));
noi2 = zeros(size(x));


k = linspace(1,h_adjusted,h_adjusted);   % index used in vectorizing the code; see the solver

%% Solver

rng shuffle                      % change the seed for every realization
for t = DeltaT:DeltaT:endTime    % time marching 
%     clear b
    R = randn(N+1,1);            % random numbers for each grid point (changes every time)
    gx_f = gx.*R;                % use the gx matrix as per the model
    g = sum(gx_f);               % Q-Wiener process definition as outlined in Grun et al, Diez et al, Lord et al
    g = [g(end-1); g(end); g'; g(1); g(2)];         % periodicity also in noise
    noi1 = (g(3:h_adjusted-2) + g(4:h_adjusted-1))./2;      % for the main domain
    noi2 = (g(2:h_adjusted-3) + g(3:h_adjusted-2))./2;      % for the main domain
    
    %% generate the pentagonal band in the sparse matrix
    
    % schemes outlined in Diez et al (2000), Grun et al (2004)[in appendix]
    % to discretize h^3 term
    h1_r = 2*h(k(3:h_adjusted-2)).^2.*h(k(4:h_adjusted-1)).^2./(h(k(3:h_adjusted-2))+h(k(4:h_adjusted-1)));          % see the discretization

    h1_l = 2*h(k(3:h_adjusted-2)).^2.*h(k(2:h_adjusted-3)).^2./(h(k(3:h_adjusted-2))+h(k(2:h_adjusted-3)));          % see the discretization

    h2_r =  0.5.*(1./h(k(3:h_adjusted-2)) + 1./h(k(4:h_adjusted-1))).*(h(k(4:h_adjusted-1)) - h(k(3:h_adjusted-2))); % see the discretization

    h2_l =  0.5.*(1./h(k(3:h_adjusted-2)) + 1./h(k(2:h_adjusted-3))).*(h(k(3:h_adjusted-2)) - h(k(2:h_adjusted-3))); % see the discretization

    A(p) =[h1_l*p2; -(h1_r+3*h1_l)*p2; 3*(h1_r+h1_l)*p2+1; -(3*h1_r+h1_l)*p2; h1_r*p2];    % see the discretization 
    
    %% generate the b-vector in Ax = b

    b = [0; 0; h(3:h_adjusted-2)-p1*(h2_r - h2_l) + p3.*(sqrt(h1_r).*noi1 - sqrt(h1_l).*noi2); 0; 0];                  % see the discretization

    h(:) = A\b;                                                                                                           % solve
    if min(h(:)) <= 0.01    % if the film height goes below a certain height stop the realization; it is insensitive to a value below 0.1. Rapid thinning
        break;
    end
%======================================================================
%               Save data after every seN time steps
%======================================================================
    Flag = Flag + 1;
    if t == DeltaT                                                          % Conditional statement to check if it is the first iteration
        h_previous = h;
    elseif (min(h_previous)-min(h))*(100/min(h_previous))<=save_threshold
        %fprintf('Difference percent = %f; Save counter = %d\n',(min(h_previous)-min(h))*(100/min(h_previous)),save_counter);
        h_previous = h;
        if(Flag==seN)
            Flag=0;
            H_Profile(:,SaveCounter) = h(:);
            T_Steps(SaveCounter) = t;
            [Sk(:,SaveCounter),f(:,SaveCounter)] = CalculateSk(h(:),L,N);
            SaveCounter = SaveCounter + 1;
        end
    else
        %fprintf('Difference percent = %f; Save counter = %d\n',(min(h_previous)-min(h))*(100/min(h_previous)),save_counter);
        if(Flag==seN)
            Flag=0;
        end
        H_Profile(:,SaveCounter) = h(:);
        T_Steps(SaveCounter) = t;
        SaveCounter = SaveCounter + 1;
        h_previous = h;
    end
        
%  some check to ensure if any realization has given singular matrix      
%         if(rcond(full(A)) < 1e-12)
%             break
%         else
%             continue
%         end

end
% save_counter
% d1 = digits
t_rupt = t;         % rupture time obtained from this simulation
end

