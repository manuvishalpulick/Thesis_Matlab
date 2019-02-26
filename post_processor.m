function post_processor(animationSkip, x, tt, L_flat, deltaX, c, deltaT, N, endTime, het, P_het, wave_dom_lsa, e, Tmp, N_reals,strhet)

    for realization= 1:N_reals
           t_load= tic;
           filename = [strhet,'_Lf_',num2str(L_flat),'_deltaX_',num2str(deltaX),'_c_',num2str(c), '_Tmp_', num2str(Tmp),'_P_het_', num2str(P_het), '_e_', num2str(e),'rzn',num2str(realization),'.mat'];
           
           %%%%%%%%%%% Only to be used to postprocess old data where T_rupt was a part of the filename%%%%%%%%%%%%%%%%%
           %filename = [strhet,'T_rupt0','_Lf_',num2str(L_flat),'_deltaX_',num2str(deltaX),'_c_',num2str(c), '_Tmp_', num2str(Tmp),'_P_het_', num2str(P_het), '_e_', num2str(e),'rzn',num2str(realization),'.mat'];
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
           mk = strcat(strhet,'_Lf_',num2str(L_flat),'_deltaX_',num2str(deltaX),'_c_',num2str(c), '_Tmp_', num2str(Tmp),'_P_het_', num2str(P_het), '_e_', num2str(e));         
           str1 = strcat('.\',mk,'\',filename);
           load(str1);
           tElapsed_load=toc(t_load);
           fprintf('\nTime taken to load the results for realization &d: %d min %f s\n',floor(tElapsed_load/60),mod(tElapsed_load,60),realization);
           tic
           h_adjusted = N+5;
           
           %% eliminating the ghost points
           h_final = [h_save((3:end-2),:)];
           %h_final = [h((3:end-2),:)];   % If analysis is done on an old
           %set of points
           q = size(h_final);  % t(q) = endTime if the simulation was stable and ran till endTime
           % q is used to prevent errors if the simulation was not stable
           q=q(2);
           E=zeros(N_reals,q);
           E_pi=zeros(N_reals,q);
           E_st=zeros(N_reals,q);
           %%If energy has to be calculated during post processing; needs
           %%ghost points also to calculate first and last values, so
           %%h_save is passed instead of h_final
           toc    
           [E(realization,:),E_pi(realization,:),E_st(realization,:)]=energy_calc(h_save, x, L_flat, deltaX, c, deltaT, N, q, P_het, e, Tmp, realization,endTime);
           
           t=[0:deltaT:endTime];
           
        %% function calls
        %% Function for making the plot of minimum height vs time
            %hmin_evol(h_min,t,realization,q)
        %% to make an animation of the film evolution
            %makeAnimation_det(animationSkip,q,c1,x,L_flat,t_new,het,P_het,e)
        %% to make an animation of the fourier spectrum evolution    
            %makeAnimation_det_fft(animationSkip,q,c1,x,L_flat,t_new,N,P_het,wave_dom,t)
        %% Makes a subplot of both film and structure factor evolution    
            makeAnimation_comparison(animationSkip,q,h_final,x,L_flat,t,het,P_het,wave_dom_lsa,e,N,realization)
        %% to make a plot of the energy evolution and different contributions to the same
            %energy_evol(h_final,E,E_pi,E_st,x,L_flat,t,het,P_het,wave_dom_lsa,e,N,realization,q)
        %% Function to plot and find the maximum growth rate     
            %omega_max_sim = k_dom_sim_calc_save(animationSkip,q,h_final,x,deltaX,c,L_flat,t,het,P_het,wave_dom_lsa,N,realization,lin_index,k_dom_sim,e,Tmp)
        %% Function to do an error analysis on the growth rate vs k plot to determine the linear regime   
            %disp_rel(h_final,L_flat,N,lin_index,t,q,animationSkip)
        %% Function to trace the interfaces to find the pinning time; to be done only for heterogeneous analysis
        if het == 1
            %interface_1(animationSkip,q,h,x,L_flat,deltaX,t,het,P_het,wave_dom_lsa,e,N,realization,h_adjusted,t_rupt,c)
        end

    end
    %save(filename)
   
end

%% for detailed comments, check make_animation_comparison
function makeAnimation_det(animationSkip,q,h_final,x,L_flat,t,het,P_het,e)
%% to specify animation conditions

v = VideoWriter('filmEvolution');
v.FrameRate = 5;  % Default 30
v.Quality = 100;    % Default 75
open(v)
hfig = figure;
o= hfig.WindowState;
hfig.WindowState='maximize';
%l=length(c(:,1));
% x=linspace(0,L_flat,l);
% N1=length(x)
for i = 1:animationSkip:q
    Y = h_final(:,i);
    area(x',Y)
    ylim([0 10])
    xlim([0 L_flat])
    xlabel('x [-]','Fontsize',16)
    ylabel('h [-]','Fontsize',16)
%    ax.XMinorGrid = 'on'
    set(gca,'FontSize',18)
   
     if het == 1
        strhet='heterogeneous';
        hold on
        plot(x',e.*cos((2*pi*x/P_het)-pi),'-')
        xlim([0 L_flat])
        ylim([0 10])
%        legend(t_new(i),'wettability')
        legend(t(i))
        hold off
    else
         legend(num2str(t(i)))
         strhet='homogeneous';
    end
    M(i) = getframe(hfig);        % Everytime getframe captures the entire 
    % figure properties correspoding to hfig and stores it in M(i) which
    % can be written in a video afterwards
end
M = M(~cellfun(@isempty,{M.cdata}));  % Ignore all null vectors generated because of animation skipping
writeVideo(v,M)
close(v)


end

function makeAnimation_det_fft(animationSkip,q,h_final,x,L_flat,t,N,P_het,wave_dom_lsa)
%% to specify animation conditions

vid = VideoWriter('Fourier_evolution');
vid.FrameRate = 5;  % Default 30
vid.Quality = 100;    % Default 75
open(vid)
hfig = figure;
%l=length(c(:,1));
%deltaX=(L_flat/N);
o= hfig.WindowState;
hfig.WindowState='maximize';



for i = 1:animationSkip:q
    Y = h_final(:,i);
    Y_diff=Y-mean(Y);
    %% trying to anlayse in fourier domain
    hk=fft(Y_diff);  % N+1 values
%    area(x',Y)
    P2=abs(hk/(L_flat));
%    P2 = fftshift(hk2) ;
    P1 = P2(1:N/2+1); %since the middle value would be N/2 +1
    P1(2:end-1) = 2*P1(2:end-1) ;  % found in net dont know why
%    k= 2*pi.*(1./([2*deltaX:2*deltaX:L_flat]))';
    f =(2*pi.*(linspace(0,(N)/2,(N+2)/2)./L_flat))' ;
    S=P1.^2;
    
%    semilogy(f(2:end-1),(S(2:end-1)),'*','linewidth',2);
%   subplot(2,1,1)
    loglog(f(2:end-1),(S(2:end-1)),'*','linewidth',2);
%    semilogy(f,(S),'-','linewidth',2);
    xlabel('k [-]','Fontsize',16)
    ylabel('(S_k)[-]','Fontsize',16)
    xlim([0 55]);
    ylim([10e-40 10e3]);
 %    ax.XMinorGrid = 'on'
    set(gca,'FontSize',18)
   legend(num2str(t(i)))
   F(i) = getframe(hfig);        % Everytime getframe captures the entire 
    % figure properties correspoding to hfig and stores it in M(i) which
    % can be written in a video afterwards 
 
%    %% Displaying and plotting the stripe period harmonics
%        S_evolution(:,iter)= S(index_het); % Recognizing the harmonics
 %%   
%    iter= iter+1;
end

f(index_het)  % To display the harmonic wave numbers if necessary
F = F(~cellfun(@isempty,{F.cdata}));
writeVideo(vid,F)
close(vid)

end

%% Funtion to make a subplot of both film and structure factor evolution   
function makeAnimation_comparison(animationSkip,q,h_final,x,L_flat,t,het,P_het,wave_dom_lsa,e,N,realization)
%% to specify animation conditions

v = VideoWriter(strcat('Evolution_subplot',num2str(realization)));
v.FrameRate = 5;  % Default 30
v.Quality = 100;    % Default 75
open(v)
%% to make the figure full screen ; works in MATLAB 2018 and above
hfig = figure;
o= hfig.WindowState;
hfig.WindowState='maximize';

%% initializing indexes for harmonic tracing
anim_index=[1:animationSkip:q];  % index for animation time
t_anim = t(anim_index); % time for animation
S_evolution= zeros(5,max(size(t_anim)));  % Initializing structure factor evolution array for harmonic evolution
iter=1;  %index for iteration
harmonic = [0.5,1,2,3,4];  % Will account for the first 4 harmonics(heterogeneous)/ (+/-) for homogeneous 
        %generating indexes for harmonics
    if het==0              
        index_het(:) = (round((L_flat/wave_dom_lsa)+harmonic(:))+1);  % +1 to account for the first element
    else
        index_het(:) = (round((L_flat/P_het).*harmonic(:))+1);   % +1 to account for the first element
    end
%%    
S_max_index = zeros(1,max(size(t_anim))); % index for maximum structure factor evolution
rupt=0;   %  rupture index ; =1 after rupture
        
for i = 1:animationSkip:q
    Y = h_final(:,i);  % Assigning h at timestep to Y
    subplot(2,1,1)
    area(x',Y)         % Plot film evolution
    ylim([0 10])
    xlim([0 L_flat])
    xlabel('x [-]','Fontsize',16)
    ylabel('h [-]','Fontsize',16)
%    ax.XMinorGrid = 'on'
    set(gca,'FontSize',18)
   
    % Show the heterogeneities if heterogeneous
    if het == 1  % to plot the heterogeneities
        strhet='heterogeneous';
        hold on
        cos_values(:) = e.*cos((2*pi*x/P_het)-pi);
        
        hetcount = (cos_values < 0);  % generates a matrix with ones on more wettable sites
        hetcount = -1.*hetcount;
        plot(x',hetcount,'LineStyle','-','color','red','LineWidth',4)
        xlim([0 L_flat])
        ylim([0 10])
%        legend(t_new(i),'wettability')
        legend(strcat('t=',num2str(t(i))))
        hold off
    else
         legend(strcat('t=',num2str(t(i))))
         strhet='homogeneous';
    end

    [S,f] = CalculateSk_f(Y,L_flat,N); % Calculate the structure factor
%% plot the structure factor with length and make animation    
    subplot(2,1,2)
    loglog(f(2:end-1),(S(2:end-1)),'*','linewidth',2);
    xlabel('k [-]','Fontsize',16)
    ylabel('(S_k)[-]','Fontsize',16)
    xlim([0 10e2]);
    ylim([10e-30 10e3]);
    xlabel('k [-]','Fontsize',16)
    ylabel('Sk [-]','Fontsize',16)
    %ax.XMinorGrid = 'on'
    set(gca,'FontSize',18)
    legend(strcat('t=',num2str(t(i))))   
    
    if rupt == 0 && min(Y) <= 0.11  % save rupture instant
        savefig(strcat('Rupture_instant',num2str(realization)))
        saveas(hfig,strcat('Rupture_instant',num2str(realization),'.jpg'))
        rupt=1;
    end
    
    B(i) = getframe(hfig);        % Everytime getframe captures the entire 
    % figure properties correspoding to hfig and stores it in M(i) which
    % can be written in a video afterwards
    [S_max,S_max_index(1,iter)] = max(S);  % Assign the indexes of max values to plot the max frequency evolution
        
    %% Displaying and plotting the stripe period harmonics
        S_evolution(:,iter)= S(index_het); % Recognizing the harmonics and asssigning values
 %%   
    iter= iter+1;  % iteration count
end
B = B(~cellfun(@isempty,{B.cdata}));  % to neglect the values which are blank (due to animationSkip)
writeVideo(v,B)
close(v)
 
savefig(strcat('Final_Sk_snapshot',num2str(realization))) % save the final morphology snapshot
%% Plot the max frequency evolution- Gives an idea about the dominant wavenumber
maxfreq_fig=figure;
plot(t_anim(2:end),f(S_max_index(2:end)))
xlabel('time')
ylabel('frequency with max Structure factor')
title('Frequency with max fourier energy evolution')
savefig(strcat('Maximum_energy_frequency',num2str(realization)))
saveas(maxfreq_fig,strcat('Maximum_energy_frequency',num2str(realization),'.jpg'))

%% Plot the harmonic evolution and save
evolfig = figure ;

for j=1:5
    
    semilogy(t_anim,(S_evolution(j,:)),'linewidth',2)
    xlabel('time')
    ylabel('Structure factor of harmonics')
    title('Structure factor evolution of relevant wave numbers')
%    ylim([10e-6 10e03]);
%    xlim([0 30])
    hold on;
end
%legendstring = [ strcat('h1',f(index_het(1))) , strcat('h2',f(index_het(2))), strcat('h3',f(index_het(3))), strcat('h4',f(index_het(4))), strcat('h5',f(index_het(5))) ]
legend('0.5 harmonic','1st harmonic','2nd harmonic','3rd harmonic','4th harmonic','location','southeast')
savefig(strcat('harmonic_evolution',num2str(realization)))
saveas(evolfig,strcat('harmonic_evolution',num2str(realization),'.jpg'))

end



