function post_processor_saved(animationSkip, x, tt, L_flat, deltaX, c, deltaT, N, endTime, t_ruptavg, het, P_het, wave_dom_lsa, e, Tmp, N_reals,strhet)

    for realization= 1:N_reals
           t_load= tic;
            filename = [strhet,'T_rupt',num2str(t_ruptavg),'_Lf_',num2str(L_flat),'_deltaX_',num2str(deltaX),'_c_',num2str(c), '_Tmp_', num2str(Tmp),'_P_het_', num2str(P_het), '_e_', num2str(e),'rzn',num2str(realization),'.mat'];
           mk = strcat(strhet,'_Lf_',num2str(L_flat),'_deltaX_',num2str(deltaX),'_c_',num2str(c), '_Tmp_', num2str(Tmp),'_P_het_', num2str(P_het), '_e_', num2str(e));           %mk = strcat(strhet,'T_rupt',num2str(t_ruptavg),'_Lf_',num2str(L_flat),'_N_',num2str(N), '_Tmp_', num2str(Tmp),'P_het', num2str(P_het), 'e', num2str(e))
           str1 = strcat('.\',mk,'\',filename);
           load(str1);
           tElapsed_load=toc(t_load);
           fprintf('\nTime taken to complete flatFilms solver: %d min %f s\n',floor(tElapsed_load/60),mod(tElapsed_load,60));
%            ho=ones(size(x))+0.001*sin(6.*(x-(L_flat/2)).^2);  % Defining the intitial condition once more
           h_adjusted = N+5;
%            ho = [ho(end-1), ho(end), ho, ho(1), ho(2)]';   % add two ghost points each side,
%            Eo = deltaX.*sum(1.5.*((ho(4:h_adjusted-1,1)-ho(2:h_adjusted-3,1))./(2.*deltaX)).^2+ (2./(ho(3:h_adjusted-2,1)).^2) - (0.1./(3.*(ho(3:h_adjusted-2,1)).^2)));
           %c(:,all(c==0)) = [];                        
           h_final = [h_save((3:end-2),:)];
           %h_final = [h((3:end-2),:)];
%          E = [Eo,E];
           q = size(h_final);
           q=q(2);
           
           t=[0:deltaT:endTime];
           E_main(:,realization) = E(:);
           E_pi_main(:,realization) = E_pi(:);
           E_st_main(:,realization) = E_st(:);
           %q = max(size(t));
        %% function call
        if Tmp == 0
            %makeAnimation_det(animationSkip,q,c1,x,L_flat,t_new,het,P_het,e)
            %makeAnimation_det_fft(animationSkip,q,c1,x,L_flat,t_new,N,P_het,wave_dom,t)
            [S_evolution(:,:,realization)]=makeAnimation_det_comparison(animationSkip,q,h_final,x,L_flat,t,het,P_het,wave_dom_lsa,e,N,realization)
            %energy_evol(h_final,E,E_st,E_pi,x,L_flat,t,het,P_het,wave_dom_lsa,e,N,realization,q)
            %omega_max_sim = k_dom_sim_calc_save(animationSkip,q,h_final,x,deltaX,c,L_flat,t,het,P_het,wave_dom_lsa,N,realization,lin_index,k_dom_sim,e,Tmp)
            %interface_1(animationSkip,q,h,x,L_flat,deltaX,t,het,P_het,wave_dom_lsa,e,N,realization,h_adjusted,t_rupt,c)
            %disp_rel(h_final,L_flat,N,lin_index,t,q,animationSkip)
        else
            makeAnimation_det_comparison(animationSkip,q,h_final,x,L_flat,t,het,P_het,wave_dom_lsa,e,N,realization)
            %energy_evol(h_final,E,E_st,E_pi,x,L_flat,t,het,P_het,wave_dom_lsa,e,N,realization,q)
            %omega_max_sim(1,realization) = k_dom_sim_calc_save(animationSkip,q,h_final,x,deltaX,c,L_flat,t,het,P_het,wave_dom_lsa,N,realization,lin_index,k_dom_sime,Tmp)
        end
    
        
    end
    figure
    shadedErrorBar(t,E_main,{@mean,@std},'lineprops',{'r-o','markerfacecolor','r'});
    %save(filename)
    %savefig('Energy_plotfig')
%% If the function mkdir_move_postprocessing is not used    
%  mk = strcat(strhet,'T_rupt',num2str(t_ruptavg),'_Lf_',num2str(L_flat),'_N_',num2str(N), '_Tmp_', num2str(Tmp),'P_het', num2str(P_het), 'e', num2str(e));  % name your video save folder
%  mk2 = mkdir(mk);                                  % make its directory
%       movefile('Evolution_subplot*',mk)                 % move all the video files to that directory
%       movefile('filmEvolution*',mk)                 % move all the video files to that directory
%       movefile('Fourier_evolution*',mk)                 % move all the video files to that directory
%       movefile('Harmonic_evolution*',mk)                 % move all the image files to that directory
%       movefile('Final_Sk_snapshot.fig*',mk)                 % move all the image files to that directory
%       movefile('Maximum_energy_frequency*',mk)                 % move all the image files to that directory
%       
    %%
end

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

function [S_evolution]=makeAnimation_det_comparison(animationSkip,q,h_final,x,L_flat,t,het,P_het,wave_dom_lsa,e,N,realization)
%% to specify animation conditions

v = VideoWriter(strcat('Evolution_subplot',num2str(realization)));
v.FrameRate = 5;  % Default 30
v.Quality = 100;    % Default 75
open(v)
hfig = figure;
o= hfig.WindowState;
hfig.WindowState='maximize';
%l=length(c(:,1));
% x=linspace(0,L_flat,l);
% N1=length(x)
%% initializing indexes for harmonic tracing
anim_index=[1:animationSkip:q];
t_anim = t(anim_index); % time for animation
S_evolution= zeros(5,max(size(t_anim)));  % Initializing structure factor for harmonic evolution
iter=1;  %index for iteration
harmonic = [0.5,1,2,3,4];  % Will account for the first 4 harmonics
        %generating indexes for harmonics
    if het==0              
        index_het(:) = (round((L_flat/wave_dom_lsa)+harmonic(:))+1);  % +1 to account for the first element
    else
        index_het(:) = (round((L_flat/P_het).*harmonic(:))+1);   % +1 to account for the first element
    end
%%    
S_max_index = zeros(1,max(size(t_anim)));
rupt=0;   %  rupture index 
        
for i = 1:animationSkip:q
    Y = h_final(:,i);
    subplot(2,1,1)
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
        plot(x',e.*cos((2*pi*x/P_het)-pi),'LineStyle','-','color','red','LineWidth',2.5)
        xlim([0 L_flat])
        ylim([0 10])
%        legend(t_new(i),'wettability')
        legend(strcat('t=',num2str(t(i))))
        hold off
    else
         legend(strcat('t=',num2str(t(i))))
         strhet='homogeneous';
    end

    [S,f] = CalculateSk_f(Y,L_flat,N);
    subplot(2,1,2)
%    semilogy((S),'-','linewidth',2);
    loglog(f(2:end-1),(S(2:end-1)),'*','linewidth',2);
 %   semilogy(f,(S),'-','linewidth',2);
    xlabel('k [-]','Fontsize',16)
    ylabel('(S_k)[-]','Fontsize',16)
    xlim([0 10e2]);
    ylim([10e-30 10e3]);
    xlabel('k [-]','Fontsize',16)
    ylabel('Sk [-]','Fontsize',16)
    %    ax.XMinorGrid = 'on'
    set(gca,'FontSize',18)
    legend(strcat('t=',num2str(t(i))))   
    
    if rupt == 0 && min(Y) <= 0.11
        savefig(strcat('Rupture_instant',num2str(realization)))
        saveas(hfig,strcat('Rupture_instant',num2str(realization),'.jpg'))
        rupt=1;
    end
    
    B(i) = getframe(hfig);        % Everytime getframe captures the entire 
    % figure properties correspoding to hfig and stores it in M(i) which
    % can be written in a video afterwards
    [S_max,S_max_index(1,iter)] = max(S);  % Assign the indexes of max values
    
    
    %% Displaying and plotting the stripe period harmonics
        S_evolution(:,iter)= S(index_het); % Recognizing the harmonics
 %%   
    iter= iter+1;
end
B = B(~cellfun(@isempty,{B.cdata}));
writeVideo(v,B)
close(v)
 
savefig(strcat('Final_Sk_snapshot',num2str(realization)))

maxfreq_fig=figure;
plot(t_anim(2:end),f(S_max_index(2:end)))
xlabel('time')
ylabel('frequency with max Structure factor')
title('Frequency with max fourier energy evolution')
savefig(strcat('Maximum_energy_frequency',num2str(realization)))
saveas(maxfreq_fig,strcat('Maximum_energy_frequency',num2str(realization),'.jpg'))

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

function energy_evol(h_final,E,E_st,E_pi,x,L_flat,t,het,P_het,wave_dom_lsa,e,N,realization,q)

    energyplot=figure;
    E_tot = subplot(2,2,[3,4])
    plot(t(1:q),E)
    xlabel('time of evolution')
    ylabel('E-E_0')
    title('Energy evolution')
    E_piplot = subplot(2,1,2);
    plot(t(1:q),E_st)
    xlabel('time of evolution')
    ylabel('E_st')
    title('Surface tension energy evolution')
    E_piplot = subplot(2,1,1);
    plot(t(1:q),E_pi)
    xlabel('time of evolution')
    ylabel('E_pi')
    title('Disjoining pressure Energy contribution')
    

end

