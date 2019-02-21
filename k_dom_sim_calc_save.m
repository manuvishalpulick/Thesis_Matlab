function  omega_max_sim = k_dom_sim_calc_save(animationSkip,q,h_final,x,deltaX,c,L_flat,t,het,P_het,wav_dom_lsa,N,realization, lin_index,k_dom_sim,e,Tmp)

find = 1;
lin_index_1 =0;
%% initializing indexes for wave structure factor tracing
anim_index=[1:animationSkip:q];
t_anim = t(anim_index); % time for animation
S_evolution= zeros(1,max(size(t_anim)));  % Initializing structure factor for harmonic evolution
S_max_index = zeros(1,max(size(t_anim)));
k_dom_index=1;   %  rupture index 
h_flat = ones(max(size(x)),1);    % Flat film initialization for linear regime determination
wave_dom_sim = 2*pi./k_dom_sim;
index_dom = (round((L_flat/wave_dom_sim))+1);  % +1 to account for the first element            
f =(2*pi.*(linspace(0,(N)/2,(N+2)/2)./L_flat))' ;                                     
for iter = 1:1:q
        Y = h_final(:,iter);
        [S,f] = CalculateSk_f(Y,L_flat,N)
 %% Displaying and plotting the stripe period harmonics
            
            S_evolution(1,iter)= S(index_dom); % Recognizing the harmonics
            if find == 1
                lin_index_1 =0;
                h_min = min(h_final(:,iter));
                h_max = max(h_final(:,iter));
                if (h_min <= 0.95 || h_max >= 1.05) && lin_index_1 == 0
                    lin_index_1 = iter;
                    lin_time = t(iter);
                    [S_max,S_max_index(1,iter)] = max(S);  % Assign the indexes of max values
                    k_dom_index = S_max_index(1,iter);
                end
            else    
                if lin_index == iter 

                     [S_max,S_max_index(1,iter)] = max(S);  % Assign the indexes of max values
                     k_dom_index = S_max_index(1,iter);
                end
            end
    iter= iter+1;
end


      
    wave_num_fig=figure;
    plot(t(1:lin_index),log(sqrt(S_evolution(1,1:lin_index))),'linewidth',2)
    xlabel('time')
    ylabel('Dominant wave number energy in the linear')
    hold on
    p = polyfit(t(1:lin_index),log(sqrt(S_evolution(1,(1:lin_index)))),1);
    plot(t(1:lin_index),(polyval(p,t(1:lin_index))))
    legend('structure factor evolution','linear fit','Location','southeast')
    hold off
    omega_max_sim = p(1);
    str2 = strcat('wave_num_fig','_Lf_',num2str(L_flat),'_deltaX_',num2str(deltaX),'_c_',num2str(c), '_Tmp_', num2str(Tmp),'_P_het_', num2str(P_het), '_e_', num2str(e),'rzn',num2str(realization),'.fig');
    savefig(str2)
    
    %growth_rate_error(omega_max_sim, k_dom_sim ,h_flat ,L_flat ,h_final ,t ,q)
   

end