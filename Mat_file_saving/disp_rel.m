function disp_rel(h_final,L_flat,N,lin_index,t,q,animationSkip)

%%LSA Prediction
k_dom_lsa = 0.6586;   % LSA prediction of dom wave number  
ho=1;
omega = @(k) (-ho.^3.*k.^4) + (k.^2.*(1/ho)) - (0.4/3.*k.^2.*(1/ho.^2));    %dispersion relation
f =(2*pi.*(linspace(0,(N)/2,(N+2)/2)./L_flat)) ;
f_lsa = [0:0.001:(2*pi.*(N/2)./L_flat)];
omega_lsa(:) = omega(f_lsa(:));
omega_lsa_err(:) = omega(f(2:end));
%% initializing indexes
lin=linspace(1,lin_index,5000);

for j = 2:100:5000%floor(lin_index/500):-4000:1
        max_iter = j;
        S_evolution_1= zeros(max(size(f)),max_iter);  % Initializing structure factor for harmonic evolution
        %loop_skip = lin_index./
        %S_max_index = zeros(1,max(size(t_anim)));
        lin_time = t(max_iter);
        tic
        for iter = 1:max_iter
            Y = h_final(:,iter);
            [S,f] = CalculateSk_f(Y,L_flat,N);
            S_evolution_1(:,iter)= S(:);  
            
        end
        toc
        omega_sim(:) = (log(sqrt(S_evolution_1(:,end)))-log(sqrt(S_evolution_1(:,1))))/(t(max_iter));
        
        for f_iter = 2:max(size(f))
            p = polyfit(t(1:max_iter),log(sqrt(S_evolution_1(f_iter,(1:max_iter)))),1);
            omega_sim_fit(f_iter,1) = p(1);
        end
        
           


        err_omega = omega_sim_fit(1:25)- omega_lsa_err(1:25);
        err_norm = norm(err_omega)
        
        
        v = VideoWriter(strcat('omega_evolution'));
        v.FrameRate = 3;  % Default 30
        v.Quality = 100;    % Default 75
        open(v)
        wavefig = figure(j);
        %hfig = figure;
        o= wavefig.WindowState;
        wavefig.WindowState='maximize';
       
        plot(f(2:end),omega_sim(2:end),'o-','LineWidth',2,'MarkerSize',5)
        hold on
        plot(f_lsa,omega_lsa,'-','LineWidth',3)
        plot(f(2:end),omega_sim_fit(2:end),'*-','LineWidth',2,'MarkerSize',5)
        plot(f,zeros(size(f)),'k')
        
        xlabel('k','Fontsize',16)
        ylabel('\omega','Fontsize',16)
        ylim([-20 2])
        xlim([0 30])
        title('LSA PREDICTION vs SIMULATION (\omega vs k)')
        legend('\omega from Simulation','\omega from LSA','\omega by least squares regression')
%        savefig('omegafig1.fig')
  B(j) = getframe(wavefig);        % Everytime getframe captures the entire 
    % figure properties correspoding to hfig and stores it in M(i) which
    % can be written in a video afterwards     
        if err_norm <= 1
            lin = 1;
            break
        end
end

B = B(~cellfun(@isempty,{B.cdata}));
writeVideo(v,B)
close(v) 

if lin ==0
    fprintf('Still not in the linear regime')
end
        plot(f,omega_sim,'o-')
        hold on
        plot(f_lsa,omega_lsa,'-')
        %plot(f,zeros(size(f)),'k','--')
        plot(f,omega_sim_fit,'*-')
        xlabel('k')
        ylabel('omega')
        ylim([-20 2])
        xlim([0 30])
        title('LSA PREDICTION vs SIMULATION (\omega vs k)')
        legend('\omega from Simulation','\omega from LSA')
        savefig('omegafig1.fig')
        
fprintf('The norm-2 error is %f\n',err_norm)
%err_norm 

fprintf('The linear time is %f\n',lin_time)
%lin_time

fprintf('The max value of h at linear regime is %f\n\n',max(h_final(:,max_iter)))
fprintf('The min value of h at linear regime is %f\n',min(h_final(:,max_iter)))

   

end