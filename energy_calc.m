function [E,E_pi,E_st]=energy_calc(h_save, x, L_flat, deltaX, c, deltaT, N, q, P_het, e, Tmp, realization, endTime)

animationSkip = 1;
h_adjusted = N+5;
N_nodes_het = round(N/(L_flat/P_het));  % No of nodes in one heterogeneous stripe
if e == 0
    m=0;
else
    m = 2.*pi./N_nodes_het; 
end
iter =1;
for evol_index = 1:animationSkip:q
    h = h_save(:,iter);
    h=h';
    hetero(:) = (1+e*cos(m.*(1:1:N+1))); % change due to heterogeneity
    %% calculation of surface tension and disjoining pressure contributions
    dE_st = (1.5.*((h(4:end-1)-h(2:end-3))./(2.*deltaX)).^2); % surface tension contribution
    E_st(iter) = trapz(x(1:end),dE_st);
    % Disjoining potential = - integral(disjoining pressure)
    datt_pot=(1./(2.*(h(3:end-2)).^2)); % attraction potential
    drep_pot= (0.1./(3.*(h(3:end-2)).^3)); % repulsion potential
    dE_pi = -hetero(:).*(datt_pot(:) - drep_pot(:));
    E_pi(iter) = trapz(x(1:end),dE_pi);
    % Total energy
    E(iter) = E_st(iter) + E_pi(iter);
    iter = iter+1;
end
t=[0:deltaT:endTime];
energy_evol_det(animationSkip,E,E_pi,E_st,t,realization,q)

end

%% Function for plotting energy evolution with time with contributions of surface tension
%and disjoining pressures to the same

function energy_evol_det(animationSkip,E,E_pi,E_st,t,realization,q)
    
    energy_piplot=subplot(2,2,1);
    plot(t(1:animationSkip:q),E_pi,'b','LineWidth',3)
    xlabel('time of evolution')
    ylabel('E_pi')
    title('Disjoining potential')
    
    energy_stplot=subplot(2,2,3);
    plot(t(1:animationSkip:q),E_st,'g','LineWidth',3)
    xlabel('time of evolution')
    ylabel('$E_{st}$')
    title('Surface tension contribution')
    
    energyplot=subplot(2,2,[2,4]);
    plot(t(1:animationSkip:q),E,'k','LineWidth',3)
    xlabel('time of evolution')
    ylabel('E')
    title('Energy evolution')
    hold on
    savefig(strcat('Energy_evolutionfig_rzn',num2str(realization)))
    

end