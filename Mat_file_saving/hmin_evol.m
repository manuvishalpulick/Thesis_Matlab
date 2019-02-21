function hmin_evol(h_min,t,realization,q)
figure
plot(t,h_min,'b')
xlabel('time','Fontsize',10)
ylabel('h_m_i_n','Fontsize',10)
title('Min. height vs time','Fontsize',10)
savefig(strcat('h_minfig_rzn',num2str(realization)))

end