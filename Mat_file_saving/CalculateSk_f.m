function [S,f] = CalculateSk_f(Y,L_flat,N)

    Y_diff=Y-mean(Y);
    %% trying to anlayse in fourier domain
    hk=fft(Y_diff);  % N+1 values
    P2=abs(hk/(L_flat));
    P1 = P2(1:N/2+1); %since the middle value would be N/2 +1
    P1(2:end-1) = 2*P1(2:end-1) ;  % as we have neglected half the energy magnitude
    f =(2*pi.*(linspace(0,(N)/2,(N+2)/2)./L_flat))' ;
    S=P1.^2;

end
