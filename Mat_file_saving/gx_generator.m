function gx = gx_generator(N,L,x)

for i = 1:N+1;
    for k = 1:(N+1)                 
        m = k - (N/2 + 1);                            
        if m < 0
            gx(k,i) = sqrt(2/L)*sin((2*pi*m*x(i))/L); % to account for negative indices
        elseif m == 0
            gx(k,i) = sqrt(1/L);                      % the zeroth index
        else
            gx(k,i) = sqrt(2/L)*cos((2*pi*m*x(i))/L); % the positive indices (symmetric to the negative ones)
        end
    end
end

%% possibility to vectorize the code, not needed because it is executed just once

% gx= zeros(N);
% K = repmat( ((N-1)/2:-1:-(N-1)/2)',1,N );
% X = repmat( linspace(0,L,N), N,1);
% amp = sqrt(2/L);
% pi2L = pi*2/L;
% 
% gx(1:(N-1)/2, 1:N) = amp*cos(pi2L*K(1:(N-1)/2, 1:N).*X(1:(N-1)/2, 1:N));
% gx((N+1)/2, 1:N) = sqrt(1/L);
% gx((N+3)/2:end, 1:N) = amp*sin( pi2L*K((N+3)/2:end, 1:N).*X((N+3)/2:end,1:N));

%% possibility to add spatial correlations

% lc = 1;
% alp = (L/(2*lc))^2;
% for i = 1:length(x);
%     for k = 1:((N/2))                 % k going from 1 to N/2
%         g1(k,i) = -sqrt(2/L)*sin((2*pi*k*x(i))/L); % to account for negative indices
%         lmb(k) = besseli(k,alp)./besseli(0,alp); % since lambda{k} is expressed as the ratio of two modified bessel functions
%         gx1(k,i) = g1(k,i).*lmb(k);  %lambda{k} g(k}
%     end
%     gx2 = flipud(gx1);  % need to convert gx1 into a column to use flipud
%     gx0(i) = sqrt(1/L);     % for 0th mode
%     for k = 1:N/2      % now actual positive k's
%         g2(k,i) = sqrt(2/L)*cos((2*pi*(k)*x(i))/L); % cosine implementation
%         lmv(k) = besseli((k),alp)./besseli(0,alp); 
%         gx3(k,i) = g2(k,i)*lmv(k); 
%     end
%     gx4 = flipud(gx3);
%     gx = [gx2; gx0; gx4];
% end

end
