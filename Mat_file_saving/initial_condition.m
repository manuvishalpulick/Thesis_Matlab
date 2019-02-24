function h0 = initial_condition(k1,k2,x)
    h0 = ones(size(x))+0.001*sin(k1.*x)+0.001*sin(k2.*x);                   %0.001*sin(6.*(x-30).^2);
end