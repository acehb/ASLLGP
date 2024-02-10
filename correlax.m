function rho=correlax(x,x_prime,theta)
%Input theta: correlation parameter
%Outputs rho: correlation matrix that gives the correlations between x and x_prime.
[m, ~] = size(x);
d=size(x,2);
if isequal(x,x_prime)
    rho=ones(m,m);
    for i = 1:d
     p=abs((x(:,i)-x(:,i)'))./theta(i);rho =rho.*(exp(-p).*(p+1));
    end
else
    n=size(x_prime,1);
    rho=ones(m,n);
    for i = 1:d
     p=abs((x(:,i)-x_prime(:,i)'))./theta(i);rho =rho.*(exp(-p).*(p+1));
    end
end
