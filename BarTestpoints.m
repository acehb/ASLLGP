% generate the 200 test points for Example 1 and compute the true EQL at those points
clear;
load('BarOutput.mat');
Nx_e=10;[x_e,w_noise]=lgwt(Nx_e,0,1);
% Gauss-Legendre quadrature nodes and weights
w=zeros(1,Nx_e);
for i=1:Nx_e
     fx2 =1;
        w(i)=w_noise(i)*fx2;
 end
nt=200;
xctest=linspace(0,1,nt)';qtest=zeros(nt,1);

parfor i=1:nt
    qtest(i)=BarTrueEQL(xctest(i),stress,mass,w,x_e);
    
end
save('200 test points for the bar example.mat');
