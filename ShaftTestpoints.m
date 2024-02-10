% generate the 200 test points for Example 1 
clear;
load('ShaftOutput.mat');
Nx_e=10;[x_e,w_noise]=lgwt(Nx_e,0,1);
% Legendre-Gauss weights
sigma_noise=1/6;w=zeros(1,Nx_e);
for i=1:Nx_e
     fx2 =normpdf(x_e(i),0.5,sigma_noise)/(normcdf(1,0.5,sigma_noise)-normcdf(0,0.5,sigma_noise));
        w(i)=w_noise(i)*fx2;
 end
nt=200;
xctest=linspace(0,1,nt)';qtest=zeros(nt,1);

parfor i=1:nt
    qtest(i)=ShaftTrueEQL(xctest(i),stress,mass,w,x_e);
    
end
save('200 test points for the shaft example.mat');
