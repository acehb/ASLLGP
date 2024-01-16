clear;
load('output.mat');
[xcgrid,xegrid] = meshgrid(0:0.1:1);
ngrid=11;yhgrid=reshape(stress(:,2),[],ngrid);
ylgrid=reshape(stress(:,1),[],ngrid);
ymtemp=reshape(mass(:,1),[],ngrid);
%
r=0.2;
ymgrid=ymtemp(1,:)-((r-0.02*2).^2+(r-0.02*2)*0.02*4+pi*(0.02.^2))*0.2*7800;

Nx_e=10;[x_e,w_noise]=lgwt(Nx_e,0,1);
% Legendre-Gauss weights
sigma_noise=1/6;w=zeros(1,Nx_e);
for i=1:Nx_e
     fx2 =normpdf(x_e(i),0.5,sigma_noise)/(normcdf(1,0.5,sigma_noise)-normcdf(0,0.5,sigma_noise));
        w(i)=w_noise(i)*fx2;
 end
% 200 test points equally spaced over [0,1] and the interpolated values for the corresponding simulator outputs
nt=200;
xctest=linspace(0,1,nt)';

ys_test=zeros(Nx_e,nt);
ym_test=zeros(Nx_e,nt);

for i=1:nt
    for j=1:Nx_e
    ys_test(j,i)=interp2(xcgrid,xegrid,yhgrid,xctest(i),x_e(j));
    ym_test(j,i)=interp1(xcgrid(1,:),ymgrid(1,:),xctest(i));
        
    end
end 
 
save('200 testpoints.mat');
