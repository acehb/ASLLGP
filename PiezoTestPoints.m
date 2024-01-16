clear;
tic
[ri_noise,ci_noise]=lgwt(6,0,1);
w_noise=zeros(1,length(ri_noise)*length(ri_noise));
k=1;
sigma_noise=1/6;
for i=1:length(ri_noise)
    for j=1:length(ri_noise)
        fx3 =normpdf(ri_noise(j),0.5,sigma_noise)/(normcdf(1,0.5,sigma_noise)-normcdf(0,0.5,sigma_noise));
        fx4 =normpdf(ri_noise(i),0.5,sigma_noise)/(normcdf(1,0.5,sigma_noise)-normcdf(0,0.5,sigma_noise));
        w_noise(k)=ci_noise(j)*ci_noise(i)*fx3*fx4;
        k=k+1;
    end
    
end
temp1=repmat(ri_noise,1,length(ri_noise));
x_e=[reshape(temp1,[],1),reshape(temp1',[],1)];
nt=400;Nx_e=k-1;
xctest=lhsdesign(nt,2);
yd_test=zeros(Nx_e,nt);
ym_test=zeros(Nx_e,nt);
Lent=100e-3+xctest(:,1)*100e-3;
Heit=0.6e-3+xctest(:,2)*0.6e-3;
d31t=1.5e-11+x_e(:,1)*1.5e-11;
vt=90+x_e(:,2)*20;
tic
for i=1:nt
    parfor j=1:Nx_e
        
        [yd_test(j,i),ym_test(j,i)]=PiezoelectricActuator(Lent(i),Heit(i),d31t(j),vt(j),2);
        
    end
end
toc

[xc1,xc2]=meshgrid(0:0.005:1,0:0.005:1);nx=201;
Lent=100e-3+xc1*100e-3;
Heit=0.6e-3+xc2*0.6e-3;
d31t=1.5e-11+x_e(:,1)*1.5e-11;
vt=90+x_e(:,2)*20;
yd_grid=zeros(nx,nx,Nx_e);
ym_grid=zeros(nx,nx,Nx_e);
tic
for ii=1:nx
    for jj=1:nx
        parfor kk=1:Nx_e
            [yd_grid(ii,jj,kk),ym_grid(ii,jj,kk)]=PiezoelectricActuator(Lent(ii,jj),Heit(ii,jj),d31t(kk),vt(kk),2);
            
        end
    end
end
toc
save('400 test points and 2D grid points.mat');