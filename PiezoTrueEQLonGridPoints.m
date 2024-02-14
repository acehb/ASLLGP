% compute the true EQL at all x_c¡Ê{0,0.005,¡­,1}^2 for Example 2 
clear
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
Nx_e=k-1;
temp1=repmat(ri_noise,1,length(ri_noise));
x_e=[reshape(temp1,[],1),reshape(temp1',[],1)];

[xc1,xc2]=meshgrid(0:0.005:1,0:0.005:1);nx=201;
qgrid=zeros(nx,nx);

parfor ii=1:nx
    for jj=1:nx     
            [qgrid(ii,jj)]=PiezoTrueEQL([xc1(ii,jj),xc2(ii,jj)],x_e,w_noise);

    end
end

save('Piezo True EQL on grid points.mat');
