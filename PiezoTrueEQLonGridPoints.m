% compute the true EQL at all x_c in {0,0.005,...­,1}^2 for Example 2 and Example 3
clear all;
[ri_noise,ci_noise]=lgwt(7,0,1);Nx_e=length(ri_noise)*length(ri_noise); 
w_noise=zeros(1,Nx_e);

temp1=repmat(ri_noise,1,length(ri_noise));

x_e=[reshape(temp1,[],1),reshape(temp1',[],1)];
% （x_c ∈ {0,0.005,...,1}^2）
nx=201;
[xc1,xc2]=meshgrid(0:0.005:1,0:0.005:1);
x_c_grid=[xc1(:),xc2(:)];
total_grid_points=nx*nx;
ydtemp_all=zeros(total_grid_points,Nx_e);
ymtemp_all=zeros(total_grid_points,Nx_e);

for idx=1:total_grid_points
    xc1_val=x_c_grid(idx,1);
    xc2_val=x_c_grid(idx,2);
    Lent=100e-3 + xc1_val*100e-3;
    Heit=0.6e-3 + xc2_val*0.6e-3;
    d31t=1.5e-11 + x_e(:,1)*1.5e-11;
    vt=90 + x_e(:,2)*20;
    ydtemp=zeros(Nx_e,1);
    ymtemp=zeros(Nx_e,1);

    parfor k=1:Nx_e
        [ydtemp(k),ymtemp(k)]=PiezoelectricActuator(Lent,Heit,d31t(k),vt(k),2);
    end

    ydtemp_all(idx,:)=ydtemp';
    ymtemp_all(idx,:)=ymtemp';
end
%save('Piezo True EQL on grid points yd ym.mat');
yd_grid=reshape(ydtemp_all', [Nx_e, nx, nx]);
yd_grid=permute(yd_grid, [2,3,1]); 
ym_grid=reshape(ymtemp_all', [Nx_e, nx, nx]);
ym_grid=permute(ym_grid, [2,3,1]); 
k=1;
betap1=5;betap2=5;
for i=1:length(ri_noise)
    for j=1:length(ri_noise)
        fx3 =betapdf(ri_noise(j),betap1,betap2);
        fx4 =betapdf(ri_noise(i),betap1,betap2);
        w_noise(k)=ci_noise(j)*ci_noise(i)*fx3*fx4;
        k=k+1;
    end

end
temp1=repmat(ri_noise,1,length(ri_noise));
x_e=[reshape(temp1,[],1),reshape(temp1',[],1)];
qgrid=zeros(nx,nx);maxl=0;
for ii=1:nx
    for jj=1:nx    
ydtemp = permute(yd_grid(ii,jj,:), [3,2,1]);
ymtemp = permute(ym_grid(ii,jj,:), [3,2,1]);
        losstemps=zeros(Nx_e,1);
        losstemp=zeros(Nx_e,1);
        losstempm=zeros(Nx_e,1);
        for k=1:Nx_e
            [losstemp(k),losstemps(k),losstempm(k)]=lossPiezo(ydtemp(k),ymtemp(k));
        if losstemp(k)>maxl
            maxl=losstemp(k);
        end
        end
        loss_1=w_noise*losstemps;qgrid(ii,jj)=loss_1+losstempm(1);
    end
end
b=round(maxl/100,3,'significant');
qgridBounded=zeros(nx,nx);
parfor ii=1:nx
    for jj=1:nx    
ydtemp = permute(yd_grid(ii,jj,:), [3,2,1]);
ymtemp = permute(ym_grid(ii,jj,:), [3,2,1]);
       
       losstempsBounded=zeros(Nx_e,1);
        losstempBounded=zeros(Nx_e,1);
        losstempmBounded=zeros(Nx_e,1);
        for k=1:Nx_e
            [losstempBounded(k),losstempsBounded(k),losstempmBounded(k)]=lossPiezo(ydtemp(k),ymtemp(k),b);
        end
        qgridBounded(ii,jj)=w_noise*losstempBounded;
    end
end
save('Piezo True EQL on grid points.mat');