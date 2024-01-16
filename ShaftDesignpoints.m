clear
load('output.mat');
[xcgrid,xegrid] = meshgrid(0:0.1:1);
ngrid=11;yhgrid=reshape(stress(:,2),[],ngrid);
ylgrid=reshape(stress(:,1),[],ngrid);
ymtemp=reshape(mass(:,1),[],ngrid);
%
r=0.2;
ymgrid=ymtemp(1,:)-((r-0.02*2).^2+(r-0.02*2)*0.02*4+pi*(0.02.^2))*0.2*7800;

m=8;n=16;d=2;
ndesign=100;
xlall=zeros(ndesign,n*d);
yhs_all=zeros(m,ndesign);yls_all=zeros(n,ndesign);
ym_all=zeros(n,ndesign);
parfor k=1:ndesign
xltemp=NestedLHD(n,m,d);xtemp=xltemp(1:m,:);
xlall(k,:)=reshape(xltemp',1,[]);
yhs_all(:,k)=interp2(xcgrid,xegrid,yhgrid,xtemp(:,1),xtemp(:,2));
yls_all(:,k)=interp2(xcgrid,xegrid,ylgrid,xltemp(:,1),xltemp(:,2));
ym_all(:,k)=interp1(xcgrid(1,:),ymgrid,xltemp(:,1)); 

end

nonefi=m+ceil(n/90.4*5.6);
xallonefi=zeros(ndesign,nonefi*d);
ys_allonefi=zeros(nonefi,ndesign);
ym_allonefi=zeros(nonefi,ndesign);
for k=1:ndesign
xtemponefi=lhsdesign(nonefi,d);
xallonefi(k,:)=reshape(xtemponefi',1,[]);
     ys_allonefi(:,k)=interp2(xcgrid,xegrid,yhgrid,xtemponefi(:,1),xtemponefi(:,2));
    ym_allonefi(:,k)=interp1(xcgrid(1,:),ymgrid,xtemponefi(:,1)); 

end
save('100 designs.mat');
