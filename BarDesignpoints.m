% Generates 100 maximin nested Latin hypercube designs and the corresponding outputs
clear
load('BarOutput.mat');
m=8;n=16;d=2;
ndesign=100;
xlall=zeros(ndesign,n*d);
yhs_all=zeros(m,ndesign);yls_all=zeros(n,ndesign);
ym_all=zeros(n,ndesign);
parfor k=1:ndesign
xltemp=NestedLHD(n,m,d);xtemp=xltemp(1:m,:);
xlall(k,:)=reshape(xltemp',1,[]);
[yhs_all(:,k),~]=BarInterp(stress,mass,xtemp,2);
[yls_all(:,k),ym_all(:,k)]=BarInterp(stress,mass,xltemp,1);
end
% Generates 100 maximin Latin hypercube designs and the corresponding outputs
nonefi=m+ceil(n/90.4*5.6);
xallonefi=zeros(ndesign,nonefi*d);
ys_allonefi=zeros(nonefi,ndesign);
ym_allonefi=zeros(nonefi,ndesign);
parfor k=1:ndesign
BestDist=0; BestDesign=[ ];
for idx=1:500
Designonefi=lhsdesign(nonefi,d,'Criterion','none');
Dist=min(pdist(Designonefi));
    if(Dist>BestDist)
        BestDist=Dist;
        BestDesign=Designonefi;
    end
end
xtemponefi=BestDesign;
xallonefi(k,:)=reshape(xtemponefi',1,[]);
[ys_allonefi(:,k),ym_allonefi(:,k)]=BarInterp(stress,mass,xtemponefi,2);

end
save('Designs for the bar example.mat');
