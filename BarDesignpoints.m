% Generates 500 maximin nested Latin hypercube designs and the corresponding outputs
%Inputs n:the number of runs for the low-fidelity experiment;m:the number of runs for the high-fidelity experiment;d:dimension of the inputs.
clear all;
load('BarOutput.mat');
rng(10);
for replicate=1:2
m=8;n=16;d=2;
ndesign=500;
xlall=zeros(ndesign,n*d);
yhs_all=zeros(m,ndesign);yls_all=zeros(n,ndesign);
ym_all=zeros(n,ndesign);
for k=1:ndesign
%xltemp:the maximin nested Latin hypercube design
pin=zeros(d,n);tao=zeros(d,m);
for j=1:d
    pim=randperm(m);
    t=n/m;
    for i=1:m
        tao(j,i)=randi([(pim(i)-1)*t+1,pim(i)*t],1);
    end
    pin(j,:)=[tao(j,:), setdiff(randperm(n),tao(j,:),'stable') ];
end
xltemp=(pin-rand(d,n))'./n;
%maxdist = mindist(xltemp);
[~,dist] = knnsearch(xltemp,xltemp,'k',2);maxdist = min(dist(:,2));
for i=2:500
    %x = nestedsample(n,m,d);
    pin=zeros(d,n);
    tao=zeros(d,m);
    for j=1:d
        pim=randperm(m);
        t=n/m;
        for ii=1:m
            tao(j,ii)=randi([(pim(ii)-1)*t+1,pim(ii)*t],1);
        end
        pin(j,:)=[tao(j,:), setdiff(randperm(n),tao(j,:),'stable') ];
    end
    x=(pin-rand(d,n))'./n;
    [~,dist0] = knnsearch(x,x,'k',2);newdist =min(dist0(:,2));
    %newdist = mindist(x);
    if newdist > maxdist
        xltemp = x;
        maxdist = newdist;
    end
end
xtemp=xltemp(1:m,:);
xlall(k,:)=reshape(xltemp',1,[]);
[yhs_all(:,k),~]=BarInterp(stress,mass,xtemp,2);
[yls_all(:,k),ym_all(:,k)]=BarInterp(stress,mass,xltemp,1);
end
% Generates 500 maximin Latin hypercube designs and the corresponding outputs
nonefi=m+round(n/90.5*5.4);
xallonefi=zeros(ndesign,nonefi*d);
ys_allonefi=zeros(nonefi,ndesign);
ym_allonefi=zeros(nonefi,ndesign);
for k=1:ndesign
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
if replicate==1
save('Designs for the bar example.mat', 'mass','d','m','n','ndesign','nonefi','stress','xallonefi','xlall','yhs_all','yls_all','ym_all','ym_allonefi','ys_allonefi');
else
save('Designs for the bar example-replicate.mat', 'mass','d','m','n','ndesign','nonefi','stress','xallonefi','xlall','yhs_all','yls_all','ym_all','ym_allonefi','ys_allonefi');
end
end