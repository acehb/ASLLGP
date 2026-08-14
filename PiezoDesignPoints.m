%Generates 500 designs and outputs at the design points for Example 2 
clear all;
m=16;n=32;d=4;
rng(5);
for replicate=1:2
ndesign=500;
xlall=zeros(ndesign,n*d);
yhd_all=zeros(m,ndesign);yld_all=zeros(n,ndesign);
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
xlall(k,:)=reshape(xltemp',1,[]);
end
% Generates 500 maximin Latin hypercube designs and the corresponding outputs
nonefi=m+round(n/3.999*0.625);
xallonefi=zeros(ndesign,nonefi*d);
yd_allonefi=zeros(nonefi,ndesign);
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

end
if replicate==1
save('Designs for the piezoelectric actuator example.mat','d','m','n','ndesign','nonefi','xallonefi','xlall');
else
save('Designs for the piezoelectric actuator example-replicate.mat','d','m','n','ndesign','nonefi','xallonefi','xlall');
end
end

for replicate=1:2
    if replicate==1
load('Designs for the piezoelectric actuator example.mat');
else
load('Designs for the piezoelectric actuator example-replicate.mat');
end
parfor k=1:ndesign
xltemp=reshape(xlall(k,:),d,[])';xtemp=xltemp(1:m,:);
xlall(k,:)=reshape(xltemp',1,[]);
Len=100e-3+xltemp(:,1)*100e-3;
Hei=0.6e-3+xltemp(:,2)*0.6e-3;
d31=1.5e-11+xltemp(:,3)*1.5e-11;
v=90+xltemp(:,4)*20;
for i=1:m
[yhd_all(i,k)]=PiezoelectricActuator(Len(i),Hei(i),d31(i),v(i),2);
end
for j=1:n
[yld_all(j,k),ym_all(j,k)]=PiezoelectricActuator(Len(j),Hei(j),d31(j),v(j),1);    
end
xtemponefi=reshape(xallonefi(k,:),d,[])';
Lenonefi=100e-3+xtemponefi(:,1)*100e-3;
Heionefi=0.6e-3+xtemponefi(:,2)*0.6e-3;
d31onefi=1.5e-11+xtemponefi(:,3)*1.5e-11;
vonefi=90+xtemponefi(:,4)*20;
for i=1:nonefi
[yd_allonefi(i,k),ym_allonefi(i,k)]=PiezoelectricActuator(Lenonefi(i),Heionefi(i),d31onefi(i),vonefi(i),2);
end
end
if replicate==1
save('Designs for the piezoelectric actuator example.mat', 'yd_allonefi','d','m','n','ndesign','nonefi','xallonefi','xlall','yhd_all','yld_all','ym_all','ym_allonefi');
else
save('Designs for the piezoelectric actuator example-replicate.mat', 'yd_allonefi','d','m','n','ndesign','nonefi','xallonefi','xlall','yhd_all','yld_all','ym_all','ym_allonefi');
end
end