clear
tic
m=16;n=32;d=4;
ndesign=100;
xlall=zeros(ndesign,n*d);
yhd_all=zeros(m,ndesign);yld_all=zeros(n,ndesign);
ym_all=zeros(n,ndesign);
parfor k=1:ndesign
xltemp=NestedLHD(n,m,d);xtemp=xltemp(1:m,:);
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

end
toc
save('100 designs for the piezoelectric actuator example.mat');
