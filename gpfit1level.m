function [thetal,fvall,betal,tao2l,irResl,irxl]=gpfit1level(xin,yin)
%Inputs xin:the design used in the low-fidelity experiment;yin:outputs data from the low-fidelity experiment
%Outputs thetal:maximum likelihood estimate of correlation parameters;fvall:-2*(log likelihood);betal:maximum likelihood estimate of the regression coefficients;tao2l:maximum likelihood estimate of variance parameter;irResl,irxl are the vector and matrix to be used by gppredict to compute point predictions and prediction intervals
format long g
% global G x
G=yin;
x=xin;
d=size(x,2);
nstart=100;nc_start=5000*d;% 
p = sobolset(d,'Skip',1e3,'Leap',1e2);X0 = net(p,nc_start);lb=0.01;ub=25;
options=optimoptions(@patternsearch,'MaxIter',10^6,'Display','off');
can_start=lb+(ub-lb)*X0;candi=zeros(nc_start,1);
for j=1:nc_start
    [candi(j)]=omle(can_start(j,:),G,x);
end
[tempc,inds]=sort(candi,'ascend');
par=zeros(nstart,d); fvaln=zeros(nstart,1);
parfor i=1:nstart
   [par(i,:), fvaln(i)]=patternsearch(@(thetal)omle(thetal,G,x),can_start(inds(i),:),[],[],[],[],lb*ones(1,d) ,ub*ones(1,d),[],options);
end
[fvall, index]=min(fvaln);
thetal=par(index,:);
n = size(x,1);
rx=correlax(x,x,thetal);

vec1=[ones(n,1)];
irxl=invandlogdet(rx);
irvec1=irxl*vec1;ivirvec1=invandlogdet(vec1'*irvec1);
betal=ivirvec1*(irvec1'*G);
Res=G-vec1*betal;
irResl=irxl*Res;
tao2l=Res'*irResl/n;

function mle=omle(thetal,G,x)
% global G x
n = length(x);
rx=correlax(x,x,thetal);
vec1=[ones(n,1)];
[irx, ldetrx]=invandlogdet(rx);
irvec1=irx*vec1;
[ivirvec1, ldet1r1]=invandlogdet(vec1'*irvec1);
betal=ivirvec1*(irvec1'*G);
Res=G-vec1*betal;
tao2l=Res'*(irx*Res)/n;
% tao2l=Res'*(irx*Res)/(n-1);mle=ldetrx+(n-1)*log(max(tao2l,0))+ldet1r1;
mle=ldetrx+n*log(max(tao2l,0));

