function [theta,fval,rou,beta,tao2,irRes,irx]=gpfit2level(xin,hlin,yin)
%Inputs xin:the design used in the high-fidelity experiment;yin:outputs data from the high-fidelity experiment
%Outputs theta:maximum likelihood estimate of correlation parameters;fval:-2*(log likelihood);beta:maximum likelihood estimate of the regression coefficients;tao2:maximum likelihood estimate of variance parameter;irRes,irx are several vectors and matrices to be used by gppredict to compute point predictions and prediction intervals
format long g
% global G x hl

G=yin;
x=xin;
hl=hlin;
d=size(x,2);
nstart=100;nc_start=5000*d;% 
p = sobolset(d,'Skip',1e3,'Leap',1e2);X0 = net(p,nc_start);lb=0.01;ub=25;
options=optimoptions(@patternsearch,'MaxIter',10^6,'Display','off');
can_start=lb+(ub-lb)*X0;candi=zeros(nc_start,1);
for j=1:nc_start
    [candi(j)]=omle(can_start(j,:),G,x,hl);
end
[tempc,inds]=sort(candi,'ascend');
par=zeros(nstart,d); fvaln=zeros(nstart,1);
parfor i=1:nstart
    [par(i,:), fvaln(i)]=patternsearch(@(theta)omle(theta,G,x,hl),can_start(inds(i),:),[],[],[],[],lb*ones(1,d) ,ub*ones(1,d),[],options);
end
[fval, index]=min(fvaln);
theta=par(index,:);
m = size(x,1);
rx=correlax(x,x,theta);

H=[hl,ones(m,1)];
irx=invandlogdet(rx);
irH=irx*H;ivHirH=invandlogdet(H'*irH);
gamma=ivHirH*(irH'*G);
Res=G-H*gamma;p=length(gamma);
irRes=irx*Res;
rou=gamma(1);
beta=gamma(2:p);
% tao2=Res'*irRes/(n-p);
tao2=Res'*irRes/m;

function mle=omle(theta,G,x,hl)
% global G x hl
m =  size(x,1);
rx=correlax(x,x,theta);
H=[hl,ones(m,1)];
[irx, ldetrx]=invandlogdet(rx);
irH=irx*H;[ivHirH,ldetHirH]=invandlogdet(H'*irH);
gamma=ivHirH*(irH'*G);
Res=G-H*gamma;p=length(gamma);
irRes=irx*Res;
tao2=Res'*irRes/m;
% tao2=Res'*irRes/(n-p);mle=ldetrx+(n-p)*log(max(tao2,0))+ldetHirH;
mle=ldetrx+m*log(max(tao2,0));
