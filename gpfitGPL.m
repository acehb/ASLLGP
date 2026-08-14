function [theta,fval,beta,tao,irRes,irx,RIf,ifRIf]=gpfitGPL(xin,yin)
%Inputs xin:the design used in the experiment;yin:outputs data from the experiment
%Outputs theta:maximum integrated likelihood estimate of correlation parameters;fval:-2*(log integrated likelihood);beta:maximum likelihood estimate of the regression coefficients;%tao:maximum integrated likelihood estimate of variance parameter;irRes,irx,RIf,ifRIf are several vectors and matrices to be used by gppredict to compute point predictions and prediction intervals%loomrlengthci,loocoverage,loo95rae,loo95rlengthci:performance measures basesd on leave-one-out cross-validation 
format long g
%global z x lx

z=yin;
x=xin;
lx=size(x,2);

options=optimoptions(@patternsearch,'MaxIter',10^6,'Display','off','MaxFunctionEvaluations', 10^6);

nstart=100;nc_start=5000*lx;% 
p = sobolset(lx,'Skip',1e3,'Leap',1e2);X0 = net(p,nc_start);lb=0.01;ub=15;
can_start=lb+(ub-lb)*X0;candi=zeros(nc_start,1);
for j=1:nc_start
    [candi(j)]=omle(can_start(j,:),z,x);
end
[tempc,inds]=sort(candi,'ascend');
par=zeros(nstart,lx); fvaln=zeros(nstart,1);
parfor i=1:nstart
   [par(i,:), fvaln(i)]=patternsearch(@(par)omle(par,z,x),can_start(inds(i),:),[],[],[],[],[lb*ones(1,lx)] ,ub*ones(1,lx),[],options);
end

[fval, index]=min(fvaln);
theta=par(index,:);
nx = length(x);
rx=correlax(x,x,theta);
f=[ones(nx,1)];
irx=invandlogdet(rx);
RIf=irx*f;ifRIf=invandlogdet(f'*RIf);
beta=ifRIf*(RIf'*z);
Res=z-f*beta;pt=length(beta);
irRes=irx*Res;
tao=Res'*irRes/(nx-pt);

function mile=omle(par,z,x)
%global z x lx
lx=size(x,2);
theta=par(1:lx);
nx = length(x);
rx=correlax(x,x,theta);
f=[ones(nx,1)];
[irx, ldetx]=invandlogdet(rx);
RIX=irx*f;
XRIX=f'*RIX;
[invXRIX, ldetxrx]=invandlogdet(XRIX);
beta=invXRIX*(RIX'*z);
Res=z-f*beta;pt=length(beta);
RIRes=irx*Res;
tao=Res'*(RIRes)/(nx-pt);
mile=ldetx+(nx-pt)*log(max(tao,0))+ldetxrx;
