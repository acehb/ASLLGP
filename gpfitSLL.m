function [theta,fval,beta,tao,irRes,irx,RIf,ifRIf,transpar]=gpfitSLL(xin,yin,ubound)
%Inputs xin:the design used in the experiment;yin:outputs data from the experiment;ubound:upper bound of the loss function(omit for unbounded loss function).
%Outputs theta:maximum integrated likelihood estimate of correlation parameters;fval:-2*(log integrated likelihood);beta:maximum likelihood estimate of the regression coefficients;%tao:maximum integrated likelihood estimate of variance parameter;
%irRes,irx,RIf,ifRIf are several vectors and matrices to be used by gppredict to compute point predictions and prediction intervals;transpar:maximum integrated likelihood estimate of a and epsilon
format long g
%global z x lx Sign

z=yin;
x=xin;
lx=size(x,2);
options=optimoptions(@patternsearch,'MaxIter',10^6,'Display','off','MaxFunctionEvaluations', 10^6);
mzl=max(z);lb=0.01;ub=15;lbe=0;ube=10*max(z);
epfactor=(ube-lbe)/(ub-lb);
if(nargin==2)   
nstart=100;nc_start=5000*(lx+1);%
p = sobolset(lx+1,'Skip',1e3,'Leap',1e2);X0 = net(p,nc_start);
can_start=[lb+(ub-lb)*X0]; 
    candi=zeros(nc_start,1);
    Sign=1;
parfor j=1:nc_start
    [candi(j)]=omle(can_start(j,:),z, x, Sign, epfactor, lbe, lb);
end
[tempc,inds]=sort(candi,'ascend'); 
    par=zeros(nstart,lx+1); fvaln=zeros(nstart,1);
    parfor i=1:nstart
        [par(i,:), fvaln(i)]=patternsearch(@(par)omle(par,z, x, Sign,epfactor, lbe, lb),can_start(inds(i),:),[],[],[],[],[lb*ones(1,lx+1)],[ub*ones(1,lx+1)],[],options);
    end

    [fval, index]=min(fvaln);
    paropt=par(index,:);
    theta=paropt(1:lx);
    transpar=[(paropt(lx+1)-lb)*epfactor+lbe 1];
    
else

lbe2=ubound;ube2=10*ubound;epfactor2=(ube2-lbe2)/(ub-lb); 
    nstart=200;nc_start=2*5000*(lx+1);% 
p= sobolset(lx+1,'Skip',1e3,'Leap',1e2);X0 = net(p,nc_start/2);
can_start=[lb+(ub-lb)*X0;lb+(ub-lb)*X0]; 
    candi=zeros(nc_start,1);
for j=1:nc_start/2
    Sign=1;
    [candi(j)]=omle(can_start(j,:),z, x, Sign, epfactor, lbe, lb);
end
for j=nc_start/2+1:nc_start
    Sign=-1;
    [candi(j)]=omle(can_start(j,:),z, x, Sign,epfactor2, lbe2, lb);
end
[tempc0,inds0]=sort(candi(1:nc_start/2),'ascend');
[tempc1,inds1]=sort(candi(nc_start/2+1:nc_start),'ascend');
par=zeros(nstart,lx+1); fvaln=zeros(nstart,1);
    parfor i=1:nstart
        if(i<=nstart/2)
            Sign=1;
            [par(i,:), fvaln(i)]=patternsearch(@(par)omle(par,z, x, Sign, epfactor, lbe, lb),can_start(inds0(i),:),[],[],[],[],[lb*ones(1,lx+1) ],[ub*ones(1,lx+1)],[],options);
        else
            Sign=-1;
            [par(i,:), fvaln(i)]=patternsearch(@(par)omle(par,z, x, Sign, epfactor2, lbe2, lb),can_start(inds1(i-nstart/2)+nc_start/2,:),[],[],[],[],[lb*ones(1,lx+1)],[ub*ones(1,lx+1)],[],options);
        end
        
    end
    [fval, index]=min(fvaln);
    paropt=par(index,:);
    theta=paropt(1:lx);
    if(index<=nstart/2)
        transpar=[(paropt(lx+1)-lb)*epfactor+lbe 1];
    else
        transpar=[(paropt(lx+1)-lb)*epfactor2+lbe2 -1];
    end
end
nx = length(x);
rx=correlax(x,x,theta);
f=[ones(nx,1)];
irx=invandlogdet(rx);
RIf=irx*f;ifRIf=invandlogdet(f'*RIf);
zl2=log(transpar(2)*(z)+transpar(1));
beta=ifRIf*(RIf'*zl2);
Res=zl2-f*beta;pt=length(beta);
irRes=irx*Res;
tao=Res'*irRes/(nx-pt);

function mle=omle(par,z, x, Sign,epfac, lbeps, lb)
%global z x lx Sign
lx=size(x,2);
theta=par(1:lx);
shift=(par((lx+1))-lb)*epfac+lbeps;
nx = length(x);
rx=correlax(x,x,theta);
f=[ones(nx,1)];
[irx, ldetx]=invandlogdet(rx);
RIX=irx*f;
XRIX=f'*RIX;
[invXRIX, ldetxrx]=invandlogdet(XRIX);
z2=log(Sign*(z)+shift);
beta=invXRIX*(RIX'*z2);
Res=z2-f*beta;pt=length(beta);
RIRes=irx*Res;
taol=Res'*(RIRes)/(nx-pt);
mle=ldetx+(nx-pt)*log(max(taol,0))+ldetxrx+2*sum(z2);
