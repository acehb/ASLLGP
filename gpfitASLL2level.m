function [theta,fval,rou,beta,tao2,irRes,irx,transpar]=gpfitASLL2level(xin,hlin,yin,ubound)
%Inputs xin:the design used in the high-fidelity experiment;yin:outputs data from the high-fidelity experiment;
%hl: the first m outputs data from the low-fidelity experiment after the shifted log transform with the estimates of a_l and epsilon_l;ubound:upper bound of the loss function.
%Outputs theta:maximum likelihood estimate of correlation parameters;fval:-2*(log likelihood);beta:maximum likelihood estimate of the regression coefficients;%tao2:maximum likelihood estimate of variance parameter;
%irRes,irx are several vectors and matrices to be used by gppredict to compute point predictions and prediction intervals %transpar: maximum likelihood estimate of a and epsilon;
format long g
%global G x d Sign hl

G=yin;
x=xin;
hl=hlin;
maxg=max(G);d=size(x,2);lb=0.01;ub=25;lbe=0;ube=10*maxg;
options=optimoptions(@patternsearch,'MaxIter',10^6,'Display','off');

if(nargin==3)
nstart=100;nc_start=5000*(d+1);%
p = sobolset(d+1,'Skip',1e3,'Leap',1e2);X0 = net(p,nc_start);
    can_start=[lb+(ub-lb)*X0(:,1:d),ube*X0(:,d+1)]; 
    candi=zeros(nc_start,1);
    Sign=1;
for j=1:nc_start
    [candi(j)]=omile(can_start(j,:),G, x, d, Sign, hl);
end
[tempc,inds]=sort(candi,'ascend');
    par=zeros(nstart,d+1); fvaln=zeros(nstart,1);
%     Sign=1;
    parfor i=1:nstart
        [par(i,:), fvaln(i)]=patternsearch(@(par)omile(par,G, x, d, Sign, hl),can_start(inds(i),:),[],[],[],[],[lb*ones(1,d) lbe],[ub*ones(1,d) ube],[],options);
    end
    [fval, index]=min(fvaln);
    paropt=par(index,:);
    theta=paropt(1:d);
    transpar=[paropt((d+1)) 1];
    
else
    nstart=200;nc_start=2*5000*(d+1);%
p = sobolset(d+1,'Skip',1e3,'Leap',1e2);X0 = net(p,nc_start/2);
    can_start=[lb+(ub-lb)*X0(:,1:d),ube*X0(:,d+1);lb+(ub-lb)*X0(:,1:d),ubound+9*ubound*X0(:,d+1)]; 
    candi=zeros(nc_start,1);
for j=1:nc_start/2
    [candi(j)]=omile(can_start(j,:),G, x, d, 1, hl);
end
for j=nc_start/2+1:nc_start
    [candi(j)]=omile(can_start(j,:),G, x, d, -1, hl);
end
[tempc0,inds0]=sort(candi(1:nc_start/2),'ascend');
[tempc1,inds1]=sort(candi(nc_start/2+1:nc_start),'ascend');
    par=zeros(nstart,d+1); fvaln=zeros(nstart,1);
    parfor i=1:nstart
        if(i<=nstart/2)
            Sign=1;
            [par(i,:), fvaln(i)]=patternsearch(@(par)omile(par,G, x, d, Sign, hl),can_start(inds0(i),:),[],[],[],[],[lb*ones(1,d) 0],[ub*ones(1,d) ube],[],options);
        else
            Sign=-1;
            [par(i,:), fvaln(i)]=patternsearch(@(par)omile(par,G, x, d, Sign, hl),can_start(inds1(i)+nc_start/2,:),[],[],[],[],[lb*ones(1,d) ubound],[ub*ones(1,d) 10*ubound],[],options);
        end
        
    end
    [fval, index]=min(fvaln);
    paropt=par(index,:);
    theta=paropt(1:d);
    if(index<=nstart/2)
        transpar=[paropt((d+1)) 1];
    else
        transpar=[paropt((d+1)) -1];
    end
end
m = size(x,1);
rx=correlax(x,x,theta);
vec1=[ones(m,1)];
irx=invandlogdet(rx);
H=[hl,vec1];
irH=irx*H;ivHirH=invandlogdet(H'*irH);
S=log(transpar(2)*(G)+transpar(1));
gamma=ivHirH*(irH'*S);
Res=S-H*gamma;
pt=length(gamma);
irRes=irx*Res;
rou=gamma(1);
beta=gamma(2:pt);
% tao2=Res'*irRes/(n-pt);
tao2=Res'*irRes/m;

function mle=omile(par,G, x, d, Sign, hl)
% global G x d Sign hl
theta=par(1:d);
epsilon=par((d+1));
m = size(x,1);
rx=correlax(x,x,theta);
% rx=correlax(x,x,theta)+1e-6*eye(n);
vec1=[ones(m,1)];
[irx, ldetrx]=invandlogdet(rx);
H=[hl,vec1];
irH=irx*H;
hRIg=H'*irH;
[invhRIg,ldethRIg]=invandlogdet(hRIg);
S=log(Sign*(G)+epsilon);
gamma=invhRIg*(irH'*S);
Res=S-H*gamma;pt=length(gamma);
tao2=Res'*(irx*Res)/m;
%
mle=ldetrx+m*log(max(tao2,0))+2*sum(S);
% tao2=Res'*(irx*Res)/(n-pt);mile=ldetrx+(n-pt)*log(max(tao2,0))+ldethRIg+2*sum(S);
% if isnan(mle)
if ~isreal(mle)
    error(message('NaN'));
end
