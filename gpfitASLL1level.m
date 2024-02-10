function [thetal,fvall,betal,tao2l,irResl,irxl,transparl]=gpfitASLL1level(xin,yin,ubound)
%Inputs xin:the design used in the low-fidelity experiment;yin:outputs data from the low-fidelity experiment;ubound:the upper bound of loss function(omit for unbounded loss function);
%Outputs theta:maximum likelihood estimate of correlation parameters;fval:-2*(log likelihood)-n;beta:maximum likelihood estimate of the regression coefficients;
%tao2:maximum likelihood estimate of variance parameter;irRes,irx are the vector and matrix to be used by gppredict to compute point predictions and prediction intervals %transparl: maximum likelihood estimate of a_l and epsilon_l;
format long g
% global G x d Sign

G=yin;
x=xin;
d=size(x,2);
maxg=max(G);d=size(x,2);lb=0.01;ub=25;lbe=0;ube=10*maxg;
options=optimoptions(@patternsearch,'MaxIter',10^6,'Display','off');
if(nargin==2)
nstart=100;nc_start=5000*(d+1);%
p = sobolset(d+1,'Skip',1e3,'Leap',1e2);X0 = net(p,nc_start);
can_start=[lb+(ub-lb)*X0(:,1:d),ube*X0(:,d+1)]; 
    candi=zeros(nc_start,1);
    Sign=1;
for j=1:nc_start
    [candi(j)]=omle(can_start(j,:),G, x, d, Sign);
end
[tempc,inds]=sort(candi,'ascend'); 
    par=zeros(nstart,d+1); fvaln=zeros(nstart,1);
    parfor i=1:nstart
        [par(i,:), fvaln(i)]=patternsearch(@(par)omle(par,G, x, d, Sign),can_start(inds(i),:),[],[],[],[],[lb*ones(1,d) lbe],[ub*ones(1,d) ube],[],options);
    end
    [fvall, index]=min(fvaln);
    paropt=par(index,:);
    thetal=paropt(1:d);
    transparl=[paropt((d+1)) 1];
    
else
 nstart=200;nc_start=2*5000*(d+1);% 
p= sobolset(d+1,'Skip',1e3,'Leap',1e2);X0 = net(p,nc_start/2);
can_start=[lb+(ub-lb)*X0(:,1:d),ube*X0(:,d+1);lb+(ub-lb)*X0(:,1:d),ubound+9*ubound*X0(:,d+1)]; 
    candi=zeros(nc_start,1);
for j=1:nc_start/2
    [candi(j)]=omle(can_start(j,:),G, x, d, 1);
end
for j=nc_start/2+1:nc_start
    [candi(j)]=omle(can_start(j,:),G, x, d, -1);
end
[tempc0,inds0]=sort(candi(1:nc_start/2),'ascend');
[tempc1,inds1]=sort(candi(nc_start/2+1:nc_start),'ascend');
par=zeros(nstart,d+1); fvaln=zeros(nstart,1);
    parfor i=1:nstart
        if(i<=nstart/2)
            Sign=1;
            [par(i,:), fvaln(i)]=patternsearch(@(par)omle(par,G, x, d, Sign),can_start(inds0(i),:),[],[],[],[],[lb*ones(1,d) lbe],[ub*ones(1,d) ube],[],options);
        else
            Sign=-1;
            [par(i,:), fvaln(i)]=patternsearch(@(par)omle(par,G, x, d, Sign),can_start(inds1(i-nc_start/2)+nc_start/2,:),[],[],[],[],[lb*ones(1,d) ubound],[ub*ones(1,d) 10*ubound],[],options);
        end
        
    end
    [fvall, index]=min(fvaln);
    paropt=par(index,:);
    thetal=paropt(1:d);
    if(index<=nstart/2)
        transparl=[paropt((d+1)) 1];
    else
        transparl=[paropt((d+1)) -1];
    end
end
n = size(x,1);
rx=correlax(x,x,thetal);
vec1=[ones(n,1)];
irxl=invandlogdet(rx);
irvec1=irxl*vec1;ivirvec1=invandlogdet(vec1'*irvec1);
S=log(transparl(2)*(G)+transparl(1));
betal=ivirvec1*(irvec1'*S);
Res=S-vec1*betal;
irResl=irxl*Res;
% tao2=Res'*irRes/(n-1);
tao2l=Res'*irResl/n;

function mle=omle(par,G, x, d, Sign)
% global G x d Sign
theta=par(1:d);
epsilon=par((d+1));
n = size(x,1);
rx=correlax(x,x,theta);
vec1=[ones(n,1)];
[irx, ldetrx]=invandlogdet(rx);
irvec1=irx*vec1;
[ivirvec1, ldet1r1]=invandlogdet(vec1'*irvec1);
S=log(Sign*(G)+epsilon);
beta=ivirvec1*(irvec1'*S);
Res=S-vec1*beta;
RIRes=irx*Res;
tao2=Res'*(RIRes)/n;
mle=ldetrx+n*log(max(tao2,0))+2*sum(S);

