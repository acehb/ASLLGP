function [EQL_ALL,lb_ALL,ub_ALL]=ALLexpectedloss(x_e,w,xl,thetal,xc,betal,tao2l,irResl,irxl,x,theta,rou,beta,tao2,irRes,irx)
%Inputs x_e:quadrature nodes;w:quadrature weights times xi(x_e);xl:the design used in the low-fidelity experiment;xc:control factor at which estimation of the expected loss is concerned;
%[thetal,betal,tao2l,irResl,irxl]=the outputs of gpfit1level.m;x:the design used in the high-fidelity experiment;[theta,rou,beta,tao2,irRes,irx]=the outputs of gpfit2level.m
%Outputs EQL_ALL:posterior mean of EQL given by the ALL(nargin==16)/the lognormal loss(nargin==9) model;lb_ALL:95% Lower credible limit of EQL given by the ALL(nargin==16)/the lognormal loss(nargin==9) model;ub_ALL:95% Upper credible limit of EQL given by the ALL(nargin==16)/the lognormal loss(nargin==9) model
Nx_e=size(x_e,1);
x_test=[repmat(xc,Nx_e,1),x_e];
if(nargin==9)
[mu,sigma2]=gppredict(xl,thetal,x_test,betal,tao2l,irResl,irxl);
else
[mu,sigma2]=gppredict(xl,thetal,x_test,betal,tao2l,irResl,irxl,x,theta,rou,beta,tao2,irRes,irx);    
end
mu=mu+log(w');
M1=sum(exp(mu+diag(sigma2)/2));
temp1=exp(mu+diag(sigma2)/2);
M2=sum(sum(temp1*temp1'.*exp(sigma2)));

mu_lambda=2*log(M1)-0.5*log(M2);
sigma2_lambda=log(M2)-2*log(M1);
EQL_ALL=exp(mu_lambda+sigma2_lambda/2);
alpha=0.05;
sigma_lambda=sqrt(sigma2_lambda);
Z_1 = fsolve(@(temp2)normcdf(-2.*sigma_lambda-temp2)-normcdf(temp2)-1+alpha,1e-5,optimset('display','off'));
lb_ALL=exp(mu_lambda+sigma_lambda.*Z_1);
ub_ALL=exp(mu_lambda+sigma_lambda.*(-2.*sigma_lambda-Z_1));
