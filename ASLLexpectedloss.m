function [EQL_ASLL,lb_ASLL,ub_ASLL]=ASLLexpectedloss(x_e,w,transpar,xl,thetal,xc,betal,tao2l,irResl,irxl,x,theta,rou,beta,tao2,irRes,irx)
%Inputs x_e:quadrature nodes;w:quadrature weights times xi(x_e);xl:the design used in the low-fidelity experiment;xc:control factor at which estimation of the expected loss is concerned;
%[thetal,betal,tao2l,irResl,irxl]=the outputs of gpfitASLL1level.m;x:the design used in the high-fidelity experiment;transpar: maximum likelihood estimate of a_l and epsilon_l(nargin==10) or maximum likelihood estimate of a and epsilon(nargin==17);[theta,rou,beta,tao2,irRes,irx]=the outputs of gpfitASLL2level.m
%Outputs EQL_ASLL:posterior mean of EQL given by the ASLLGP/SLLGP model;lb_ASLL:95% Lower credible limit of EQL given by the ASLLGP/SLLGP model;ub_ASLL:95% Upper credible limit of EQL given by the ASLLGP/SLLGP model
Nx_e=size(x_e,1);
x_test=[repmat(xc,Nx_e,1),x_e];
if(nargin==10)
[mu_ASLL,c_ASLL]=gppredict(xl,thetal,x_test,betal,tao2l,irResl,irxl);
else
[mu_ASLL,c_ASLL]=gppredict(xl,thetal,x_test,betal,tao2l,irResl,irxl,x,theta,rou,beta,tao2,irRes,irx);
end
mu_ASLL=mu_ASLL+log(w');
M1=sum(exp(mu_ASLL+diag(c_ASLL)/2));
temp1=exp(mu_ASLL+diag(c_ASLL)/2);
M2=sum(sum(temp1*temp1'.*exp(c_ASLL)));

mu_lambda=2*log(M1)-0.5*log(M2);
sigma2_lambda=log(M2)-2*log(M1);
EQL_ASLL=((exp(mu_lambda+sigma2_lambda/2)-sum(w)*transpar(1))/transpar(2));
alpha=0.05;
sigma_lambda=sqrt(sigma2_lambda);
Z_1 = fsolve(@(temp2)normcdf(-2.*sigma_lambda-temp2)-normcdf(temp2)-1+alpha,1e-5,optimset('display','off'));
if(transpar(2)==1)
    lb_ASLL=((exp(mu_lambda+sigma_lambda.*Z_1)-sum(w)*transpar(1))/transpar(2));
    ub_ASLL=((exp(mu_lambda+sigma_lambda.*(-2.*sigma_lambda-Z_1))-sum(w)*transpar(1))/transpar(2));
else
    ub_ASLL=((exp(mu_lambda+sigma_lambda.*Z_1)-sum(w)*transpar(1))/transpar(2));
    lb_ASLL=((exp(mu_lambda+sigma_lambda.*(-2.*sigma_lambda-Z_1))-sum(w)*transpar(1))/transpar(2));
end


