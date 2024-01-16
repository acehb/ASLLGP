function [EQL_AGPL,lb_AGPL,ub_AGPL]=AGPLexpectedloss(x_e,w,xl,thetal,xc,betal,tao2l,irResl,irxl,x,theta,rou,beta,tao2,irRes,irx)
%Inputs x_e:quadrature nodes;w:quadrature weights times p(x_e);xc:control factor at which estimation of the expected loss is concerned;
%[thetal,betal,tao2l,irResl,irxl]=the outputs of gpfit1level.m;x:the design used in the high-fidelity experiment;[theta,rou,beta,tao2,irRes,irx]=the outputs of gpfit2level.m;
%Outputs EQL_AGPL:posterior mean of EQL given by the AGPL model;lb_AGPL:95% Lower credible limit of EQL given by the AGPL model;ub_AGPL:95% Upper credible limit of EQL given by the AGPL model

Nx_e=size(x_e,1);
x_test=[repmat(xc,Nx_e,1),x_e];
if(nargin==9)
[mu_AGPL,c_AGPL]=gppredict(xl,thetal,x_test,betal,tao2l,irResl,irxl);
else
[mu_AGPL,c_AGPL]=gppredict(xl,thetal,x_test,betal,tao2l,irResl,irxl,x,theta,rou,beta,tao2,irRes,irx);
end   
EQL_AGPL=w*mu_AGPL;
Lvar=w*c_AGPL*w';
lb_AGPL=EQL_AGPL-sqrt(Lvar).*norminv(0.975);
ub_AGPL=EQL_AGPL+sqrt(Lvar).*norminv(0.975);





