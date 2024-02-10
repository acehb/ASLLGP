function [mu,var,mul,varl]=gppredict(xl,thetal,x_test,betal,tao2l,irResl,irxl,x,theta,rou,beta,tao2,irRes,irx)
%Inputs xl:the design used in the low-fidelity experiment;x_test:input setting at which the prediction of EQL is desired;
%x:the design used in the high-fidelity experiment;[thetal,betal,tao2l,irResl,irxl]=the outputs of gpfit1level/gpfitASLL1level;[theta,rou,beta,tao2,irRes,irx]=the outputs of gpfit2level/gpfitASLL2level
%Outputs mu:posterior mean of the ASLLGP/ALL/AGPL emulator at x_test;var:posterior covariance of the ASLLGP/ALL/AGPL emulator at x_test;
%mul:posterior mean of the SLLGP/LL/GPL emulator at x_test;varl:posterior covariance of the SLLGP/LL/GPL emulator at x_test;
ntest=size(x_test,1);
vec1_test=[ones(ntest,1)];
rl= correlax(xl,x_test,thetal);
mul=vec1_test*betal+rl'*irResl;

varl=max(tao2l*(correlax(x_test,x_test,thetal)-rl'*irxl*rl),0);

if(nargin==7)
    mu=mul;var=varl;
    return
else
 h_test=[mul,ones(ntest,1)]; gamma=[rou;beta];
r= correlax(x,x_test,theta);
mu=h_test*gamma+r'*irRes;

var=max(rou^2*varl+tao2*(correlax(x_test,x_test,theta)-r'*irx*r),0);
end
end
