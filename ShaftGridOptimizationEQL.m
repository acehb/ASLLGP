clear
load('Designs for the shaft example.mat');
Nx_e=10;[x_e,w_noise]=lgwt(Nx_e,0,1);
% Legendre-Gauss weights
sigma_noise=1/6;w=zeros(1,Nx_e);
for i=1:Nx_e
    fx2 =normpdf(x_e(i),0.5,sigma_noise)/(normcdf(1,0.5,sigma_noise)-normcdf(0,0.5,sigma_noise));
    w(i)=w_noise(i)*fx2;
end
nx=1001;x_c=linspace(0,1,nx)';
% compute the true EQL at the grid {0,0.001,бн,1}.
qgrid=zeros(nx,1);
parfor i=1:nx
    qgrid(i)=ShaftTrueEQL(x_c(i),stress,mass,w,x_e);
end
tic

% choose one design
choose_d=12
xl=reshape(xlall(choose_d,:),d,[])';
yls=yls_all(:,choose_d);yhs=yhs_all(:,choose_d);ym=ym_all(:,choose_d);
x=xl(1:m,:);
Ll=zeros(n,1);Ll_s=zeros(n,1);Lh=zeros(m,1);Lh_s=zeros(m,1);
for i=1:n
    [Ll(i),Ll_s(i)]=lossShaft(yls(i),ym(i));
end
for i=1:m
    [Lh(i),Lh_s(i)]=lossShaft(yhs(i),ym(i));
end
%fit 3 models

[ASLLthetal,ASLLfvall,ASLLbetal,ASLLtao2l,ASLLirResl,ASLLirxl,transparl]=gpfitASLL1level(xl,Ll);
ASLLhl=log(transparl(2)*Ll(1:m)+transparl(1));
[ASLLtheta,ASLLfval,ASLLrou,ASLLbeta,ASLLtao2,ASLLirRes,ASLLirx,transpar]=gpfitASLL2level(x,ASLLhl,Lh);

[AGPLthetal,AGPLfvall,AGPLbetal,AGPLtao2l,AGPLirResl,AGPLirxl]=gpfit1level(xl,Ll);
Lhl=Ll(1:m);
[AGPLtheta,AGPLfval,AGPLrou,AGPLbeta,AGPLtao2,AGPLirRes,AGPLirx]=gpfit2level(x,Lhl,Lh);

ALLLh=log(Lh);ALLLl=log(Ll); ALLhl=ALLLl(1:m);
[ALLthetal,ALLfvall,ALLbetal,ALLtao2l,ALLirResl,ALLirxl]=gpfit1level(xl,ALLLl);
[ALLtheta,ALLfval,ALLrou,ALLbeta,ALLtao2,ALLirRes,ALLirx]=gpfit2level(x,ALLhl,ALLLh);
%predict EQL
lossLCL_ASLL=zeros(nx,1);
lossUCL_ASLL=zeros(nx,1);
EQLASLL=zeros(nx,1);
parfor i=1:nx
    
    [EQLASLL(i),lossLCL_ASLL(i),lossUCL_ASLL(i)]=ASLLexpectedloss(x_e,w,transpar,xl,ASLLthetal,x_c(i),ASLLbetal,ASLLtao2l,ASLLirResl,ASLLirxl...
        ,x,ASLLtheta,ASLLrou,ASLLbeta,ASLLtao2,ASLLirRes,ASLLirx);
    
end
lossLCL_AGPL=zeros(nx,1);
lossUCL_AGPL=zeros(nx,1);
EQLAGPL=zeros(nx,1);
parfor i=1:nx
    [EQLAGPL(i),lossLCL_AGPL(i),lossUCL_AGPL(i)]=AGPLexpectedloss(x_e,w,xl,AGPLthetal,x_c(i),AGPLbetal,AGPLtao2l,AGPLirResl,AGPLirxl,...
        x,AGPLtheta,AGPLrou,AGPLbeta,AGPLtao2,AGPLirRes,AGPLirx);
    
end

lossLCL_ALL=zeros(nx,1);
lossUCL_ALL=zeros(nx,1);
EQLALL=zeros(nx,1);
parfor i=1:nx
    [EQLALL(i),lossLCL_ALL(i),lossUCL_ALL(i)]=ALLexpectedloss(x_e,w,xl,ALLthetal,x_c(i),ALLbetal,ALLtao2l,ALLirResl,ALLirxl,...
        x,ALLtheta,ALLrou,ALLbeta,ALLtao2,ALLirRes,ALLirx);
    
end
% the estimated optimal control-factor setting and the minimizer of the posterior mean of the EQL, given by the ASLLGP, ALL and AGPL models.
[minq,mintemp]=min(qgrid);xcmin=x_c(mintemp);
[minALL,minALLtemp]=min(EQLALL);xcminALL=x_c(minALLtemp);
[minASLL,minASLLtemp]=min(EQLASLL);xcminASLL=x_c(minASLLtemp);
[minAGPL,minAGPLtemp]=min(EQLAGPL);xcminAGPL=x_c(minAGPLtemp);
Optxc=[xcmin,xcminASLL,xcminALL,xcminAGPL]
% minEQLesti=[minq,minASLL,minALL,minAGPL]
%The true EQL at the estimated optimal control-factor setting given by the ASLLGP, ALL and AGPL models.
minALLEQL=ShaftTrueEQL(xcminALL,stress,mass,w,x_e);
minASLLEQL=ShaftTrueEQL(xcminASLL,stress,mass,w,x_e);
minAGPLEQL=ShaftTrueEQL(xcminAGPL,stress,mass,w,x_e);
minEQLtrue=[minq,minASLLEQL,minALLEQL,minAGPLEQL]
%% plot of the posterior mean, and lower and upper 95% credible limits for the EQL given by the ASLLGP, ALL and AGPL models versus x_c.
transq=logscale(qgrid);transEQLALL=logscale(EQLALL);translossUCL_ALL=logscale(lossUCL_ALL);translossLCL_ALL=logscale(lossLCL_ALL);
transEQLAGPL=logscale(EQLAGPL);translossLCL_AGPL=logscale(lossLCL_AGPL);translossUCL_AGPL=logscale(lossUCL_AGPL);
transEQLASLL=logscale(EQLASLL);translossLCL_ASLL=logscale(lossLCL_ASLL);translossUCL_ASLL=logscale(lossUCL_ASLL);

choosep=[301:20:460,460:30:1001];chooseb=[301:10:460,470:20:1000];
%choosep=[61:4:92,92:6:201];chooseb=[61:2:92,94:4:200];
x_c31=x_c(choosep);transq31=transq(choosep);
x_cb=x_c(chooseb);translossUCL_ALL31=translossUCL_ALL(chooseb);translossLCL_ALL31=translossLCL_ALL(chooseb);transEQLALL31=transEQLALL(choosep);
transEQLAGPL31=transEQLAGPL(choosep);translossLCL_AGPL31=translossLCL_AGPL(chooseb);translossUCL_AGPL31=translossUCL_AGPL(chooseb);
transEQLASLL31=transEQLASLL(choosep);translossLCL_ASLL31=translossLCL_ASLL(chooseb);translossUCL_ASLL31=translossUCL_ASLL(chooseb);
ylbgpl=-0.14;yubgpl=0.3;ftick5=25;ygpl=[-0.1:0.05:0.25];
ylbll=-0.0005;yubll=0.23;yll=[0:0.04:0.2];ms=20;lw=4;fx5=0.04;
wx=0.28;wy=0.88;wdx=0.05;fy=0.085;ft=21;

figure(5)
subplot(1,3,2);subplot('position',[fx5+wdx+wx*1 fy,wx,wy]);
plot(x_c31,transq31,'-.r','LineWidth',lw);hold on;plot(x_c31,transEQLALL31,'-','Color',[.8,.3,0],'LineWidth',lw);
hold on;plot(x_cb,translossLCL_ALL31,'k:','LineWidth',lw,'MarkerSize',ms);hold on;f1=plot(x_cb,translossUCL_ALL31,'k:','LineWidth',lw,'MarkerSize',ms);set(f1,'handlevisibility','off');
hold off;title('(b)','FontSize',ft,'FontWeight','bold');
leg2=legend('True EQL',[' ',sprintf('\n'),'Posterior mean of the EQL ',sprintf('\n'),'given by the ALL model'],[ ' ',sprintf('\n'),'95% upper/lower credible',sprintf('\n'),'interval limit for the EQL',sprintf('\n'),'given by the ALL model']);
xlabel('Control factor x_c','FontSize',ft,'FontWeight','bold');
xlim([0.3,1]);xticks([0.3:0.1:1]);a2 = get(gca,'XTickLabel');set(gca,'XTickLabel',a2,'fontsize',ftick5,'fontweight','bold')
ylabel('EQL','FontSize',ft,'FontWeight','bold');
ylim([ylbll,yubll]);set(gca,'yTick',[yll]);set(gca,'YTickLabel',{'0','10','20','32','45','58'});
set(gca,'fontsize',ft,'fontweight','bold');
subplot(1,3,3);subplot('position',[fx5+wdx*2+wx*2, fy,wx,wy]);
plot(x_c31,transq31,'-.r','LineWidth',lw);hold on;plot(x_c31,transEQLAGPL31,'-','Color',[0 0.4 0],'LineWidth',lw);
hold on;plot(x_cb,translossLCL_AGPL31,'k:','LineWidth',lw,'MarkerSize',ms);hold on;plot(x_cb,translossUCL_AGPL31,'k:','LineWidth',lw,'MarkerSize',ms);
hold off;title('(c)','FontSize',ft,'FontWeight','bold');
legend('True EQL',[' ',sprintf('\n'),'Posterior mean of the EQL ',sprintf('\n'),'given by the AGPL model'],[' ',sprintf('\n'),'95% upper/lower credible',sprintf('\n'),'interval limit for the EQL',sprintf('\n'),'given by the AGPL model']);
xlabel('Control factor x_c','FontSize',ft,'FontWeight','bold');
xlim([0.3,1]);xticks([0.3:0.1:1]);a3 = get(gca,'XTickLabel');set(gca,'XTickLabel',a3,'fontsize',ftick5,'fontweight','bold')
ylabel('EQL','FontSize',ft,'FontWeight','bold');
ylim([ylbgpl,yubgpl]);set(gca,'yTick',[ygpl]);set(gca,'YTickLabel',{'-26','-12','0','12','26','41','58','78'});
set(gca,'fontsize',ft,'fontweight','bold');
subplot(1,3,1);subplot('position',[fx5 fy,wx,wy]);
plot(x_c31,transq31,'-.r','LineWidth',lw);hold on;plot(x_c31,transEQLASLL31,'b-','LineWidth',lw);hold on;
plot(x_cb,translossLCL_ASLL31,'k:','LineWidth',lw,'MarkerSize',ms);hold on;plot(x_cb,translossUCL_ASLL31,'k:','LineWidth',lw,'MarkerSize',ms);
hold off;title('(a)','FontSize',ft,'FontWeight','bold');
legend('True EQL',[' ',sprintf('\n'),'Posterior mean of the EQL',sprintf('\n'),'given by the ASLLGP model'],[' ',sprintf('\n'),'95% upper/lower credible',sprintf('\n'),'interval limit for the EQL',sprintf('\n'),'given by the ASLLGP model']);
xlabel('Control factor x_c','FontSize',ft,'FontWeight','bold');
xlim([0.3,1]);xticks([0.3:0.1:1]);a1 = get(gca,'XTickLabel');set(gca,'XTickLabel',a1,'fontsize',ftick5,'fontweight','bold')
ylabel('EQL','FontSize',ft,'FontWeight','bold');
ylim([ylbll,yubll]);set(gca,'yTick',[yll]);set(gca,'YTickLabel',{'0','10','20','32','45','58'});
set(gca,'fontsize',ft,'fontweight','bold');

function logx=logscale(x)
% modified log scale to the y-axis
C=2;
logx = sign(x).*(log10(1+abs(x)./(10^C)));
end