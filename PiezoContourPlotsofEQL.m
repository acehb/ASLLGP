clear;
load('Designs for the piezoelectric actuator example.mat');
load('Piezo True EQL on grid points.mat');
[indextrue1,indextrue2]=find(qgrid==min(min(qgrid)));xcmin=[xc1(indextrue1,indextrue2),xc2(indextrue1,indextrue2)];

% choose one design
choose_d=86
xl=reshape(xlall(choose_d,:),d,[])';
yld=yld_all(:,choose_d);yhd=yhd_all(:,choose_d);ym=ym_all(:,choose_d);
x=xl(1:m,:);
Ll=zeros(n,1);Ll_d=zeros(n,1);Lh=zeros(m,1);Lh_d=zeros(m,1);
for i=1:n
    [Ll(i),Ll_d(i)]=lossPiezo(yld(i),ym(i));
end
for i=1:m
    [Lh(i),Lh_d(i)]=lossPiezo(yhd(i),ym(i));
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
expectedlossASLL=zeros(nx,nx);
parfor i=1:nx
    for j=1:nx
        [expectedlossASLL(i,j),~,~]=ASLLexpectedloss(x_e,w_noise,transpar,xl,ASLLthetal,[xc1(i,j),xc2(i,j)],ASLLbetal,ASLLtao2l,ASLLirResl,ASLLirxl...
            ,x,ASLLtheta,ASLLrou,ASLLbeta,ASLLtao2,ASLLirRes,ASLLirx);
    end
end
[indexASLL1,indexASLL2]=find(expectedlossASLL==min(min(expectedlossASLL)));xcminASLL=[xc1(indexASLL1,indexASLL2),xc2(indexASLL1,indexASLL2)];

expectedlossAGPL=zeros(nx,nx);
parfor i=1:nx
    for j=1:nx
        [expectedlossAGPL(i,j),~,~]=AGPLexpectedloss(x_e,w_noise,xl,AGPLthetal,[xc1(i,j),xc2(i,j)],AGPLbetal,AGPLtao2l,AGPLirResl,AGPLirxl,...
            x,AGPLtheta,AGPLrou,AGPLbeta,AGPLtao2,AGPLirRes,AGPLirx);
    end
end
[indexAGPL1,indexAGPL2]=find(expectedlossAGPL==min(min(expectedlossAGPL)));xcminAGPL=[xc1(indexAGPL1,indexAGPL2),xc2(indexAGPL1,indexAGPL2)];

expectedlossALL=zeros(nx,nx);
parfor i=1:nx
    for j=1:nx
        [ expectedlossALL(i,j),~,~]=ALLexpectedloss(x_e,w_noise,xl,ALLthetal,[xc1(i,j),xc2(i,j)],ALLbetal,ALLtao2l,ALLirResl,ALLirxl,...
            x,ALLtheta,ALLrou,ALLbeta,ALLtao2,ALLirRes,ALLirx);
    end
end
[indexALL1,indexALL2]=find(expectedlossALL==min(min(expectedlossALL)));xcminALL=[xc1(indexALL1,indexALL2),xc2(indexALL1,indexALL2)];

tEQL=min(min(qgrid))
minEQL=[min(min(expectedlossASLL)),min(min(expectedlossALL)),min(min(expectedlossAGPL))]
%% Figure 5: contour plot of the posterior mean, and lower and upper 95% credible limits for the EQL given by the ASLLGP, ALL and AGPL models versus x_c=(x_1,x_2).
range1=[1,2,4,6,8];rangeGPL=[0,1,2,4,6,8];
labs=21;lw=4;ls2=699;ls1=260;ls3=260;%ls3=240;%ls3=250;
wx=0.28;wy=0.88;wdx=0.05;fy=0.085;fx=0.05;fleg=20;ftick=21;ft=21;
figure(4)
subplot(1,3,1);subplot('position',[fx,fy,wx-0.01,wy]);
[casll1,hasll1]=contour(xc1,xc2,qgrid*1e5,range1,'--r');hold on;% [casll1,hasll1]=contour(xc1,xc2,qgrid,range1,'--r','ShowText','on');hold on;
h1 = plot(NaN, '-.r','LineWidth',lw);
[casll2,hasll2]=contour(xc1,xc2,expectedlossASLL*1e5,range1,'b','ShowText','on');
h2 = plot(NaN, '-b','LineWidth',lw);
title('(a)','FontSize',ft,'FontWeight','normal');
xlabel('Control factor x_1','FontSize',labs,'FontWeight','normal');ylabel('Control factor x_2','FontSize',labs,'FontWeight','normal');
clabel(casll2,hasll2,'FontSize',labs,'FontWeight','normal');hasll1.LineWidth = lw;hasll2.LineWidth = lw;
hasll2.LabelSpacing=ls2;
clabel(casll1,hasll1,'FontSize',labs,'FontWeight','normal');hasll1.LabelSpacing=ls1;
leg1=legend([h1 h2],['True',sprintf('\n'),'EQL{\times}10^{5}'],[' ',sprintf('\n'),'Posterior',sprintf('\n'),'mean of',sprintf('\n'), 'the EQL',sprintf('\n'), 'given by',sprintf('\n'), 'the ASLLGP',sprintf('\n'), 'model{\times}10^{5}']);
set(leg1,'FontSize',fleg,'FontWeight','bold');
xlim([0,1]);xticks([0:0.2:1]);ca1 = get(gca,'XTickLabel');set(gca,'XTickLabel',ca1,'fontsize',ftick,'fontweight','normal');
yticks([0:0.2:1]);cay1 = get(gca,'YTickLabel');set(gca,'YTickLabel',cay1,'fontsize',ftick,'fontweight','normal');
subplot(1,3,2);subplot('position',[fx+wdx+wx*1,fy,wx,wy]);
[call1,hall1]=contour(xc1,xc2,qgrid*1e5,range1,'--r','ShowText','on');hold on;
h3 = plot(NaN, '-.r','LineWidth',lw);
[call2,hall2]=contour(xc1,xc2,expectedlossALL*1e5,range1,'ShowText','on','EdgeColor',[.8,.3,0]);
h4 = plot(NaN, '-','Color',[.8,.3,0],'LineWidth',lw);
leg2=legend([h3 h4],['True',sprintf('\n'),'EQL{\times}10^{5}'],[' ',sprintf('\n'),'Posterior',sprintf('\n'),'mean of',sprintf('\n'), 'the EQL',sprintf('\n'), 'given by',sprintf('\n'), 'the ALL',sprintf('\n'), 'model{\times}10^{5}']);
set(leg2,'FontSize',fleg,'FontWeight','bold');title('(b)','FontSize',ft,'FontWeight','normal');
clabel(call1,hall1,'FontSize',labs,'FontWeight','normal');clabel(call2,hall2,'FontSize',labs,'FontWeight','normal');hall1.LineWidth = lw;hall2.LineWidth = lw;
xlabel('Control factor x_1','FontSize',labs,'FontWeight','normal');ylabel('Control factor x_2','FontSize',labs,'FontWeight','normal');
hall1.LabelSpacing=ls3;hall2.LabelSpacing=ls2;
xlim([0,1]);xticks([0:0.2:1]);ca2 = get(gca,'XTickLabel');set(gca,'XTickLabel',ca2,'fontsize',ftick,'fontweight','normal');
yticks([0:0.2:1]);cay2 = get(gca,'YTickLabel');set(gca,'YTickLabel',cay2,'fontsize',ftick,'fontweight','normal');
subplot(1,3,3);subplot('position',[fx+wdx*2+wx*2+0.01, fy,wx-0.01,wy]);
[cagpl1,hagpl1]=contour(xc1,xc2,qgrid*1e5,rangeGPL,'--r','ShowText','on');hold on;
h5 = plot(NaN, '-.r','LineWidth',lw);
[cagpl2,hagpl2]=contour(xc1,xc2,expectedlossAGPL*1e5,rangeGPL,'ShowText','on','EdgeColor',[0 0.4 0]);
h6 = plot(NaN, '-','Color',[0 0.4 0],'LineWidth',lw);
leg3=legend([h5 h6],['True',sprintf('\n'),'EQL{\times}10^{5}'],[' ',sprintf('\n'),'Posterior',sprintf('\n'),'mean of',sprintf('\n'), 'the EQL',sprintf('\n'), 'given by',sprintf('\n'), 'the AGPL',sprintf('\n'), 'model{\times}10^{5}']);
set(leg3,'FontSize',fleg,'FontWeight','bold');title('(c)','FontSize',ft,'FontWeight','normal');
clabel(cagpl1,hagpl1,'FontSize',labs,'FontWeight','normal');clabel(cagpl2,hagpl2,'FontSize',labs,'FontWeight','normal');hagpl1.LineWidth = lw;hagpl2.LineWidth = lw;
xlabel('Control factor x_1','FontSize',labs,'FontWeight','normal');ylabel('Control factor x_2','FontSize',labs,'FontWeight','normal');
hagpl1.LabelSpacing=ls1;hagpl2.LabelSpacing=ls2;
xlim([0,1]);xticks([0:0.2:1]);ca3 = get(gca,'XTickLabel');set(gca,'XTickLabel',ca3,'fontsize',ftick,'fontweight','normal');
yticks([0:0.2:1]);cay3= get(gca,'YTickLabel');set(gca,'YTickLabel',cay3,'fontsize',ftick,'fontweight','normal');
