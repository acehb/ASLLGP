clear all;
load('Designs for the bar example.mat');
load('200 test points for the bar example.mat');
%BarTestpoints;
Lhall=zeros(m,ndesign);Lhs_all=zeros(m,ndesign);
Llall=zeros(n,ndesign);Lls_all=zeros(n,ndesign);
for k=1:ndesign
    for i=1:m
        [Lhall(i,k),Lhs_all(i,k)]=lossBar(yhs_all(i,k),ym_all(i,k));
    end
    for i=1:n
        [Llall(i,k),Lls_all(i,k)]=lossBar(yls_all(i,k),ym_all(i,k));
    end
end
nx=1001;x_c=linspace(0,1,nx)';
% compute the true EQL at the grid {0,0.001,бн,1}.
qgrid=zeros(nx,1);
parfor i=1:nx
    qgrid(i)=BarTrueEQL(x_c(i),stress,mass,w,x_e);
end


MARE=zeros(ndesign,3);
MRCIL=zeros(ndesign,3);
cover=zeros(ndesign,3);
ARE95=zeros(ndesign,3);
RCIL95=zeros(ndesign,3);
xcminASLL=zeros(ndesign,1);minASLL=zeros(ndesign,1);minASLLEQL=zeros(ndesign,1);
xcminAGPL=zeros(ndesign,1);minAGPL=zeros(ndesign,1);minAGPLEQL=zeros(ndesign,1);
xcminALL=zeros(ndesign,1);minALL=zeros(ndesign,1);minALLEQL=zeros(ndesign,1);
xcminASLLind=zeros(ndesign,1);xcminALLind=zeros(ndesign,1);xcminAGPLind=zeros(ndesign,1);
para_ASLL=zeros(ndesign, 8);
ASLLbetalall=zeros(ndesign,1);ASLLtao2lall=zeros(ndesign,1);ASLLirReslall=zeros(n,ndesign);
ASLLbetaall=zeros(ndesign,1);ASLLrouall=zeros(ndesign,1);ASLLtao2all=zeros(ndesign,1);ASLLirResall=zeros(m,ndesign);
para_AGPL=zeros(ndesign, 4);AGPLbetalall=zeros(ndesign,1);AGPLtao2lall=zeros(ndesign,1);AGPLirReslall=zeros(n,ndesign);
AGPLbetaall=zeros(ndesign,1);AGPLrouall=zeros(ndesign,1);AGPLtao2all=zeros(ndesign,1);AGPLirResall=zeros(m,ndesign);

para_ALL=zeros(ndesign, 4);ALLbetalall=zeros(ndesign,1);ALLtao2lall=zeros(ndesign,1);ALLirReslall=zeros(n,ndesign);
ALLbetaall=zeros(ndesign,1);ALLrouall=zeros(ndesign,1);ALLtao2all=zeros(ndesign,1);ALLirResall=zeros(m,ndesign);

for k=1:ndesign
    %fit 3 models
    xl=reshape(xlall(k,:),d,[])';x=xl(1:m,:);
    Lh=Lhall(:,k);Ll=Llall(:,k);Lhl=Llall(1:m,k);
    [ASLLthetal,ASLLfvall,ASLLbetal,ASLLtao2l,ASLLirResl,ASLLirxl,transparl]=gpfitASLL1level(xl,Ll);
    ASLLbetalall(k)=ASLLbetal;ASLLtao2lall(k)=ASLLtao2l;ASLLirReslall(:,k)=ASLLirResl;
    ASLLhl=log(transparl(2)*Lhl+transparl(1));
    [ASLLtheta,ASLLfval,ASLLrou,ASLLbeta,ASLLtao2,ASLLirRes,ASLLirx,transpar]=gpfitASLL2level(x,ASLLhl,Lh);
    para_ASLL(k,1:2)=ASLLthetal;  para_ASLL(k,3:4)=ASLLtheta;para_ASLL(k,5:8)=[transparl(2),transpar(2),transparl(1),transpar(1)];
    ASLLbetaall(k)=ASLLbeta;ASLLrouall(k)=ASLLrou;ASLLtao2all(k)=ASLLtao2;ASLLirResall(:,k)=ASLLirRes;
    [AGPLthetal,AGPLfvall,AGPLbetal,AGPLtao2l,AGPLirResl,AGPLirxl]=gpfit1level(xl,Ll);
    [AGPLtheta,AGPLfval,AGPLrou,AGPLbeta,AGPLtao2,AGPLirRes,AGPLirx]=gpfit2level(x,Lhl,Lh);
    para_AGPL(k,:)=[AGPLthetal,AGPLtheta];AGPLbetalall(k)=AGPLbetal;AGPLtao2lall(k)=AGPLtao2l;AGPLirReslall(:,k)=AGPLirResl;
    AGPLbetaall(k)=AGPLbeta;AGPLrouall(k)=AGPLrou;AGPLtao2all(k)=AGPLtao2;AGPLirResall(:,k)=AGPLirRes;

    ALLLh=log(Lh);ALLLl=log(Ll); ALLhl=ALLLl(1:m);
    [ALLthetal,ALLfvall,ALLbetal,ALLtao2l,ALLirResl,ALLirxl]=gpfit1level(xl,ALLLl);
    [ALLtheta,ALLfval,ALLrou,ALLbeta,ALLtao2,ALLirRes,ALLirx]=gpfit2level(x,ALLhl,ALLLh);
    para_ALL(k,:)=[ALLthetal,ALLtheta];ALLbetalall(k)=ALLbetal;ALLtao2lall(k)=ALLtao2l;ALLirReslall(:,k)=ALLirResl;
    ALLbetaall(k)=ALLbeta;ALLrouall(k)=ALLrou;ALLtao2all(k)=ALLtao2;ALLirResall(:,k)=ALLirRes;
    % get the summaries of the five EQL prediction performance measures for the ASLLGP, ALL, and AGPL models that are presented in Table 2 and the paired sample t-statistics in Table 2.
    lossLCL_ASLL=zeros(nt,1);
    lossUCL_ASLL=zeros(nt,1);
    ASLLEQL=zeros(nt,1);
    for ii=1:nt
        [ASLLEQL(ii),lossLCL_ASLL(ii),lossUCL_ASLL(ii)]=ASLLexpectedloss(x_e,w,transpar,xl,ASLLthetal,xctest(ii,:),ASLLbetal,ASLLtao2l,ASLLirResl,ASLLirxl...
            ,x,ASLLtheta,ASLLrou,ASLLbeta,ASLLtao2,ASLLirRes,ASLLirx);
    end
    MAREASLL=median(abs(qtest-ASLLEQL)./qtest);
    MRCILASLL=median((lossUCL_ASLL-lossLCL_ASLL)./qtest);
    coverASLL=mean((qtest>=lossLCL_ASLL).*(qtest<=lossUCL_ASLL)*100);
    ARE95ASLL=prctile(abs(qtest-ASLLEQL)./qtest,95);
    RCIL95ASLL=prctile((lossUCL_ASLL-lossLCL_ASLL)./qtest,95);

    lossLCL_AGPL=zeros(nt,1);
    lossUCL_AGPL=zeros(nt,1);
    AGPLEQL=zeros(nt,1);
    for ii=1:nt
        [AGPLEQL(ii),lossLCL_AGPL(ii),lossUCL_AGPL(ii)]=AGPLexpectedloss(x_e,w,xl,AGPLthetal,xctest(ii,:),AGPLbetal,AGPLtao2l,AGPLirResl,AGPLirxl,...
            x,AGPLtheta,AGPLrou,AGPLbeta,AGPLtao2,AGPLirRes,AGPLirx);
    end
    MAREAGPL=median(abs(qtest-AGPLEQL)./qtest);
    MRCILAGPL=median((lossUCL_AGPL-lossLCL_AGPL)./qtest);
    coverAGPL=mean((qtest>=lossLCL_AGPL).*(qtest<=lossUCL_AGPL)*100);
    ARE95AGPL=prctile(abs(qtest-AGPLEQL)./qtest,95);
    RCIL95AGPL=prctile((lossUCL_AGPL-lossLCL_AGPL)./qtest,95);

    lossLCL_ALL=zeros(nt,1);
    lossUCL_ALL=zeros(nt,1);
    ALLEQL=zeros(nt,1);
    for ii=1:nt
        [ALLEQL(ii),lossLCL_ALL(ii),lossUCL_ALL(ii)]=ALLexpectedloss(x_e,w,xl,ALLthetal,xctest(ii,:),ALLbetal,ALLtao2l,ALLirResl,ALLirxl,...
            x,ALLtheta,ALLrou,ALLbeta,ALLtao2,ALLirRes,ALLirx);
    end
    MAREALL=median(abs(qtest-ALLEQL)./qtest);
    MRCILALL=median((lossUCL_ALL-lossLCL_ALL)./qtest);
    coverALL=mean((qtest>=lossLCL_ALL).*(qtest<=lossUCL_ALL)*100);
    ARE95ALL=prctile(abs(qtest-ALLEQL)./qtest,95);
    RCIL95ALL=prctile((lossUCL_ALL-lossLCL_ALL)./qtest,95);

    MARE(k,:)=[MAREASLL,MAREALL,MAREAGPL];
    MRCIL(k,:)=[MRCILASLL,MRCILALL,MRCILAGPL];
    cover(k,:)=[coverASLL,coverALL,coverAGPL];
    ARE95(k,:)=[ARE95ASLL,ARE95ALL,ARE95AGPL];
    RCIL95(k,:)=[RCIL95ASLL,RCIL95ALL,RCIL95AGPL];
    %% get the summaries of the true EQL at the estimated optimal control-factor setting given by each of the ASLLGP, ALL and AGPL models
    EQLASLLgrid=zeros(nx,1);
    parfor i=1:nx
        [EQLASLLgrid(i),~,~]=ASLLexpectedloss(x_e,w,transpar,xl,ASLLthetal,x_c(i),ASLLbetal,ASLLtao2l,ASLLirResl,ASLLirxl...
            ,x,ASLLtheta,ASLLrou,ASLLbeta,ASLLtao2,ASLLirRes,ASLLirx);

    end

    [minASLL(k),minASLLtemp]=min(EQLASLLgrid);xcminASLL(k)=x_c(minASLLtemp);
    xcminASLLind(k,:)=minASLLtemp;minASLLEQL(k)=qgrid(minASLLtemp);

    EQLAGPLgrid=zeros(nx,1);
    parfor i=1:nx
        [EQLAGPLgrid(i),~,~]=AGPLexpectedloss(x_e,w,xl,AGPLthetal,x_c(i),AGPLbetal,AGPLtao2l,AGPLirResl,AGPLirxl,...
            x,AGPLtheta,AGPLrou,AGPLbeta,AGPLtao2,AGPLirRes,AGPLirx);

    end
    [minAGPL(k),minAGPLtemp]=min(EQLAGPLgrid);xcminAGPL(k)=x_c(minAGPLtemp);
    xcminAGPLind(k,:)=minAGPLtemp;minAGPLEQL(k)=qgrid(minAGPLtemp);

    EQLALLgrid=zeros(nx,1);
    parfor i=1:nx
        [EQLALLgrid(i),~,~]=ALLexpectedloss(x_e,w,xl,ALLthetal,x_c(i),ALLbetal,ALLtao2l,ALLirResl,ALLirxl,...
            x,ALLtheta,ALLrou,ALLbeta,ALLtao2,ALLirRes,ALLirx);

    end
    [minALL(k),minALLtemp]=min(EQLALLgrid);xcminALL(k)=x_c(minALLtemp);
    xcminALLind(k,:)=minALLtemp; minALLEQL(k)=qgrid(minALLtemp);
end
parfor k=1:ndesign
    minALLEQL(k)=BarTrueEQL(xcminALL(k),stress,mass,w,x_e);
    minASLLEQL(k)=BarTrueEQL(xcminASLL(k),stress,mass,w,x_e);
    minAGPLEQL(k)=BarTrueEQL(xcminAGPL(k),stress,mass,w,x_e);
end
disp('Minimum true EQL');
[minq,mintemp]=min(qgrid);
disp(minq);
disp('True optimal control-factor setting');
xcmin=x_c(mintemp)
Nx_e0=1;[x_e0,w0]=lgwt(Nx_e0,0,1);
parfor i=1:nx
    qgrid0(i)=BarTrueEQL(x_c(i),stress,mass,w0,x_e0);
end
[minq0,mintemp0]=min(qgrid0);
disp('Optimal control-factor setting when x_e is fixed at 0.5');
xcmin0=x_c(mintemp0)
disp('True EQL at the optimal control-factor setting when x_e is fixed at 0.5');
minqX0=BarTrueEQL(xcmin0,stress,mass,w,x_e)
disp('Table 2');
mMARE=mean(MARE)
disp('Sample means for MARE given by ASLLGP, ALL, AGPL models');
stderrMARE=std(MARE)./sqrt(ndesign)
disp('Standard errors for MARE given by ASLLGP, ALL, AGPL models');
[~,~,~,s11]=ttest(MARE(:,2)-MARE(:,1));[~,~,~,s12]=ttest(MARE(:,3)-MARE(:,1));
disp('Paired t-statistics for MARE given by ASLLGP, ALL, AGPL models');
tMARE=[0,s11.tstat,s12.tstat]
disp('Sample means for ARE_0.95 given by ASLLGP, ALL, AGPL models');
mARE95=mean(ARE95)
disp('Standard errors for ARE_0.95 given by ASLLGP, ALL, AGPL models');
stderrARE95=std(ARE95)./sqrt(ndesign)
disp('Paired t-statistics for ARE_0.95 given by ASLLGP, ALL, AGPL models');
[~,~,~,s41]=ttest(ARE95(:,2)-ARE95(:,1));[~,~,~,s42]=ttest(ARE95(:,3)-ARE95(:,1));
tARE95=[0,s41.tstat,s42.tstat]
disp('Sample means for MRCIL given by ASLLGP, ALL, AGPL models');
mMRCIL=mean(MRCIL)
disp('Standard errors for MRCIL given by ASLLGP, ALL, AGPL models');
stderrMRCIL=std(MRCIL)./sqrt(ndesign)
disp('Paired t-statistics for MRCIL given by ASLLGP, ALL, AGPL models');
[~,~,~,s21]=ttest(MRCIL(:,2)-MRCIL(:,1));[~,~,~,s22]=ttest(MRCIL(:,3)-MRCIL(:,1));
tMRCIL=[0,s21.tstat,s22.tstat]
disp('Sample means for RCIL_0.95 given by ASLLGP, ALL, AGPL models');
mRCIL95=mean(RCIL95)
disp('Standard errors for RCIL_0.95 given by ASLLGP, ALL, AGPL models');
stderrRCIL95=std(RCIL95)./sqrt(ndesign)
disp('Paired t-statistics for RCIL_0.95 given by ASLLGP, ALL, AGPL models');
[~,~,~,s51]=ttest(RCIL95(:,2)-RCIL95(:,1));[~,~,~,s52]=ttest(RCIL95(:,3)-RCIL95(:,1));
tRCIL95=[0,s51.tstat,s52.tstat]
disp('Sample means for EC-95 given by ASLLGP, ALL, AGPL models');
mcover=mean(cover)
disp('Standard errors for EC-95 given by ASLLGP, ALL, AGPL models');
stderrcover=std(cover)./sqrt(ndesign)
[~,~,~,s31]=ttest(cover(:,2)-cover(:,1));[~,~,~,s32]=ttest(cover(:,3)-cover(:,1));
disp('Paired t-statistics for EC-95 given by ASLLGP, ALL, AGPL models');
tcover=[0,s31.tstat,s32.tstat]

%% get the summaries of the five EQL prediction performance measures for the single-fidelity GP models and the two sample t-statistics in Table 3.
Lallonefi=zeros(nonefi,ndesign);Ls_allonefi=zeros(nonefi,ndesign);
for k=1:ndesign
    for i=1:nonefi
        [Lallonefi(i,k),Ls_allonefi(i,k)]=lossBar(ys_allonefi(i,k),ym_allonefi(i,k));
    end
end

MAREonefi=zeros(ndesign,3);
MRCILonefi=zeros(ndesign,3);
coveronefi=zeros(ndesign,3);
ARE95onefi=zeros(ndesign,3);
RCIL95onefi=zeros(ndesign,3);
transpar_SLL=zeros(ndesign,2);
xcminSLLonefiMILE=zeros(ndesign,1);minSLLonefiMILE=zeros(ndesign,1);minSLLEQLonefi=zeros(ndesign,1);
xcminGPLonefiMILE=zeros(ndesign,1);minGPLonefiMILE=zeros(ndesign,1);minGPLEQLonefi=zeros(ndesign,1);
xcminLLonefiMILE=zeros(ndesign,1);minLLonefiMILE=zeros(ndesign,1);minLLEQLonefi=zeros(ndesign,1);
xcminSLLonefiMILEind=zeros(ndesign,1);xcminLLonefiMILEind=zeros(ndesign,1);xcminGPLonefiMILEind=zeros(ndesign,1);

para_SLL=zeros(ndesign, 4);SLLbetaall=zeros(ndesign,1);SLLtao2all=zeros(ndesign,1);SLLirResall=zeros(nonefi,ndesign);
para_GPL=zeros(ndesign, 2);GPLbetalall=zeros(ndesign,1);GPLtao2lall=zeros(ndesign,1);GPLirReslall=zeros(nonefi,ndesign);

para_LL=zeros(ndesign, 2);LLbetalall=zeros(ndesign,1);LLtao2lall=zeros(ndesign,1);LLirReslall=zeros(nonefi,ndesign);

for k=1:ndesign
    xonefi=reshape(xallonefi(k,:),d,[])';
    Lonefi=Lallonefi(:,k);
    [SLLtheta_MILE,SLLfvall_MILE,SLLbeta_MILE,SLLtao2_MILE,SLLirRes_MILE,SLLirx_MILE,~,~,transpar_MILE]=gpfitSLL(xonefi,Lonefi);
    transpar_SLL(k,:)=transpar_MILE;para_SLL(k,1:2)=SLLtheta_MILE; para_SLL(k,3:4)=[transpar_MILE(2),transpar_MILE(1)];
    SLLbetaall(k)=SLLbeta_MILE;SLLtao2all(k)=SLLtao2_MILE;SLLirResall(:,k)=SLLirRes_MILE;

    lossLCL_SLLMILE=zeros(nt,1);
    lossUCL_SLLMILE=zeros(nt,1);
    SLLEQLMILE=zeros(nt,1);
    for ii=1:nt
        [SLLEQLMILE(ii),lossLCL_SLLMILE(ii),lossUCL_SLLMILE(ii)]=ASLLexpectedloss(x_e,w,transpar_MILE,xonefi,SLLtheta_MILE,xctest(ii,:),SLLbeta_MILE,SLLtao2_MILE,SLLirRes_MILE,SLLirx_MILE);

    end
    MARESLLMILE=median(abs(qtest-SLLEQLMILE)./qtest);
    MRCILSLLMILE=median((lossUCL_SLLMILE-lossLCL_SLLMILE)./qtest);
    coverSLLMILE=mean((qtest>=lossLCL_SLLMILE).*(qtest<=lossUCL_SLLMILE)*100);
    ARE95SLLMILE=prctile(abs(qtest-SLLEQLMILE)./qtest,95);
    RCIL95SLLMILE=prctile((lossUCL_SLLMILE-lossLCL_SLLMILE)./qtest,95);

    [GPLtheta_MILE,GPLfval_MILE,GPLbeta_MILE,GPLtao2_MILE,GPLirRes_MILE,GPLirx_MILE,~,~]=gpfitGPL(xonefi,Lonefi);
    para_GPL(k,:)=GPLtheta_MILE;GPLbetalall(k)=GPLbeta_MILE;GPLtao2lall(k)=GPLtao2_MILE;GPLirReslall(:,k)=GPLirRes_MILE;
    lossLCL_GPLMILE=zeros(nt,1);
    lossUCL_GPLMILE=zeros(nt,1);
    GPLEQLMILE=zeros(nt,1);
    for ii=1:nt
        [GPLEQLMILE(ii),lossLCL_GPLMILE(ii),lossUCL_GPLMILE(ii)]=AGPLexpectedloss(x_e,w,xonefi,GPLtheta_MILE,xctest(ii,:),GPLbeta_MILE,GPLtao2_MILE,GPLirRes_MILE,GPLirx_MILE);
    end
    MAREGPLMILE=median(abs(qtest-GPLEQLMILE)./qtest);
    MRCILGPLMILE=median((lossUCL_GPLMILE-lossLCL_GPLMILE)./qtest);
    coverGPLMILE=mean((qtest>=lossLCL_GPLMILE).*(qtest<=lossUCL_GPLMILE)*100);
    ARE95GPLMILE=prctile(abs(qtest-GPLEQLMILE)./qtest,95);
    RCIL95GPLMILE=prctile((lossUCL_GPLMILE-lossLCL_GPLMILE)./qtest,95);

    [LLtheta_MILE,LLfval_MILE,LLbeta_MILE,LLtao2_MILE,LLirRes_MILE,LLirx_MILE,~,~]=gpfitLL(xonefi,Lonefi);
    para_LL(k,:)=LLtheta_MILE;LLbetalall(k)=LLbeta_MILE;LLtao2lall(k)=LLtao2_MILE;LLirReslall(:,k)=LLirRes_MILE;
    lossLCL_LLMILE=zeros(nt,1);
    lossUCL_LLMILE=zeros(nt,1);
    LLEQLMILE=zeros(nt,1);
    for ii=1:nt
        [LLEQLMILE(ii),lossLCL_LLMILE(ii),lossUCL_LLMILE(ii)]=ALLexpectedloss(x_e,w,xonefi,LLtheta_MILE,xctest(ii,:),LLbeta_MILE,LLtao2_MILE,LLirRes_MILE,LLirx_MILE);
    end
    MARELLMILE=median(abs(qtest-LLEQLMILE)./qtest);
    MRCILLLMILE=median((lossUCL_LLMILE-lossLCL_LLMILE)./qtest);
    coverLLMILE=mean((qtest>=lossLCL_LLMILE).*(qtest<=lossUCL_LLMILE)*100);
    ARE95LLMILE=prctile(abs(qtest-LLEQLMILE)./qtest,95);
    RCIL95LLMILE=prctile((lossUCL_LLMILE-lossLCL_LLMILE)./qtest,95);

    MAREonefi(k,:)=[MARESLLMILE,MARELLMILE,MAREGPLMILE];
    MRCILonefi(k,:)=[MRCILSLLMILE,MRCILLLMILE,MRCILGPLMILE];
    coveronefi(k,:)=[coverSLLMILE,coverLLMILE,coverGPLMILE];
    ARE95onefi(k,:)=[ARE95SLLMILE,ARE95LLMILE,ARE95GPLMILE];
    RCIL95onefi(k,:)=[RCIL95SLLMILE,RCIL95LLMILE,RCIL95GPLMILE];

    %% get the summaries of the robust optimization performance of shifted log loss GP, lognormal loss, and GP for loss models.
    EQLSLLgridMILE=zeros(nx,1);
    parfor i=1:nx
        [EQLSLLgridMILE(i),~,~]=ASLLexpectedloss(x_e,w,transpar_MILE,xonefi,SLLtheta_MILE,x_c(i),SLLbeta_MILE,SLLtao2_MILE,SLLirRes_MILE,SLLirx_MILE);

    end

    [minSLLonefiMILE(k),minSLLtempMILE]=min(EQLSLLgridMILE);xcminSLLonefiMILE(k)=x_c(minSLLtempMILE);
    xcminSLLonefiMILEind(k,:)=minSLLtempMILE;minSLLEQLonefi(k)=qgrid(minSLLtempMILE);

    EQLGPLgridMILE=zeros(nx,1);
    parfor i=1:nx
        [EQLGPLgridMILE(i),~,~]=AGPLexpectedloss(x_e,w,xonefi,GPLtheta_MILE,x_c(i),GPLbeta_MILE,GPLtao2_MILE,GPLirRes_MILE,GPLirx_MILE);
    end
    [minGPLonefiMILE(k),minGPLtempMILE]=min(EQLGPLgridMILE);xcminGPLonefiMILE(k)=x_c(minGPLtempMILE);
    xcminGPLonefiMILEind(k,:)=minGPLtempMILE;minGPLEQLonefi(k)=qgrid(minGPLtempMILE);
    EQLLLgridMILE=zeros(nx,1);
    parfor i=1:nx
        [EQLLLgridMILE(i),~,~]=ALLexpectedloss(x_e,w,xonefi,LLtheta_MILE,x_c(i),LLbeta_MILE,LLtao2_MILE,LLirRes_MILE,LLirx_MILE);
    end
    [minLLonefiMILE(k),minLLtempMILE]=min(EQLLLgridMILE);xcminLLonefiMILE(k)=x_c(minLLtempMILE);
    xcminLLonefiMILEind(k,:)=minLLtempMILE;minLLEQLonefi(k)=qgrid(minLLtempMILE);
end
parfor k=1:ndesign
    minLLEQLonefi(k)=BarTrueEQL(xcminLLonefiMILE(k),stress,mass,w,x_e);
    minSLLEQLonefi(k)=BarTrueEQL(xcminSLLonefiMILE(k),stress,mass,w,x_e);
    minGPLEQLonefi(k)=BarTrueEQL(xcminGPLonefiMILE(k),stress,mass,w,x_e);
end
% five performance measures and two sample t-statistic
disp('Table 3');
disp('Sample means for MARE given by shifted log loss GP, lognormal loss, and GP for loss models');
mMAREonefi=mean(MAREonefi)
disp('Standard errors for MARE given by shifted log loss GP, lognormal loss, and GP for loss models');
stderrMAREonefi=std(MAREonefi)./sqrt(ndesign)
[~,~,~,stats11]=ttest2(MAREonefi(:,1),MARE(:,1),'Vartype','unequal');[~,~,~,stats12]=ttest2(MAREonefi(:,2),MARE(:,1),'Vartype','unequal');[~,~,~,stats13]=ttest2(MAREonefi(:,3),MARE(:,1),'Vartype','unequal');
disp('Two sample t-statistics for MARE given by shifted log loss GP, lognormal loss, and GP for loss models');
tMAREdifffi=[stats11.tstat,stats12.tstat,stats13.tstat]
tMAREdf=[stats11.df,stats12.df,stats13.df];
disp('Sample means for ARE_0.95 given by shifted log loss GP, lognormal loss, and GP for loss models');
mARE95onefi=mean(ARE95onefi)
disp('Standard errors for ARE_0.95 given by shifted log loss GP, lognormal loss, and GP for loss models');
stderrARE95onefi=std(ARE95onefi)./sqrt(ndesign)
[~,~,~,stats41]=ttest2(ARE95onefi(:,1),ARE95(:,1),'Vartype','unequal');[~,~,~,stats42]=ttest2(ARE95onefi(:,2),ARE95(:,1),'Vartype','unequal');[~,~,~,stats43]=ttest2(ARE95onefi(:,3),ARE95(:,1),'Vartype','unequal');
tARE95difffi=[stats41.tstat,stats42.tstat,stats43.tstat]
disp('Two sample t-statistics for ARE_0.95 given by shifted log loss GP, lognormal loss, and GP for loss models');
tARE95df=[stats41.df,stats42.df,stats43.df];
disp('Sample means for MRCIL given by shifted log loss GP, lognormal loss, and GP for loss models');
mMRCILonefi=mean(MRCILonefi)
disp('Standard errors for MRCIL given by shifted log loss GP, lognormal loss, and GP for loss models');
stderrMRCILonefi=std(MRCILonefi)./sqrt(ndesign)
disp('Two sample t-statistics for MRCIL given by shifted log loss GP, lognormal loss, and GP for loss models');
[~,~,~,stats21]=ttest2(MRCILonefi(:,1),MRCIL(:,1),'Vartype','unequal');[~,~,~,stats22]=ttest2(MRCILonefi(:,2),MRCIL(:,1),'Vartype','unequal');[~,~,~,stats23]=ttest2(MRCILonefi(:,3),MRCIL(:,1),'Vartype','unequal');
tMRCILdifffi=[stats21.tstat,stats22.tstat,stats23.tstat]
tMRCILdf=[stats21.df,stats22.df,stats23.df];
disp('Sample means for RCIL_0.95 given by shifted log loss GP, lognormal loss, and GP for loss models');
mRCIL95onefi=mean(RCIL95onefi)
disp('Standard errors for RCIL_0.95 given by shifted log loss GP, lognormal loss, and GP for loss models');
stderrRCIL95onefi=std(RCIL95onefi)./sqrt(ndesign)
disp('Two sample t-statistics for RCIL_0.95 given by shifted log loss GP, lognormal loss, and GP for loss models');
[~,~,~,stats51]=ttest2(RCIL95onefi(:,1),RCIL95(:,1),'Vartype','unequal');[~,~,~,stats52]=ttest2(RCIL95onefi(:,2),RCIL95(:,1),'Vartype','unequal');[~,~,~,stats53]=ttest2(RCIL95onefi(:,3),RCIL95(:,1),'Vartype','unequal');
tRCIL95difffi=[stats51.tstat,stats52.tstat,stats53.tstat]
tRCIL95df=[stats51.df,stats52.df,stats53.df];
disp('Sample means for EC-95 given by shifted log loss GP, lognormal loss, and GP for loss models');
mcover_onefi=mean(coveronefi)
disp('Standard errors for EC-95 given by shifted log loss GP, lognormal loss, and GP for loss models');
stderrcover_onefi=std(coveronefi)./sqrt(ndesign)
[~,~,~,stats31]=ttest2(coveronefi(:,1),cover(:,1),'Vartype','unequal');[~,~,~,stats32]=ttest2(coveronefi(:,2),cover(:,1),'Vartype','unequal');[~,~,~,stats33]=ttest2(coveronefi(:,3),cover(:,1),'Vartype','unequal');
tcoverdifffi=[stats31.tstat,stats32.tstat,stats33.tstat]
disp('Two sample t-statistics for EC-95 given by shifted log loss GP, lognormal loss, and GP for loss models');
tcoverdf=[stats31.df,stats32.df,stats33.df];
disp('Table 4');
disp('Sample means for the true EQL at the estimated optimal control-factor setting given by ASLLGP, ALL, AGPL models');
minEQL=[minASLLEQL,minALLEQL,minAGPLEQL];
MEQL=mean(minEQL)
disp('Standard errors for the true EQL at the estimated optimal control-factor setting given by ASLLGP, ALL, AGPL models');
STDEQL=std(minEQL)./sqrt(ndesign)
[~,~,~,s1]=ttest(minEQL(:,2)-minEQL(:,1));[~,~,~,s2]=ttest(minEQL(:,3)-minEQL(:,1));
disp('Paired t-statistics for the true EQL at the estimated optimal control-factor setting given by ASLLGP, ALL, AGPL models');
tEQL=[0,s1.tstat,s2.tstat]
disp('Sample means for the true EQL at the estimated optimal control-factor setting given by shifted log loss GP, lognormal loss, and GP for loss models');
minEQLonefi=[minSLLEQLonefi,minLLEQLonefi,minGPLEQLonefi];
MEQLonefi=mean(minEQLonefi)
disp('Standard errors for the true EQL at the estimated optimal control-factor setting given by shifted log loss GP, lognormal loss, and GP for loss models');
STDEQLonefi=std(minEQLonefi)./sqrt(ndesign)
disp('Two sample t-statistics for the true EQL at the estimated optimal control-factor setting given by shifted log loss GP, lognormal loss, and GP for loss models');
[~,~,~,s0onefi]=ttest2(minEQLonefi(:,1),minEQL(:,1),'Vartype','unequal');[~,~,~,s1onefi]=ttest2(minEQLonefi(:,2),minEQL(:,1),'Vartype','unequal');[~,~,~,s2onefi]=ttest2(minEQLonefi(:,3),minEQL(:,1),'Vartype','unequal');
tEQLonefi=[s0onefi.tstat,s1onefi.tstat,s2onefi.tstat]

mindf=min([min(tMAREdf),min(tMRCILdf),min(tARE95df),min(tRCIL95df),min(tcoverdf)]);
maxdf=max([max(tMAREdf),max(tMRCILdf),max(tARE95df),max(tRCIL95df),max(tcoverdf)]);
tEQLdf=[s0onefi.df,s1onefi.df,s2onefi.df];

