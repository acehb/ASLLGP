clear all;
load('Designs for the piezoelectric actuator example.mat');
load('Piezo True EQL on grid points.mat')
load('400 test points for the piezoelectric actuator example.mat');
Lhall=zeros(m,ndesign);Lhd_all=zeros(m,ndesign);
Llall=zeros(n,ndesign);Lld_all=zeros(n,ndesign);
for k=1:ndesign
    for i=1:m
        [Lhall(i,k),Lhd_all(i,k)]=lossPiezo(yhd_all(i,k),ym_all(i,k));
    end
    for i=1:n
        [Llall(i,k),Lld_all(i,k)]=lossPiezo(yld_all(i,k),ym_all(i,k));
    end
end

MARE=zeros(ndesign,3);
MRCIL=zeros(ndesign,3);
cover=zeros(ndesign,3);
ARE95=zeros(ndesign,3);
RCIL95=zeros(ndesign,3);
temp1=repmat(ri_noise,1,length(ri_noise));
points=[reshape(temp1,[],1),reshape(temp1',[],1)];

[indextrue1,indextrue2]=find(qgrid==min(min(qgrid)));xcmin=[xc1(indextrue1,indextrue2),xc2(indextrue1,indextrue2)];%
xcminASLL=zeros(ndesign,2);
xcminALL=zeros(ndesign,2);xcminAGPL=zeros(ndesign,2);
minALLEQL=zeros(ndesign,1);minASLLEQL=zeros(ndesign,1);minAGPLEQL=zeros(ndesign,1);
xcminASLLind=zeros(ndesign,2);xcminALLind=zeros(ndesign,2);xcminAGPLind=zeros(ndesign,2);

para_ASLL=zeros(ndesign,12);ASLLbetalall=zeros(ndesign,1);ASLLtao2lall=zeros(ndesign,1);ASLLirReslall=zeros(n,ndesign);
ASLLbetaall=zeros(ndesign,1);ASLLrouall=zeros(ndesign,1);ASLLtao2all=zeros(ndesign,1);ASLLirResall=zeros(m,ndesign);
para_AGPL=zeros(ndesign, 8);AGPLbetalall=zeros(ndesign,1);AGPLtao2lall=zeros(ndesign,1);AGPLirReslall=zeros(n,ndesign);
AGPLbetaall=zeros(ndesign,1);AGPLrouall=zeros(ndesign,1);AGPLtao2all=zeros(ndesign,1);AGPLirResall=zeros(m,ndesign);

para_ALL=zeros(ndesign, 8);ALLbetalall=zeros(ndesign,1);ALLtao2lall=zeros(ndesign,1);ALLirReslall=zeros(n,ndesign);
ALLbetaall=zeros(ndesign,1);ALLrouall=zeros(ndesign,1);ALLtao2all=zeros(ndesign,1);ALLirResall=zeros(m,ndesign);

for k=1:ndesign
    %fit 3 models
    xl=reshape(xlall(k,:),d,[])';x=xl(1:m,:);%
    Lh=Lhall(:,k);Ll=Llall(:,k);Lhl=Llall(1:m,k);
    [ASLLthetal,ASLLfvall,ASLLbetal,ASLLtao2l,ASLLirResl,ASLLirxl,transparl]=gpfitASLL1level(xl,Ll);
    ASLLbetalall(k)=ASLLbetal;ASLLtao2lall(k)=ASLLtao2l;ASLLirReslall(:,k)=ASLLirResl;
    ASLLhl=log(transparl(2)*Lhl+transparl(1));
    [ASLLtheta,ASLLfval,ASLLrou,ASLLbeta,ASLLtao2,ASLLirRes,ASLLirx,transpar]=gpfitASLL2level(x,ASLLhl,Lh);
     para_ASLL(k,1:4)=ASLLthetal;  para_ASLL(k,5:8)=ASLLtheta; para_ASLL(k,9:12)=[transparl(2),transpar(2),transparl(1),transpar(1)];   
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
%get the summaries of the five EQL prediction performance measures for the ASLLGP, ALL, and AGPL models           
    lossLCL_ASLL=zeros(nt,1);
    lossUCL_ASLL=zeros(nt,1);
    ASLLEQL=zeros(nt,1);
    parfor ii=1:nt
        [ASLLEQL(ii),lossLCL_ASLL(ii),lossUCL_ASLL(ii)]=ASLLexpectedloss(x_e,w_noise,transpar,xl,ASLLthetal,xctest(ii,:),ASLLbetal,ASLLtao2l,ASLLirResl,ASLLirxl...
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
    parfor ii=1:nt
        [AGPLEQL(ii),lossLCL_AGPL(ii),lossUCL_AGPL(ii)]=AGPLexpectedloss(x_e,w_noise,xl,AGPLthetal,xctest(ii,:),AGPLbetal,AGPLtao2l,AGPLirResl,AGPLirxl,...
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
    parfor ii=1:nt
        [ALLEQL(ii),lossLCL_ALL(ii),lossUCL_ALL(ii)]=ALLexpectedloss(x_e,w_noise,xl,ALLthetal,xctest(ii,:),ALLbetal,ALLtao2l,ALLirResl,ALLirxl,...
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
%% get the summaries of the robust optimization performance of the ASLLGP, ALL, and AGPL models.
expectedlossASLL=zeros(nx,nx);
parfor i=1:nx
    for j=1:nx
        [expectedlossASLL(i,j),~,~]=ASLLexpectedloss(points,w_noise,transpar,xl,ASLLthetal,[xc1(i,j),xc2(i,j)],ASLLbetal,ASLLtao2l,ASLLirResl,ASLLirxl...
            ,x,ASLLtheta,ASLLrou,ASLLbeta,ASLLtao2,ASLLirRes,ASLLirx);
    end
end
[indexASLL1,indexASLL2]=find(expectedlossASLL==min(min(expectedlossASLL)));xcminASLL(k,:)=[xc1(indexASLL1,indexASLL2),xc2(indexASLL1,indexASLL2)];
xcminASLLind(k,:)=[indexASLL1,indexASLL2];minASLLEQL(k)=qgrid(indexASLL1,indexASLL2);

   expectedlossAGPL=zeros(nx,nx);
    parfor i=1:nx
        for j=1:nx
        [expectedlossAGPL(i,j),~,~]=AGPLexpectedloss(points,w_noise,xl,AGPLthetal,[xc1(i,j),xc2(i,j)],AGPLbetal,AGPLtao2l,AGPLirResl,AGPLirxl,...
            x,AGPLtheta,AGPLrou,AGPLbeta,AGPLtao2,AGPLirRes,AGPLirx);
        end
    end
    [indexAGPL1,indexAGPL2]=find(expectedlossAGPL==min(min(expectedlossAGPL)));xcminAGPL(k,:)=[xc1(indexAGPL1,indexAGPL2),xc2(indexAGPL1,indexAGPL2)];
xcminAGPLind(k,:)=[indexAGPL1,indexAGPL2];minAGPLEQL(k)=qgrid(indexAGPL1,indexAGPL2);

     expectedlossALL=zeros(nx,nx);
      parfor i=1:nx
        for j=1:nx
        [ expectedlossALL(i,j),~,~]=ALLexpectedloss(points,w_noise,xl,ALLthetal,[xc1(i,j),xc2(i,j)],ALLbetal,ALLtao2l,ALLirResl,ALLirxl,...
            x,ALLtheta,ALLrou,ALLbeta,ALLtao2,ALLirRes,ALLirx);
        end
      end
    [indexALL1,indexALL2]=find(expectedlossALL==min(min(expectedlossALL)));xcminALL(k,:)=[xc1(indexALL1,indexALL2),xc2(indexALL1,indexALL2)];
xcminALLind(k,:)=[indexALL1,indexALL2]; minALLEQL(k)=qgrid(indexALL1,indexALL2);
end

disp('Minimum true EQL');
minq=min(min(qgrid))
disp('True optimal control-factor setting');
xcmin
disp('Mean of true EQL');
expq=trapz(trapz(qgrid))*(0.005^2)
minEQL=[minASLLEQL,minALLEQL,minAGPLEQL];
disp('Table 6');
% five performance measures and the paired sample t-statistics
disp('Sample means for MARE given by ASLLGP, ALL, AGPL models');
mMARE=mean(MARE)*100
disp('Standard errors for MARE given by ASLLGP, ALL, AGPL models');
stderrMARE=std(MARE)./sqrt(ndesign)*100
disp('Paired t-statistics for MARE given by ASLLGP, ALL, AGPL models');
[~,~,~,s11]=ttest(MARE(:,2)-MARE(:,1));[~,~,~,s12]=ttest(MARE(:,3)-MARE(:,1));
tMARE=[0,s11.tstat,s12.tstat]
disp('Sample means for ARE_0.95 given by ASLLGP, ALL, AGPL models');
mARE95=mean(ARE95)*100
disp('Standard errors for ARE_0.95 given by ASLLGP, ALL, AGPL models');
stderrARE95=std(ARE95)./sqrt(ndesign)*100
disp('Paired t-statistics for ARE_0.95 given by ASLLGP, ALL, AGPL models');
[~,~,~,s41]=ttest(ARE95(:,2)-ARE95(:,1));[~,~,~,s42]=ttest(ARE95(:,3)-ARE95(:,1));
tARE95=[0,s41.tstat,s42.tstat]
disp('Sample means for MRCIL given by ASLLGP, ALL, AGPL models');
mMRCIL=mean(MRCIL)*100
disp('Standard errors for MRCIL given by ASLLGP, ALL, AGPL models');
stderrMRCIL=std(MRCIL)./sqrt(ndesign)*100
disp('Paired t-statistics for MRCIL given by ASLLGP, ALL, AGPL models');
[~,~,~,s21]=ttest(MRCIL(:,2)-MRCIL(:,1));[~,~,~,s22]=ttest(MRCIL(:,3)-MRCIL(:,1));
tMRCIL=[0,s21.tstat,s22.tstat]
disp('Sample means for RCIL_0.95 given by ASLLGP, ALL, AGPL models');
mRCIL95=mean(RCIL95)*100
disp('Standard errors for RCIL_0.95 given by ASLLGP, ALL, AGPL models');
stderrRCIL95=std(RCIL95)./sqrt(ndesign)*100
disp('Paired t-statistics for RCIL_0.95 given by ASLLGP, ALL, AGPL models');
[~,~,~,s51]=ttest(RCIL95(:,2)-RCIL95(:,1));[~,~,~,s52]=ttest(RCIL95(:,3)-RCIL95(:,1));
tRCIL95=[0,s51.tstat,s52.tstat]
disp('Sample means for EC-95 given by ASLLGP, ALL, AGPL models');
mcover=mean(cover)
disp('Standard errors for EC-95 given by ASLLGP, ALL, AGPL models');
stderrcover=std(cover)./sqrt(ndesign)
disp('Paired t-statistics for EC-95 given by ASLLGP, ALL, AGPL models');
[~,~,~,s31]=ttest(cover(:,2)-cover(:,1));[~,~,~,s32]=ttest(cover(:,3)-cover(:,1));
tcover=[0,s31.tstat,s32.tstat]
%% get the summaries of the five EQL prediction performance measures for the single-fidelity GP models and the two sample t-statistics.
Lallonefi=zeros(nonefi,ndesign);Ld_allonefi=zeros(nonefi,ndesign);
for k=1:ndesign
    for i=1:nonefi
        [Lallonefi(i,k),Ld_allonefi(i,k)]=lossPiezo(yd_allonefi(i,k),ym_allonefi(i,k));
    end
end

MAREonefi=zeros(ndesign,3);
MRCILonefi=zeros(ndesign,3);
coveronefi=zeros(ndesign,3);
ARE95onefi=zeros(ndesign,3);
RCIL95onefi=zeros(ndesign,3);
transpar_onefi=zeros(ndesign,2);

xcminSLLonefiMILE=zeros(ndesign,2);
xcminLLonefiMILE=zeros(ndesign,2);xcminGPLonefiMILE=zeros(ndesign,2);
minLLEQLonefi=zeros(ndesign,1);minSLLEQLonefi=zeros(ndesign,1);minGPLEQLonefi=zeros(ndesign,1);
xcminSLLonefiMILEind=zeros(ndesign,2);xcminLLonefiMILEind=zeros(ndesign,2);xcminGPLonefiMILEind=zeros(ndesign,2);

temp1=repmat(ri_noise,1,length(ri_noise));
points=[reshape(temp1,[],1),reshape(temp1',[],1)];

para_SLL=zeros(ndesign, 6);SLLbetaall=zeros(ndesign,1);SLLtao2all=zeros(ndesign,1);SLLirResall=zeros(nonefi,ndesign);
para_GPL=zeros(ndesign, 4);GPLbetalall=zeros(ndesign,1);GPLtao2lall=zeros(ndesign,1);GPLirReslall=zeros(nonefi,ndesign);
para_LL=zeros(ndesign, 4);LLbetalall=zeros(ndesign,1);LLtao2lall=zeros(ndesign,1);LLirReslall=zeros(nonefi,ndesign);

for k=1:ndesign
    xonefi=reshape(xallonefi(k,:),d,[])';
    Lonefi=Lallonefi(:,k);
    [SLLtheta_onefiMILE,SLLfvall_onefiMILE,SLLbeta_onefiMILE,SLLtao2_onefiMILE,SLLirRes_onefiMILE,SLLirx_onefiMILE,~,~,transparonefiMILE]=gpfitSLL(xonefi,Lonefi);

    transpar_onefi(k,:)=transparonefiMILE;para_SLL(k,1:4)=SLLtheta_onefiMILE; para_SLL(k,5:6)=[transparonefiMILE(2),transparonefiMILE(1)];
SLLbetaall(k)=SLLbeta_onefiMILE;SLLtao2all(k)=SLLtao2_onefiMILE;SLLirResall(:,k)=SLLirRes_onefiMILE;

    lossLCL_SLLonefiMILE=zeros(nt,1);
    lossUCL_SLLonefiMILE=zeros(nt,1);
    
    SLLEQLonefiMILE=zeros(nt,1);
    parfor ii=1:nt
        [SLLEQLonefiMILE(ii),lossLCL_SLLonefiMILE(ii),lossUCL_SLLonefiMILE(ii)]=ASLLexpectedloss(x_e,w_noise,transparonefiMILE,xonefi,SLLtheta_onefiMILE,xctest(ii,:),SLLbeta_onefiMILE,SLLtao2_onefiMILE,SLLirRes_onefiMILE,SLLirx_onefiMILE);
        
    end
    MARESLLonefiMILE=median(abs(qtest-SLLEQLonefiMILE)./qtest);
    MRCILSLLonefiMILE=median((lossUCL_SLLonefiMILE-lossLCL_SLLonefiMILE)./qtest);
    coverSLLonefiMILE=mean((qtest>=lossLCL_SLLonefiMILE).*(qtest<=lossUCL_SLLonefiMILE)*100);
    ARE95SLLonefiMILE=prctile(abs(qtest-SLLEQLonefiMILE)./qtest,95);
    RCIL95SLLonefiMILE=prctile((lossUCL_SLLonefiMILE-lossLCL_SLLonefiMILE)./qtest,95);
    
    [GPLtheta_onefiMILE,GPLfval_onefiMILE,GPLbeta_onefiMILE,GPLtao2_onefiMILE,GPLirRes_onefiMILE,GPLirx_onefiMILE,~,~]=gpfitGPL(xonefi,Lonefi);
    para_GPL(k,:)=GPLtheta_onefiMILE;GPLbetalall(k)=GPLbeta_onefiMILE;GPLtao2lall(k)=GPLtao2_onefiMILE;GPLirReslall(:,k)=GPLirRes_onefiMILE;
    lossLCL_GPLonefiMILE=zeros(nt,1);
    lossUCL_GPLonefiMILE=zeros(nt,1);
    GPLEQLonefiMILE=zeros(nt,1);
    parfor ii=1:nt
        [GPLEQLonefiMILE(ii),lossLCL_GPLonefiMILE(ii),lossUCL_GPLonefiMILE(ii)]=AGPLexpectedloss(x_e,w_noise,xonefi,GPLtheta_onefiMILE,xctest(ii,:),GPLbeta_onefiMILE,GPLtao2_onefiMILE,GPLirRes_onefiMILE,GPLirx_onefiMILE);
    end
    MAREGPLonefiMILE=median(abs(qtest-GPLEQLonefiMILE)./qtest);
    MRCILGPLonefiMILE=median((lossUCL_GPLonefiMILE-lossLCL_GPLonefiMILE)./qtest);
    coverGPLonefiMILE=mean((qtest>=lossLCL_GPLonefiMILE).*(qtest<=lossUCL_GPLonefiMILE)*100);
    ARE95GPLonefiMILE=prctile(abs(qtest-GPLEQLonefiMILE)./qtest,95);
    RCIL95GPLonefiMILE=prctile((lossUCL_GPLonefiMILE-lossLCL_GPLonefiMILE)./qtest,95);
    
    [LLtheta_onefiMILE,LLfval_onefiMILE,LLbeta_onefiMILE,LLtao2_onefiMILE,LLirRes_onefiMILE,LLirx_onefiMILE,~,~]=gpfitLL(xonefi,Lonefi);
    para_LL(k,:)=LLtheta_onefiMILE;LLbetalall(k)=LLbeta_onefiMILE;LLtao2lall(k)=LLtao2_onefiMILE;LLirReslall(:,k)=LLirRes_onefiMILE;
    lossLCL_LLonefiMILE=zeros(nt,1);
    lossUCL_LLonefiMILE=zeros(nt,1);
    LLEQLonefiMILE=zeros(nt,1);
    parfor ii=1:nt
        [LLEQLonefiMILE(ii),lossLCL_LLonefiMILE(ii),lossUCL_LLonefiMILE(ii)]=ALLexpectedloss(x_e,w_noise,xonefi,LLtheta_onefiMILE,xctest(ii,:),LLbeta_onefiMILE,LLtao2_onefiMILE,LLirRes_onefiMILE,LLirx_onefiMILE);
    end
    MARELLonefiMILE=median(abs(qtest-LLEQLonefiMILE)./qtest);
    MRCILLLonefiMILE=median((lossUCL_LLonefiMILE-lossLCL_LLonefiMILE)./qtest);
    coverLLonefiMILE=mean((qtest>=lossLCL_LLonefiMILE).*(qtest<=lossUCL_LLonefiMILE)*100);
    ARE95LLonefiMILE=prctile(abs(qtest-LLEQLonefiMILE)./qtest,95);
    RCIL95LLonefiMILE=prctile((lossUCL_LLonefiMILE-lossLCL_LLonefiMILE)./qtest,95);
    
    MAREonefi(k,:)=[MARESLLonefiMILE,MARELLonefiMILE,MAREGPLonefiMILE];
    MRCILonefi(k,:)=[MRCILSLLonefiMILE,MRCILLLonefiMILE,MRCILGPLonefiMILE];
    coveronefi(k,:)=[coverSLLonefiMILE,coverLLonefiMILE,coverGPLonefiMILE];
    ARE95onefi(k,:)=[ARE95SLLonefiMILE,ARE95LLonefiMILE,ARE95GPLonefiMILE];
    RCIL95onefi(k,:)=[RCIL95SLLonefiMILE,RCIL95LLonefiMILE,RCIL95GPLonefiMILE];
    %% get the summaries of the robust optimization performance of shifted log loss GP, lognormal loss, and GP for loss models.
expectedlossSLLonefiMILE=zeros(nx,nx);
parfor i=1:nx
    for j=1:nx
        [expectedlossSLLonefiMILE(i,j),~,~]=ASLLexpectedloss(points,w_noise,transparonefiMILE,xonefi,SLLtheta_onefiMILE,[xc1(i,j),xc2(i,j)],SLLbeta_onefiMILE,SLLtao2_onefiMILE,SLLirRes_onefiMILE,SLLirx_onefiMILE);
            
    end
end
[indexSLL1onefiMILE,indexSLL2onefiMILE]=find(expectedlossSLLonefiMILE==min(min(expectedlossSLLonefiMILE)));xcminSLLonefiMILE(k,:)=[xc1(indexSLL1onefiMILE,indexSLL2onefiMILE),xc2(indexSLL1onefiMILE,indexSLL2onefiMILE)];
xcminSLLonefiMILEind(k,:)=[indexSLL1onefiMILE,indexSLL2onefiMILE];minSLLEQLonefi(k)=qgrid(indexSLL1onefiMILE,indexSLL2onefiMILE);

   expectedlossGPLonefiMILE=zeros(nx,nx);
    parfor i=1:nx
        for j=1:nx
        [expectedlossGPLonefiMILE(i,j),~,~]=AGPLexpectedloss(points,w_noise,xonefi,GPLtheta_onefiMILE,[xc1(i,j),xc2(i,j)],GPLbeta_onefiMILE,GPLtao2_onefiMILE,GPLirRes_onefiMILE,GPLirx_onefiMILE);
            
        end
    end
    [indexGPL1onefiMILE,indexGPL2onefiMILE]=find(expectedlossGPLonefiMILE==min(min(expectedlossGPLonefiMILE)));xcminGPLonefiMILE(k,:)=[xc1(indexGPL1onefiMILE,indexGPL2onefiMILE),xc2(indexGPL1onefiMILE,indexGPL2onefiMILE)];
xcminGPLonefiMILEind(k,:)=[indexGPL1onefiMILE,indexGPL2onefiMILE];minGPLEQLonefi(k)=qgrid(indexGPL1onefiMILE,indexGPL2onefiMILE);

     expectedlossLLonefiMILE=zeros(nx,nx);
      parfor i=1:nx
        for j=1:nx
        [ expectedlossLLonefiMILE(i,j),~,~]=ALLexpectedloss(points,w_noise,xonefi,LLtheta_onefiMILE,[xc1(i,j),xc2(i,j)],LLbeta_onefiMILE,LLtao2_onefiMILE,LLirRes_onefiMILE,LLirx_onefiMILE);
           
        end
      end
    [indexLL1onefiMILE,indexLL2onefiMILE]=find(expectedlossLLonefiMILE==min(min(expectedlossLLonefiMILE)));xcminLLonefiMILE(k,:)=[xc1(indexLL1onefiMILE,indexLL2onefiMILE),xc2(indexLL1onefiMILE,indexLL2onefiMILE)];
xcminLLonefiMILEind(k,:)=[indexLL1onefiMILE,indexLL2onefiMILE];minLLEQLonefi(k)=qgrid(indexLL1onefiMILE,indexLL2onefiMILE);

end

disp('Table 7');
disp('Sample means for the true EQL at the estimated optimal control-factor setting given by ASLLGP, ALL, AGPL models');
MEQL=mean(minEQL)
disp('Standard errors for the true EQL at the estimated optimal control-factor setting given by ASLLGP, ALL, AGPL models');
STDEQL=std(minEQL)./sqrt(ndesign)
disp('Paired t-statistics for the true EQL at the estimated optimal control-factor setting given by ASLLGP, ALL, AGPL models');
[~,~,~,s1]=ttest(minEQL(:,2)-minEQL(:,1));[~,~,~,s2]=ttest(minEQL(:,3)-minEQL(:,1));
tEQL=[0,s1.tstat,s2.tstat]
disp('Sample means for the true EQL at the estimated optimal control-factor setting given by shifted log loss GP, lognormal loss, and GP for loss models');
minEQLonefi=[minSLLEQLonefi,minLLEQLonefi,minGPLEQLonefi];
MEQLonefi=mean(minEQLonefi)
disp('Standard errors for the true EQL at the estimated optimal control-factor setting given by shifted log loss GP, lognormal loss, and GP for loss models');
STDEQLonefi=std(minEQLonefi)./sqrt(ndesign)
[~,~,~,s0onefi]=ttest2(minEQLonefi(:,1),minEQL(:,1),'Vartype','unequal');[~,~,~,s1onefi]=ttest2(minEQLonefi(:,2),minEQL(:,1),'Vartype','unequal');[~,~,~,s2onefi]=ttest2(minEQLonefi(:,3),minEQL(:,1),'Vartype','unequal');
disp('Two sample t-statistics for the true EQL at the estimated optimal control-factor setting given by shifted log loss GP, lognormal loss, and GP for loss models');
tEQLonefi=[s0onefi.tstat,s1onefi.tstat,s2onefi.tstat]

disp('Table S3.1');
disp('Sample means for MARE given by shifted log loss GP, lognormal loss, and GP for loss models');
mMAREonefi=mean(MAREonefi)*100
disp('Standard errors for MARE given by shifted log loss GP, lognormal loss, and GP for loss models');
stderrMAREonefi=std(MAREonefi)./sqrt(ndesign)*100
disp('Two sample t-statistics for MARE given by shifted log loss GP, lognormal loss, and GP for loss models');
[~,~,~,stats11]=ttest2(MAREonefi(:,1),MARE(:,1),'Vartype','unequal');[~,~,~,stats12]=ttest2(MAREonefi(:,2),MARE(:,1),'Vartype','unequal');[~,~,~,stats13]=ttest2(MAREonefi(:,3),MARE(:,1),'Vartype','unequal');
tMAREdifffi=[stats11.tstat,stats12.tstat,stats13.tstat]
disp('Sample means for ARE_0.95 given by shifted log loss GP, lognormal loss, and GP for loss models');
mARE95onefi=mean(ARE95onefi)*100
disp('Standard errors for ARE_0.95 given by shifted log loss GP, lognormal loss, and GP for loss models');
stderrARE95onefi=std(ARE95onefi)./sqrt(ndesign)*100
disp('Two sample t-statistics for ARE_0.95 given by shifted log loss GP, lognormal loss, and GP for loss models');
[~,~,~,stats41]=ttest2(ARE95onefi(:,1),ARE95(:,1),'Vartype','unequal');[~,~,~,stats42]=ttest2(ARE95onefi(:,2),ARE95(:,1),'Vartype','unequal');[~,~,~,stats43]=ttest2(ARE95onefi(:,3),ARE95(:,1),'Vartype','unequal');
tARE95difffi=[stats41.tstat,stats42.tstat,stats43.tstat]
disp('Sample means for MRCIL given by shifted log loss GP, lognormal loss, and GP for loss models');
mMRCILonefi=mean(MRCILonefi)*100
disp('Standard errors for MRCIL given by shifted log loss GP, lognormal loss, and GP for loss models');
stderrMRCILonefi=std(MRCILonefi)./sqrt(ndesign)*100
disp('Two sample t-statistics for MRCIL given by shifted log loss GP, lognormal loss, and GP for loss models');
[~,~,~,stats21]=ttest2(MRCILonefi(:,1),MRCIL(:,1),'Vartype','unequal');[~,~,~,stats22]=ttest2(MRCILonefi(:,2),MRCIL(:,1),'Vartype','unequal');[~,~,~,stats23]=ttest2(MRCILonefi(:,3),MRCIL(:,1),'Vartype','unequal');
tMRCILdifffi=[stats21.tstat,stats22.tstat,stats23.tstat]
disp('Sample means for RCIL_0.95 given by shifted log loss GP, lognormal loss, and GP for loss models');
mRCIL95onefi=mean(RCIL95onefi)*100
disp('Standard errors for RCIL_0.95 given by shifted log loss GP, lognormal loss, and GP for loss models');
stderrRCIL95onefi=std(RCIL95onefi)./sqrt(ndesign)*100
disp('Two sample t-statistics for RCIL_0.95 given by shifted log loss GP, lognormal loss, and GP for loss models');
[~,~,~,stats51]=ttest2(RCIL95onefi(:,1),RCIL95(:,1),'Vartype','unequal');[~,~,~,stats52]=ttest2(RCIL95onefi(:,2),RCIL95(:,1),'Vartype','unequal');[~,~,~,stats53]=ttest2(RCIL95onefi(:,3),RCIL95(:,1),'Vartype','unequal');
tRCIL95difffi=[stats51.tstat,stats52.tstat,stats53.tstat]
disp('Sample means for EC-95 given by shifted log loss GP, lognormal loss, and GP for loss models');
mcover_onefi=mean(coveronefi)
disp('Standard errors for EC-95 given by shifted log loss GP, lognormal loss, and GP for loss models');
stderrcover_onefi=std(coveronefi)./sqrt(ndesign)
disp('Two sample t-statistics for EC-95 given by shifted log loss GP, lognormal loss, and GP for loss models');
[~,~,~,stats31]=ttest2(coveronefi(:,1),cover(:,1),'Vartype','unequal');[~,~,~,stats32]=ttest2(coveronefi(:,2),cover(:,1),'Vartype','unequal');[~,~,~,stats33]=ttest2(coveronefi(:,3),cover(:,1),'Vartype','unequal');
tcoverdifffi=[stats31.tstat,stats32.tstat,stats33.tstat]


