% clear;clc
clear;
load('100 designs for the shaft example.mat');% 100 designs given by designpoints.m
load('200 test points for the shaft example.mat');% 200 test points given by testpoints.m

Lhall=zeros(m,ndesign);Lhs_all=zeros(m,ndesign);
Llall=zeros(n,ndesign);Lls_all=zeros(n,ndesign);
for k=1:ndesign
    for i=1:m
        [Lhall(i,k),Lhs_all(i,k)]=lossShaft(yhs_all(i,k),ym_all(i,k));
    end
    for i=1:n
        [Llall(i,k),Lls_all(i,k)]=lossShaft(yls_all(i,k),ym_all(i,k));
    end
end

Loss_ALL=zeros(ndesign, nt);loss_ASLL=zeros(ndesign, nt);Loss_AGPL=zeros(ndesign, nt);LossLCL_AGPL=zeros(ndesign, nt);LossUCL_AGPL=zeros(ndesign, nt);
LossLCL_ALL=zeros(ndesign, nt);LossUCL_ALL=zeros(ndesign, nt);LossLCL_ASLL=zeros(ndesign, nt);LossUCL_ASLL=zeros(ndesign, nt);
trans_ASLL=zeros(ndesign, 4);

temptest=zeros(Nx_e,nt);
L_stest=zeros(Nx_e,nt);
L_mtest=zeros(Nx_e,nt);
for i=1:nt
    for j=1:Nx_e
        [temptest(j,i),L_stest(j,i),L_mtest(j,i)]=lossShaft(ys_test(j,i),ym_test(j,i));
    end
end
Ltest1=(w*L_stest)';Ltest2=L_mtest(1,:)';
qtest=Ltest2+Ltest1;

MARE=zeros(ndesign,3);
MRCIL=zeros(ndesign,3);
cover=zeros(ndesign,3);
ARE95=zeros(ndesign,3);
RCIL95=zeros(ndesign,3);
for k=1:ndesign
    %fit 3 models
    xl=reshape(xlall(k,:),d,[])';x=xl(1:m,:);
    Lh=Lhall(:,k);Ll=Llall(:,k);Lhl=Llall(1:m,k);
    [ASLLthetal,ASLLfvall,ASLLbetal,ASLLtao2l,ASLLirResl,ASLLirxl,transparl]=gpfitASLL1level(xl,Ll);
    ASLLhl=log(transparl(2)*Lhl+transparl(1));
    [ASLLtheta,ASLLfval,ASLLrou,ASLLbeta,ASLLtao2,ASLLirRes,ASLLirx,transpar]=gpfitASLL2level(x,ASLLhl,Lh);
     trans_ASLL(k,1:2)=transparl;  trans_ASLL(k,3:4)=transpar;  
    [AGPLthetal,AGPLfvall,AGPLbetal,AGPLtao2l,AGPLirResl,AGPLirxl]=gpfit1level(xl,Ll);
    [AGPLtheta,AGPLfval,AGPLrou,AGPLbeta,AGPLtao2,AGPLirRes,AGPLirx]=gpfit2level(x,Lhl,Lh);
        
    ALLLh=log(Lh);ALLLl=log(Ll); ALLhl=ALLLl(1:m);
    [ALLthetal,ALLfvall,ALLbetal,ALLtao2l,ALLirResl,ALLirxl]=gpfit1level(xl,ALLLl);
    [ALLtheta,ALLfval,ALLrou,ALLbeta,ALLtao2,ALLirRes,ALLirx]=gpfit2level(x,ALLhl,ALLLh);
        
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
end
% Table 2: five performance measures in predicting EQL obtained with the ASLLGP, ALL, and AGPL models for the shaft example.
mMARE=mean(MARE)
stderrMARE=std(MARE)./sqrt(ndesign)
[~,~,~,s11]=ttest(MARE(:,2)-MARE(:,1));[~,~,~,s12]=ttest(MARE(:,3)-MARE(:,1));
tMARE=[0,s11.tstat,s12.tstat]

mMRCIL=mean(MRCIL)
stderrMRCIL=std(MRCIL)./sqrt(ndesign)
[~,~,~,s21]=ttest(MRCIL(:,2)-MRCIL(:,1));[~,~,~,s22]=ttest(MRCIL(:,3)-MRCIL(:,1));
tMRCIL=[0,s21.tstat,s22.tstat]

mARE95=mean(ARE95)
stderrARE95=std(ARE95)./sqrt(ndesign)
[~,~,~,s41]=ttest(ARE95(:,2)-ARE95(:,1));[~,~,~,s42]=ttest(ARE95(:,3)-ARE95(:,1));
tARE95=[0,s41.tstat,s42.tstat]

mRCIL95=mean(RCIL95)
stderrRCIL95=std(RCIL95)./sqrt(ndesign)
[~,~,~,s51]=ttest(RCIL95(:,2)-RCIL95(:,1));[~,~,~,s52]=ttest(RCIL95(:,3)-RCIL95(:,1));
tRCIL95=[0,s51.tstat,s52.tstat]

mcover=mean(cover)
stderrcover=std(cover)./sqrt(ndesign)
[~,~,~,s31]=ttest(cover(:,2)-cover(:,1));[~,~,~,s32]=ttest(cover(:,3)-cover(:,1));
tcover=[0,s31.tstat,s32.tstat]
toc
%% Table 3:Sample mean and its standard error for five EQL prediction performance measures, i.e., MARE, 95%ARE, MRCIL, 95%RCIL, and EC_95, given by the shifted log loss GP, lognormal loss, and GP for loss models. Each two sample t-statistic is for testing whether the difference between the mean of a performance measure given by the first model and that given by the second model in one of the following pairs of models is zero: (1) the shifted log loss GP and ASLLGP models (first column), (2) the lognormal loss and ALL models (second column), and (3) the GP for loss and AGPL models (third column).
tic
Lallonefi=zeros(nonefi,ndesign);Ls_allonefi=zeros(nonefi,ndesign);
for k=1:ndesign
    for i=1:nonefi
        [Lallonefi(i,k),Ls_allonefi(i,k)]=lossShaft(ys_allonefi(i,k),ym_allonefi(i,k));
    end
end

MAREonefi=zeros(ndesign,3);
MRCILonefi=zeros(ndesign,3);
coveronefi=zeros(ndesign,3);
ARE95onefi=zeros(ndesign,3);
RCIL95onefi=zeros(ndesign,3);
transpar_onefi=zeros(ndesign,2);

for k=1:ndesign
    xonefi=reshape(xallonefi(k,:),d,[])';
    Lonefi=Lallonefi(:,k);
    [ASLLtheta_onefi,ASLLfvall_onefi,ASLLbeta_onefi,ASLLtao2_onefi,ASLLirRes_onefi,ASLLirx_onefi,transpar]=gpfitASLL1level(xonefi,Lonefi);
    
    transpar_onefi(k,:)=transpar;
    lossLCL_ASLLonefi=zeros(nt,1);
    lossUCL_ASLLonefi=zeros(nt,1);
    
    ASLLEQLonefi=zeros(nt,1);
    for ii=1:nt
        [ASLLEQLonefi(ii),lossLCL_ASLLonefi(ii),lossUCL_ASLLonefi(ii)]=ASLLexpectedloss(x_e,w,transpar,xonefi,ASLLtheta_onefi,xctest(ii,:),ASLLbeta_onefi,ASLLtao2_onefi,ASLLirRes_onefi,ASLLirx_onefi);
        
    end
    MAREASLLonefi=median(abs(qtest-ASLLEQLonefi)./qtest);
    MRCILASLLonefi=median((lossUCL_ASLLonefi-lossLCL_ASLLonefi)./qtest);
    coverASLLonefi=mean((qtest>=lossLCL_ASLLonefi).*(qtest<=lossUCL_ASLLonefi)*100);
    ARE95ASLLonefi=prctile(abs(qtest-ASLLEQLonefi)./qtest,95);
    RCIL95ASLLonefi=prctile((lossUCL_ASLLonefi-lossLCL_ASLLonefi)./qtest,95);
    
    [AGPLtheta_onefi,AGPLfval_onefi,AGPLbeta_onefi,AGPLtao2_onefi,AGPLirRes_onefi,AGPLirx_onefi]=gpfit1level(xonefi,Lonefi);
    lossLCL_AGPLonefi=zeros(nt,1);
    lossUCL_AGPLonefi=zeros(nt,1);
    AGPLEQLonefi=zeros(nt,1);
    for ii=1:nt
        [AGPLEQLonefi(ii),lossLCL_AGPLonefi(ii),lossUCL_AGPLonefi(ii)]=AGPLexpectedloss(x_e,w,xonefi,AGPLtheta_onefi,xctest(ii,:),AGPLbeta_onefi,AGPLtao2_onefi,AGPLirRes_onefi,AGPLirx_onefi);
    end
    MAREAGPLonefi=median(abs(qtest-AGPLEQLonefi)./qtest);
    MRCILAGPLonefi=median((lossUCL_AGPLonefi-lossLCL_AGPLonefi)./qtest);
    coverAGPLonefi=mean((qtest>=lossLCL_AGPLonefi).*(qtest<=lossUCL_AGPLonefi)*100);
    ARE95AGPLonefi=prctile(abs(qtest-AGPLEQLonefi)./qtest,95);
    RCIL95AGPLonefi=prctile((lossUCL_AGPLonefi-lossLCL_AGPLonefi)./qtest,95);
    
    ALLLonefi=log(Lonefi);
    [ALLtheta_onefi,ALLfva_onefi,ALLbeta_onefi,ALLtao2_onefi,ALLirRes_onefi,ALLirx_onefi]=gpfit1level(xonefi,ALLLonefi);
    lossLCL_ALLonefi=zeros(nt,1);
    lossUCL_ALLonefi=zeros(nt,1);
    ALLEQLonefi=zeros(nt,1);
    for ii=1:nt
        [ALLEQLonefi(ii),lossLCL_ALLonefi(ii),lossUCL_ALLonefi(ii)]=ALLexpectedloss(x_e,w,xonefi,ALLtheta_onefi,xctest(ii,:),ALLbeta_onefi,ALLtao2_onefi,ALLirRes_onefi,ALLirx_onefi);
    end
    MAREALLonefi=median(abs(qtest-ALLEQLonefi)./qtest);
    MRCILALLonefi=median((lossUCL_ALLonefi-lossLCL_ALLonefi)./qtest);
    coverALLonefi=mean((qtest>=lossLCL_ALLonefi).*(qtest<=lossUCL_ALLonefi)*100);
    ARE95ALLonefi=prctile(abs(qtest-ALLEQLonefi)./qtest,95);
    RCIL95ALLonefi=prctile((lossUCL_ALLonefi-lossLCL_ALLonefi)./qtest,95);
    
    MAREonefi(k,:)=[MAREASLLonefi,MAREALLonefi,MAREAGPLonefi];
    MRCILonefi(k,:)=[MRCILASLLonefi,MRCILALLonefi,MRCILAGPLonefi];
    coveronefi(k,:)=[coverASLLonefi,coverALLonefi,coverAGPLonefi];
    ARE95onefi(k,:)=[ARE95ASLLonefi,ARE95ALLonefi,ARE95AGPLonefi];
    RCIL95onefi(k,:)=[RCIL95ASLLonefi,RCIL95ALLonefi,RCIL95AGPLonefi];
    
end
mMAREonefi=mean(MAREonefi)
stderrMAREonefi=std(MAREonefi)./sqrt(ndesign)
[~,~,~,stats11]=ttest2(MAREonefi(:,1),MARE(:,1),'Vartype','unequal');[~,~,~,stats12]=ttest2(MAREonefi(:,2),MARE(:,2),'Vartype','unequal');[~,~,~,stats13]=ttest2(MAREonefi(:,3),MARE(:,3),'Vartype','unequal');
tMAREdifffi=[stats11.tstat,stats12.tstat,stats13.tstat]
tMAREdf=[stats11.df,stats12.df,stats13.df];

mMRCILonefi=mean(MRCILonefi)
stderrMRCILonefi=std(MRCILonefi)./sqrt(ndesign)
[~,~,~,stats21]=ttest2(MRCILonefi(:,1),MRCIL(:,1),'Vartype','unequal');[~,~,~,stats22]=ttest2(MRCILonefi(:,2),MRCIL(:,2),'Vartype','unequal');[~,~,~,stats23]=ttest2(MRCILonefi(:,3),MRCIL(:,3),'Vartype','unequal');
tMRCILdifffi=[stats21.tstat,stats22.tstat,stats23.tstat]
tMRCILdf=[stats21.df,stats22.df,stats23.df];

mARE95onefi=mean(ARE95onefi)
stderrARE95onefi=std(ARE95onefi)./sqrt(ndesign)
[~,~,~,stats41]=ttest2(ARE95onefi(:,1),ARE95(:,1),'Vartype','unequal');[~,~,~,stats42]=ttest2(ARE95onefi(:,2),ARE95(:,2),'Vartype','unequal');[~,~,~,stats43]=ttest2(ARE95onefi(:,3),ARE95(:,3),'Vartype','unequal');
tARE95difffi=[stats41.tstat,stats42.tstat,stats43.tstat]
tARE95df=[stats41.df,stats42.df,stats43.df];

mRCIL95onefi=mean(RCIL95onefi)
stderrRCIL95onefi=std(RCIL95onefi)./sqrt(ndesign)
[~,~,~,stats51]=ttest2(RCIL95onefi(:,1),RCIL95(:,1),'Vartype','unequal');[~,~,~,stats52]=ttest2(RCIL95onefi(:,2),RCIL95(:,2),'Vartype','unequal');[~,~,~,stats53]=ttest2(RCIL95onefi(:,3),RCIL95(:,3),'Vartype','unequal');
tRCIL95difffi=[stats51.tstat,stats52.tstat,stats53.tstat]
tRCIL95df=[stats51.df,stats52.df,stats53.df];

mcover_onefi=mean(coveronefi)
stderrcover_onefi=std(coveronefi)./sqrt(ndesign)
[~,~,~,stats31]=ttest2(coveronefi(:,1),cover(:,1),'Vartype','unequal');[~,~,~,stats32]=ttest2(coveronefi(:,2),cover(:,2),'Vartype','unequal');[~,~,~,stats33]=ttest2(coveronefi(:,3),cover(:,3),'Vartype','unequal');
tcoverdifffi=[stats31.tstat,stats32.tstat,stats33.tstat]
tcoverdf=[stats31.df,stats32.df,stats33.df];

 mindf=min([min(tMAREdf),min(tMRCILdf),min(tARE95df),min(tRCIL95df),min(tcoverdf)])
 maxdf=max([max(tMAREdf),max(tMRCILdf),max(tARE95df),max(tRCIL95df),max(tcoverdf)])
