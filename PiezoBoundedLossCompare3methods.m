clear all;
load('Designs for the piezoelectric actuator example.mat');
load('Piezo True EQL on grid points.mat')
load('400 test points for the piezoelectric actuator example.mat');
ub=1;
LhallBounded=zeros(m,ndesign);LhdBounded_all=zeros(m,ndesign);
LlallBounded=zeros(n,ndesign);LldBounded_all=zeros(n,ndesign);
for k=1:ndesign
    for i=1:m
        [LhallBounded(i,k),LhdBounded_all(i,k)]=lossPiezo(yhd_all(i,k),ym_all(i,k),b);
    end
    for i=1:n
        [LlallBounded(i,k),LldBounded_all(i,k)]=lossPiezo(yld_all(i,k),ym_all(i,k),b);
    end
end

temp1=repmat(ri_noise,1,length(ri_noise));
points=[reshape(temp1,[],1),reshape(temp1',[],1)];

[indextrue1Bounded,indextrue2Bounded]=find(qgridBounded==min(min(qgridBounded)));xcminBounded=[xc1(indextrue1Bounded,indextrue2Bounded),xc2(indextrue1Bounded,indextrue2Bounded)];%
xcminASLLBounded=zeros(ndesign,2);
xcminALLBounded=zeros(ndesign,2);xcminAGPLBounded=zeros(ndesign,2);
minALLEQLBounded=zeros(ndesign,1);minASLLEQLBounded=zeros(ndesign,1);minAGPLEQLBounded=zeros(ndesign,1);
xcminASLLBoundedind=zeros(ndesign,2);xcminALLBoundedind=zeros(ndesign,2);xcminAGPLBoundedind=zeros(ndesign,2);
MAREBounded=zeros(ndesign,3);
MRCILBounded=zeros(ndesign,3);
coverBounded=zeros(ndesign,3);
ARE95Bounded=zeros(ndesign,3);
RCIL95Bounded=zeros(ndesign,3);
para_ASLLBounded=zeros(ndesign,12);
ASLLbetalallBounded=zeros(ndesign,1);ASLLtao2lallBounded=zeros(ndesign,1);ASLLirReslallBounded=zeros(n,ndesign);
ASLLbetaallBounded=zeros(ndesign,1);ASLLrouallBounded=zeros(ndesign,1);ASLLtao2allBounded=zeros(ndesign,1);ASLLirResallBounded=zeros(m,ndesign);
para_AGPLBounded=zeros(ndesign, 8);AGPLbetalallBounded=zeros(ndesign,1);AGPLtao2lallBounded=zeros(ndesign,1);AGPLirReslallBounded=zeros(n,ndesign);
AGPLbetaallBounded=zeros(ndesign,1);AGPLrouallBounded=zeros(ndesign,1);AGPLtao2allBounded=zeros(ndesign,1);AGPLirResallBounded=zeros(m,ndesign);

para_ALLBounded=zeros(ndesign, 8);ALLbetalallBounded=zeros(ndesign,1);ALLtao2lallBounded=zeros(ndesign,1);ALLirReslallBounded=zeros(n,ndesign);
ALLbetaallBounded=zeros(ndesign,1);ALLrouallBounded=zeros(ndesign,1);ALLtao2allBounded=zeros(ndesign,1);ALLirResallBounded=zeros(m,ndesign);


for k=1:ndesign
    %fit 3 models
    xl=reshape(xlall(k,:),d,[])';x=xl(1:m,:);%
    LhBounded=LhallBounded(:,k);LlBounded=LlallBounded(:,k);LhlBounded=LlallBounded(1:m,k);
    [ASLLthetalBounded,ASLLfvallBounded,ASLLbetalBounded,ASLLtao2lBounded,ASLLirReslBounded,ASLLirxlBounded,transparlBounded]=gpfitASLL1level(xl,LlBounded,ub);
    ASLLhlBounded=log(transparlBounded(2)*LhlBounded+transparlBounded(1));
    [ASLLthetaBounded,ASLLfvalBounded,ASLLrouBounded,ASLLbetaBounded,ASLLtao2Bounded,ASLLirResBounded,ASLLirxBounded,transparBounded]=gpfitASLL2level(x,ASLLhlBounded,LhBounded,ub);
    para_ASLLBounded(k,1:4)=ASLLthetalBounded;  para_ASLLBounded(k,5:8)=ASLLthetaBounded; para_ASLLBounded(k,9:12)=[transparlBounded(2),transparBounded(2),transparlBounded(1),transparBounded(1)];   
ASLLbetalallBounded(k)=ASLLbetalBounded;ASLLtao2lallBounded(k)=ASLLtao2lBounded;ASLLirReslallBounded(:,k)=ASLLirReslBounded;
ASLLbetaallBounded(k)=ASLLbetaBounded;ASLLrouallBounded(k)=ASLLrouBounded;ASLLtao2allBounded(k)=ASLLtao2Bounded;ASLLirResallBounded(:,k)=ASLLirResBounded;

    [AGPLthetalBounded,AGPLfvallBounded,AGPLbetalBounded,AGPLtao2lBounded,AGPLirReslBounded,AGPLirxlBounded]=gpfit1level(xl,LlBounded);
    [AGPLthetaBounded,AGPLfvalBounded,AGPLrouBounded,AGPLbetaBounded,AGPLtao2Bounded,AGPLirResBounded,AGPLirxBounded]=gpfit2level(x,LhlBounded,LhBounded);
   para_AGPLBounded(k,:)=[AGPLthetalBounded,AGPLthetaBounded];AGPLbetalallBounded(k)=AGPLbetalBounded;AGPLtao2lallBounded(k)=AGPLtao2lBounded;AGPLirReslallBounded(:,k)=AGPLirReslBounded;
AGPLbetaallBounded(k)=AGPLbetaBounded;AGPLrouallBounded(k)=AGPLrouBounded;AGPLtao2allBounded(k)=AGPLtao2Bounded;AGPLirResallBounded(:,k)=AGPLirResBounded;

    ALLLhBounded=log(LhBounded);ALLLlBounded=log(LlBounded); ALLhlBounded=ALLLlBounded(1:m);
    [ALLthetalBounded,ALLfvallBounded,ALLbetalBounded,ALLtao2lBounded,ALLirReslBounded,ALLirxlBounded]=gpfit1level(xl,ALLLlBounded);
    [ALLthetaBounded,ALLfvalBounded,ALLrouBounded,ALLbetaBounded,ALLtao2Bounded,ALLirResBounded,ALLirxBounded]=gpfit2level(x,ALLhlBounded,ALLLhBounded);
    para_ALLBounded(k,:)=[ALLthetalBounded,ALLthetaBounded];ALLbetalallBounded(k)=ALLbetalBounded;ALLtao2lallBounded(k)=ALLtao2lBounded;ALLirReslallBounded(:,k)=ALLirReslBounded;
    ALLbetaallBounded(k)=ALLbetaBounded;ALLrouallBounded(k)=ALLrouBounded;ALLtao2allBounded(k)=ALLtao2Bounded;ALLirResallBounded(:,k)=ALLirResBounded;

    %get the summaries of the five EQL prediction performance measures for the ASLLGP, ALL, and AGPL models
    lossLCLBounded_ASLL=zeros(nt,1);
    lossUCLBounded_ASLL=zeros(nt,1);
    ASLLEQLBounded=zeros(nt,1);
    parfor ii=1:nt
        [ASLLEQLBounded(ii),lossLCLBounded_ASLL(ii),lossUCLBounded_ASLL(ii)]=ASLLexpectedloss(x_e,w_noise,transparBounded,xl,ASLLthetalBounded,xctest(ii,:),ASLLbetalBounded,ASLLtao2lBounded,ASLLirReslBounded,ASLLirxlBounded...
            ,x,ASLLthetaBounded,ASLLrouBounded,ASLLbetaBounded,ASLLtao2Bounded,ASLLirResBounded,ASLLirxBounded);
    end
    MAREBoundedASLL=median(abs(qBoundedtest-ASLLEQLBounded)./qBoundedtest);
    MRCILBoundedASLL=median((lossUCLBounded_ASLL-lossLCLBounded_ASLL)./qBoundedtest);
    coverBoundedASLL=mean((qBoundedtest>=lossLCLBounded_ASLL).*(qBoundedtest<=lossUCLBounded_ASLL)*100);
    ARE95BoundedASLL=prctile(abs(qBoundedtest-ASLLEQLBounded)./qBoundedtest,95);
    RCIL95BoundedASLL=prctile((lossUCLBounded_ASLL-lossLCLBounded_ASLL)./qBoundedtest,95);

    lossLCLBounded_AGPL=zeros(nt,1);
    lossUCLBounded_AGPL=zeros(nt,1);
    AGPLEQLBounded=zeros(nt,1);
    parfor ii=1:nt
        [AGPLEQLBounded(ii),lossLCLBounded_AGPL(ii),lossUCLBounded_AGPL(ii)]=AGPLexpectedloss(x_e,w_noise,xl,AGPLthetalBounded,xctest(ii,:),AGPLbetalBounded,AGPLtao2lBounded,AGPLirReslBounded,AGPLirxlBounded,...
            x,AGPLthetaBounded,AGPLrouBounded,AGPLbetaBounded,AGPLtao2Bounded,AGPLirResBounded,AGPLirxBounded);
    end
    MAREBoundedAGPL=median(abs(qBoundedtest-AGPLEQLBounded)./qBoundedtest);
    MRCILBoundedAGPL=median((lossUCLBounded_AGPL-lossLCLBounded_AGPL)./qBoundedtest);
    coverBoundedAGPL=mean((qBoundedtest>=lossLCLBounded_AGPL).*(qBoundedtest<=lossUCLBounded_AGPL)*100);
    ARE95BoundedAGPL=prctile(abs(qBoundedtest-AGPLEQLBounded)./qBoundedtest,95);
    RCIL95BoundedAGPL=prctile((lossUCLBounded_AGPL-lossLCLBounded_AGPL)./qBoundedtest,95);

    lossLCLBounded_ALL=zeros(nt,1);
    lossUCLBounded_ALL=zeros(nt,1);
    ALLEQLBounded=zeros(nt,1);
    parfor ii=1:nt
        [ALLEQLBounded(ii),lossLCLBounded_ALL(ii),lossUCLBounded_ALL(ii)]=ALLexpectedloss(x_e,w_noise,xl,ALLthetalBounded,xctest(ii,:),ALLbetalBounded,ALLtao2lBounded,ALLirReslBounded,ALLirxlBounded,...
            x,ALLthetaBounded,ALLrouBounded,ALLbetaBounded,ALLtao2Bounded,ALLirResBounded,ALLirxBounded);
    end
    MAREBoundedALL=median(abs(qBoundedtest-ALLEQLBounded)./qBoundedtest);
    MRCILBoundedALL=median((lossUCLBounded_ALL-lossLCLBounded_ALL)./qBoundedtest);
    coverBoundedALL=mean((qBoundedtest>=lossLCLBounded_ALL).*(qBoundedtest<=lossUCLBounded_ALL)*100);
    ARE95BoundedALL=prctile(abs(qBoundedtest-ALLEQLBounded)./qBoundedtest,95);
    RCIL95BoundedALL=prctile((lossUCLBounded_ALL-lossLCLBounded_ALL)./qBoundedtest,95);

    MAREBounded(k,:)=[MAREBoundedASLL,MAREBoundedALL,MAREBoundedAGPL];
    MRCILBounded(k,:)=[MRCILBoundedASLL,MRCILBoundedALL,MRCILBoundedAGPL];
    coverBounded(k,:)=[coverBoundedASLL,coverBoundedALL,coverBoundedAGPL];
    ARE95Bounded(k,:)=[ARE95BoundedASLL,ARE95BoundedALL,ARE95BoundedAGPL];
    RCIL95Bounded(k,:)=[RCIL95BoundedASLL,RCIL95BoundedALL,RCIL95BoundedAGPL];
    %% get the summaries of the robust optimization performance of the ASLLGP, ALL, and AGPL models.
    expectedlossASLLBounded=zeros(nx,nx);
    parfor i=1:nx
        for j=1:nx
            [expectedlossASLLBounded(i,j),~,~]=ASLLexpectedloss(points,w_noise,transparBounded,xl,ASLLthetalBounded,[xc1(i,j),xc2(i,j)],ASLLbetalBounded,ASLLtao2lBounded,ASLLirReslBounded,ASLLirxlBounded...
                ,x,ASLLthetaBounded,ASLLrouBounded,ASLLbetaBounded,ASLLtao2Bounded,ASLLirResBounded,ASLLirxBounded);
        end
    end
    [indexASLL1Bounded,indexASLL2Bounded]=find(expectedlossASLLBounded==min(min(expectedlossASLLBounded)));xcminASLLBounded(k,:)=[xc1(indexASLL1Bounded,indexASLL2Bounded),xc2(indexASLL1Bounded,indexASLL2Bounded)];
minASLLEQLBounded(k)=qgridBounded(indexASLL1Bounded,indexASLL2Bounded);xcminASLLBoundedind(k,:)=[indexASLL1Bounded,indexASLL2Bounded];
    expectedlossAGPLBounded=zeros(nx,nx);
    parfor i=1:nx
        for j=1:nx
            [expectedlossAGPLBounded(i,j),~,~]=AGPLexpectedloss(points,w_noise,xl,AGPLthetalBounded,[xc1(i,j),xc2(i,j)],AGPLbetalBounded,AGPLtao2lBounded,AGPLirReslBounded,AGPLirxlBounded,...
                x,AGPLthetaBounded,AGPLrouBounded,AGPLbetaBounded,AGPLtao2Bounded,AGPLirResBounded,AGPLirxBounded);
        end
    end
    [indexAGPL1Bounded,indexAGPL2Bounded]=find(expectedlossAGPLBounded==min(min(expectedlossAGPLBounded)));xcminAGPLBounded(k,:)=[xc1(indexAGPL1Bounded,indexAGPL2Bounded),xc2(indexAGPL1Bounded,indexAGPL2Bounded)];
minAGPLEQLBounded(k)=qgridBounded(indexAGPL1Bounded,indexAGPL2Bounded);xcminAGPLBoundedind(k,:)=[indexAGPL1Bounded,indexAGPL2Bounded];
    expectedlossALLBounded=zeros(nx,nx);
    parfor i=1:nx
        for j=1:nx
            [ expectedlossALLBounded(i,j),~,~]=ALLexpectedloss(points,w_noise,xl,ALLthetalBounded,[xc1(i,j),xc2(i,j)],ALLbetalBounded,ALLtao2lBounded,ALLirReslBounded,ALLirxlBounded,...
                x,ALLthetaBounded,ALLrouBounded,ALLbetaBounded,ALLtao2Bounded,ALLirResBounded,ALLirxBounded);
        end
    end
    [indexALL1Bounded,indexALL2Bounded]=find(expectedlossALLBounded==min(min(expectedlossALLBounded)));xcminALLBounded(k,:)=[xc1(indexALL1Bounded,indexALL2Bounded),xc2(indexALL1Bounded,indexALL2Bounded)];
minALLEQLBounded(k)=qgridBounded(indexALL1Bounded,indexALL2Bounded);xcminALLBoundedind(k,:)=[indexALL1Bounded,indexALL2Bounded];
end

disp('Minimum true EQL');
minqBounded=min(min(qgridBounded))
disp('Mean of true EQL');
expqBounded=trapz(trapz(qgridBounded))*(0.005^2)
disp('Table S3.1');
disp('Sample means for MARE given by ASLLGP, ALL, AGPL models');
mMAREBounded=mean(MAREBounded)*100
disp('Standard errors for MARE given by ASLLGP, ALL, AGPL models');
stderrMAREBounded=std(MAREBounded)./sqrt(ndesign)*100
disp('Paired t-statistics for MARE given by ASLLGP, ALL, AGPL models');
[~,~,~,s11Bounded]=ttest(MAREBounded(:,2)-MAREBounded(:,1));[~,~,~,s12Bounded]=ttest(MAREBounded(:,3)-MAREBounded(:,1));
tMAREBounded=[0,s11Bounded.tstat,s12Bounded.tstat]
disp('Sample means for ARE_0.95 given by ASLLGP, ALL, AGPL models');
mARE95Bounded=mean(ARE95Bounded)*100
disp('Standard errors for ARE_0.95 given by ASLLGP, ALL, AGPL models');
stderrARE95Bounded=std(ARE95Bounded)./sqrt(ndesign)*100
disp('Paired t-statistics for ARE_0.95 given by ASLLGP, ALL, AGPL models');
[~,~,~,s41Bounded]=ttest(ARE95Bounded(:,2)-ARE95Bounded(:,1));[~,~,~,s42Bounded]=ttest(ARE95Bounded(:,3)-ARE95Bounded(:,1));
tARE95Bounded=[0,s41Bounded.tstat,s42Bounded.tstat]
disp('Sample means for MRCIL given by ASLLGP, ALL, AGPL models');
mMRCILBounded=mean(MRCILBounded)*100
disp('Standard errors for MRCIL given by ASLLGP, ALL, AGPL models');
stderrMRCILBounded=std(MRCILBounded)./sqrt(ndesign)*100
disp('Paired t-statistics for MRCIL given by ASLLGP, ALL, AGPL models');
[~,~,~,s21Bounded]=ttest(MRCILBounded(:,2)-MRCILBounded(:,1));[~,~,~,s22Bounded]=ttest(MRCILBounded(:,3)-MRCILBounded(:,1));
tMRCILBounded=[0,s21Bounded.tstat,s22Bounded.tstat]
disp('Sample means for RCIL_0.95 given by ASLLGP, ALL, AGPL models');
mRCIL95Bounded=mean(RCIL95Bounded)*100
disp('Standard errors for RCIL_0.95 given by ASLLGP, ALL, AGPL models');
stderrRCIL95Bounded=std(RCIL95Bounded)./sqrt(ndesign)*100
disp('Paired t-statistics for RCIL_0.95 given by ASLLGP, ALL, AGPL models');
[~,~,~,s51Bounded]=ttest(RCIL95Bounded(:,2)-RCIL95Bounded(:,1));[~,~,~,s52Bounded]=ttest(RCIL95Bounded(:,3)-RCIL95Bounded(:,1));
tRCIL95Bounded=[0,s51Bounded.tstat,s52Bounded.tstat]
disp('Sample means for EC-95 given by ASLLGP, ALL, AGPL models');
mcoverBounded=mean(coverBounded)
disp('Standard errors for EC-95 given by ASLLGP, ALL, AGPL models');
stderrcoverBounded=std(coverBounded)./sqrt(ndesign)
disp('Paired t-statistics for EC-95 given by ASLLGP, ALL, AGPL models');
[~,~,~,s31Bounded]=ttest(coverBounded(:,2)-coverBounded(:,1));[~,~,~,s32Bounded]=ttest(coverBounded(:,3)-coverBounded(:,1));
tcoverBounded=[0,s31Bounded.tstat,s32Bounded.tstat]
%% get the summaries of the five EQL prediction performance measures for the single-fidelity GP models and the two sample t-statistics.
LallonefiBounded=zeros(nonefi,ndesign);Ld_allonefiBounded=zeros(nonefi,ndesign);
for k=1:ndesign
    for i=1:nonefi
        [LallonefiBounded(i,k),Ld_allonefiBounded(i,k)]=lossPiezo(yd_allonefi(i,k),ym_allonefi(i,k),b);
    end
end

MAREBoundedMILE=zeros(ndesign,3);
MRCILBoundedMILE=zeros(ndesign,3);
coverBoundedMILE=zeros(ndesign,3);
ARE95BoundedMILE=zeros(ndesign,3);
RCIL95BoundedMILE=zeros(ndesign,3);
transpar_BoundedMILE=zeros(ndesign,2);

xcminSLLBounded_MILE=zeros(ndesign,2);xcminLLBounded_MILE=zeros(ndesign,2);xcminGPLBounded_MILE=zeros(ndesign,2);
minLLEQLBounded_MILE=zeros(ndesign,1);minSLLEQLBounded_MILE=zeros(ndesign,1);minGPLEQLBounded_MILE=zeros(ndesign,1);
xcminSLLBounded_MILEind=zeros(ndesign,2);xcminLLBounded_MILEind=zeros(ndesign,2);xcminGPLBounded_MILEind=zeros(ndesign,2);

para_SLLBounded=zeros(ndesign, 6);SLLbetaallBounded=zeros(ndesign,1);SLLtao2allBounded=zeros(ndesign,1);SLLirResallBounded=zeros(nonefi,ndesign);
para_GPLBounded=zeros(ndesign, 4);GPLbetalallBounded=zeros(ndesign,1);GPLtao2lallBounded=zeros(ndesign,1);GPLirReslallBounded=zeros(nonefi,ndesign);
para_LLBounded=zeros(ndesign, 4);LLbetalallBounded=zeros(ndesign,1);LLtao2lallBounded=zeros(ndesign,1);LLirReslallBounded=zeros(nonefi,ndesign);
for k=1:ndesign
    xonefi=reshape(xallonefi(k,:),d,[])';
    LonefiBounded=LallonefiBounded(:,k);
    [SLLthetaBounded_MILE,SLLfvallBounded_MILE,SLLbetaBounded_MILE,SLLtao2Bounded_MILE,SLLirResBounded_MILE,SLLirxBounded_MILE,~,~,transparBounded_MILE]=gpfitSLL(xonefi,LonefiBounded,ub);

    transpar_BoundedMILE(k,:)=transparBounded_MILE;
   para_SLLBounded(k,1:4)=SLLthetaBounded_MILE; para_SLLBounded(k,5:6)=[transparBounded_MILE(2),transparBounded_MILE(1)];
SLLbetaallBounded(k)=SLLbetaBounded_MILE;SLLtao2allBounded(k)=SLLtao2Bounded_MILE;SLLirResallBounded(:,k)=SLLirResBounded_MILE;

    lossLCLBounded_SLLMILE=zeros(nt,1);
    lossUCLBounded_SLLMILE=zeros(nt,1);

    SLLEQLBoundedMILE=zeros(nt,1);
    for ii=1:nt
        [SLLEQLBoundedMILE(ii),lossLCLBounded_SLLMILE(ii),lossUCLBounded_SLLMILE(ii)]=ASLLexpectedloss(x_e,w_noise,transparBounded_MILE,xonefi,SLLthetaBounded_MILE,xctest(ii,:),SLLbetaBounded_MILE,SLLtao2Bounded_MILE,SLLirResBounded_MILE,SLLirxBounded_MILE);

    end
    MARESLLBoundedMILE=median(abs(qBoundedtest-SLLEQLBoundedMILE)./qBoundedtest);
    MRCILSLLBoundedMILE=median((lossUCLBounded_SLLMILE-lossLCLBounded_SLLMILE)./qBoundedtest);
    coverSLLBoundedMILE=mean((qBoundedtest>=lossLCLBounded_SLLMILE).*(qBoundedtest<=lossUCLBounded_SLLMILE)*100);
    ARE95SLLBoundedMILE=prctile(abs(qBoundedtest-SLLEQLBoundedMILE)./qBoundedtest,95);
    RCIL95SLLBoundedMILE=prctile((lossUCLBounded_SLLMILE-lossLCLBounded_SLLMILE)./qBoundedtest,95);

    [GPLthetaBounded_MILE,GPLfvalBounded_MILE,GPLbetaBounded_MILE,GPLtao2Bounded_MILE,GPLirResBounded_MILE,GPLirxBounded_MILE,~,~]=gpfitGPL(xonefi,LonefiBounded);
    para_GPLBounded(k,:)=GPLthetaBounded_MILE;GPLbetalallBounded(k)=GPLbetaBounded_MILE;GPLtao2lallBounded(k)=GPLtao2Bounded_MILE;GPLirReslallBounded(:,k)=GPLirResBounded_MILE;
    lossLCLBounded_GPLMILE=zeros(nt,1);
    lossUCLBounded_GPLMILE=zeros(nt,1);
    GPLEQLBoundedMILE=zeros(nt,1);
    for ii=1:nt
        [GPLEQLBoundedMILE(ii),lossLCLBounded_GPLMILE(ii),lossUCLBounded_GPLMILE(ii)]=AGPLexpectedloss(x_e,w_noise,xonefi,GPLthetaBounded_MILE,xctest(ii,:),GPLbetaBounded_MILE,GPLtao2Bounded_MILE,GPLirResBounded_MILE,GPLirxBounded_MILE);
    end
    MAREGPLBoundedMILE=median(abs(qBoundedtest-GPLEQLBoundedMILE)./qBoundedtest);
    MRCILGPLBoundedMILE=median((lossUCLBounded_GPLMILE-lossLCLBounded_GPLMILE)./qBoundedtest);
    coverGPLBoundedMILE=mean((qBoundedtest>=lossLCLBounded_GPLMILE).*(qBoundedtest<=lossUCLBounded_GPLMILE)*100);
    ARE95GPLBoundedMILE=prctile(abs(qBoundedtest-GPLEQLBoundedMILE)./qBoundedtest,95);
    RCIL95GPLBoundedMILE=prctile((lossUCLBounded_GPLMILE-lossLCLBounded_GPLMILE)./qBoundedtest,95);

    [LLthetaBounded_MILE,LLfvalBounded_MILE,LLbetaBounded_MILE,LLtao2Bounded_MILE,LLirResBounded_MILE,LLirxBounded_MILE,~,~]=gpfitLL(xonefi,LonefiBounded);
    para_LLBounded(k,:)=LLthetaBounded_MILE;LLbetalallBounded(k)=LLbetaBounded_MILE;LLtao2lallBounded(k)=LLtao2Bounded_MILE;LLirReslallBounded(:,k)=LLirResBounded_MILE;

    lossLCLBounded_LLMILE=zeros(nt,1);
    lossUCLBounded_LLMILE=zeros(nt,1);
    LLEQLBoundedMILE=zeros(nt,1);
    for ii=1:nt
        [LLEQLBoundedMILE(ii),lossLCLBounded_LLMILE(ii),lossUCLBounded_LLMILE(ii)]=ALLexpectedloss(x_e,w_noise,xonefi,LLthetaBounded_MILE,xctest(ii,:),LLbetaBounded_MILE,LLtao2Bounded_MILE,LLirResBounded_MILE,LLirxBounded_MILE);
    end
    MARELLBoundedMILE=median(abs(qBoundedtest-LLEQLBoundedMILE)./qBoundedtest);
    MRCILLLBoundedMILE=median((lossUCLBounded_LLMILE-lossLCLBounded_LLMILE)./qBoundedtest);
    coverLLBoundedMILE=mean((qBoundedtest>=lossLCLBounded_LLMILE).*(qBoundedtest<=lossUCLBounded_LLMILE)*100);
    ARE95LLBoundedMILE=prctile(abs(qBoundedtest-LLEQLBoundedMILE)./qBoundedtest,95);
    RCIL95LLBoundedMILE=prctile((lossUCLBounded_LLMILE-lossLCLBounded_LLMILE)./qBoundedtest,95);

    MAREBoundedMILE(k,:)=[MARESLLBoundedMILE,MARELLBoundedMILE,MAREGPLBoundedMILE];
    MRCILBoundedMILE(k,:)=[MRCILSLLBoundedMILE,MRCILLLBoundedMILE,MRCILGPLBoundedMILE];
    coverBoundedMILE(k,:)=[coverSLLBoundedMILE,coverLLBoundedMILE,coverGPLBoundedMILE];
    ARE95BoundedMILE(k,:)=[ARE95SLLBoundedMILE,ARE95LLBoundedMILE,ARE95GPLBoundedMILE];
    RCIL95BoundedMILE(k,:)=[RCIL95SLLBoundedMILE,RCIL95LLBoundedMILE,RCIL95GPLBoundedMILE];

    %% get the summaries of the robust optimization performance of the single-fidelity GP models
    expectedlossSLLBounded_MILE=zeros(nx,nx);
    parfor i=1:nx
        for j=1:nx
            [expectedlossSLLBounded_MILE(i,j),~,~]=ASLLexpectedloss(points,w_noise,transparBounded_MILE,xonefi,SLLthetaBounded_MILE,[xc1(i,j),xc2(i,j)],SLLbetaBounded_MILE,SLLtao2Bounded_MILE,SLLirResBounded_MILE,SLLirxBounded_MILE);
        end
    end
    [indexSLL1Bounded_MILE,indexSLL2Bounded_MILE]=find(expectedlossSLLBounded_MILE==min(min(expectedlossSLLBounded_MILE)));xcminSLLBounded_MILE(k,:)=[xc1(indexSLL1Bounded_MILE,indexSLL2Bounded_MILE),xc2(indexSLL1Bounded_MILE,indexSLL2Bounded_MILE)];
minSLLEQLBounded_MILE(k)=qgridBounded(indexSLL1Bounded_MILE,indexSLL2Bounded_MILE);xcminSLLBounded_MILEind(k,:)=[indexSLL1Bounded_MILE,indexSLL2Bounded_MILE];

    expectedlossGPLBounded_MILE=zeros(nx,nx);
    parfor i=1:nx
        for j=1:nx
            [expectedlossGPLBounded_MILE(i,j),~,~]=AGPLexpectedloss(points,w_noise,xonefi,GPLthetaBounded_MILE,[xc1(i,j),xc2(i,j)],GPLbetaBounded_MILE,GPLtao2Bounded_MILE,GPLirResBounded_MILE,GPLirxBounded_MILE);
        end
    end
    [indexGPL1Bounded_MILE,indexGPL2Bounded_MILE]=find(expectedlossGPLBounded_MILE==min(min(expectedlossGPLBounded_MILE)));xcminGPLBounded_MILE(k,:)=[xc1(indexGPL1Bounded_MILE,indexGPL2Bounded_MILE),xc2(indexGPL1Bounded_MILE,indexGPL2Bounded_MILE)];
minGPLEQLBounded_MILE(k)=qgridBounded(indexGPL1Bounded_MILE,indexGPL2Bounded_MILE);xcminGPLBounded_MILEind(k,:)=[indexGPL1Bounded_MILE,indexGPL2Bounded_MILE];

    expectedlossLLBounded_MILE=zeros(nx,nx);
    parfor i=1:nx
        for j=1:nx
            [ expectedlossLLBounded_MILE(i,j),~,~]=ALLexpectedloss(points,w_noise,xonefi,LLthetaBounded_MILE,[xc1(i,j),xc2(i,j)],LLbetaBounded_MILE,LLtao2Bounded_MILE,LLirResBounded_MILE,LLirxBounded_MILE);
        end
    end
    [indexLL1Bounded_MILE,indexLL2Bounded_MILE]=find(expectedlossLLBounded_MILE==min(min(expectedlossLLBounded_MILE)));xcminLLBounded_MILE(k,:)=[xc1(indexLL1Bounded_MILE,indexLL2Bounded_MILE),xc2(indexLL1Bounded_MILE,indexLL2Bounded_MILE)];
minLLEQLBounded_MILE(k)=qgridBounded(indexLL1Bounded_MILE,indexLL2Bounded_MILE);xcminLLBounded_MILEind(k,:)=[indexLL1Bounded_MILE,indexLL2Bounded_MILE];
end

disp('Table S3.2');
disp('Sample means for MARE given by shifted log loss GP, lognormal loss, and GP for loss models');
mMAREBoundedMILE=mean(MAREBoundedMILE)*100
disp('Standard errors for MARE given by shifted log loss GP, lognormal loss, and GP for loss models');
stderrMAREBoundedMILE=std(MAREBoundedMILE)./sqrt(ndesign)*100
disp('Two sample t-statistics for MARE given by shifted log loss GP, lognormal loss, and GP for loss models');
[~,~,~,stats11BoundedMILE]=ttest2(MAREBoundedMILE(:,1),MAREBounded(:,1),'Vartype','unequal');[~,~,~,stats12BoundedMILE]=ttest2(MAREBoundedMILE(:,2),MAREBounded(:,1),'Vartype','unequal');[~,~,~,stats13BoundedMILE]=ttest2(MAREBoundedMILE(:,3),MAREBounded(:,1),'Vartype','unequal');
tMAREdifffiBoundedMILE=[stats11BoundedMILE.tstat,stats12BoundedMILE.tstat,stats13BoundedMILE.tstat]
disp('Sample means for ARE_0.95 given by shifted log loss GP, lognormal loss, and GP for loss models');
mARE95BoundedMILE=mean(ARE95BoundedMILE)*100
disp('Standard errors for ARE_0.95 given by shifted log loss GP, lognormal loss, and GP for loss models');
stderrARE95BoundedMILE=std(ARE95BoundedMILE)./sqrt(ndesign)*100
[~,~,~,stats41BoundedMILE]=ttest2(ARE95BoundedMILE(:,1),ARE95Bounded(:,1),'Vartype','unequal');[~,~,~,stats42BoundedMILE]=ttest2(ARE95BoundedMILE(:,2),ARE95Bounded(:,1),'Vartype','unequal');[~,~,~,stats43BoundedMILE]=ttest2(ARE95BoundedMILE(:,3),ARE95Bounded(:,1),'Vartype','unequal');
disp('Two sample t-statistics for ARE_0.95 given by shifted log loss GP, lognormal loss, and GP for loss models');
tARE95difffiBoundedMILE=[stats41BoundedMILE.tstat,stats42BoundedMILE.tstat,stats43BoundedMILE.tstat]
disp('Sample means for MRCIL given by shifted log loss GP, lognormal loss, and GP for loss models');
mMRCILBoundedMILE=mean(MRCILBoundedMILE)*100
disp('Standard errors for MRCIL given by shifted log loss GP, lognormal loss, and GP for loss models');
stderrMRCILBoundedMILE=std(MRCILBoundedMILE)./sqrt(ndesign)*100
[~,~,~,stats21BoundedMILE]=ttest2(MRCILBoundedMILE(:,1),MRCILBounded(:,1),'Vartype','unequal');[~,~,~,stats22BoundedMILE]=ttest2(MRCILBoundedMILE(:,2),MRCILBounded(:,1),'Vartype','unequal');[~,~,~,stats23BoundedMILE]=ttest2(MRCILBoundedMILE(:,3),MRCILBounded(:,1),'Vartype','unequal');
disp('Two sample t-statistics for MRCIL given by shifted log loss GP, lognormal loss, and GP for loss models');
tMRCILdifffiBoundedMILE=[stats21BoundedMILE.tstat,stats22BoundedMILE.tstat,stats23BoundedMILE.tstat]
disp('Sample means for RCIL_0.95 given by shifted log loss GP, lognormal loss, and GP for loss models');
mRCIL95BoundedMILE=mean(RCIL95BoundedMILE)*100
disp('Standard errors for RCIL_0.95 given by shifted log loss GP, lognormal loss, and GP for loss models');
stderrRCIL95BoundedMILE=std(RCIL95BoundedMILE)./sqrt(ndesign)*100
disp('Two sample t-statistics for RCIL_0.95 given by shifted log loss GP, lognormal loss, and GP for loss models');
[~,~,~,stats51BoundedMILE]=ttest2(RCIL95BoundedMILE(:,1),RCIL95Bounded(:,1),'Vartype','unequal');[~,~,~,stats52BoundedMILE]=ttest2(RCIL95BoundedMILE(:,2),RCIL95Bounded(:,1),'Vartype','unequal');[~,~,~,stats53BoundedMILE]=ttest2(RCIL95BoundedMILE(:,3),RCIL95Bounded(:,1),'Vartype','unequal');
tRCIL95difffiBoundedMILE=[stats51BoundedMILE.tstat,stats52BoundedMILE.tstat,stats53BoundedMILE.tstat]
disp('Sample means for EC-95 given by shifted log loss GP, lognormal loss, and GP for loss models');
mcoverBounded_MILE=mean(coverBoundedMILE)
disp('Standard errors for EC-95 given by shifted log loss GP, lognormal loss, and GP for loss models');
stderrcoverBounded_MILE=std(coverBoundedMILE)./sqrt(ndesign)
disp('Two sample t-statistics for EC-95 given by shifted log loss GP, lognormal loss, and GP for loss models');
[~,~,~,stats31BoundedMILE]=ttest2(coverBoundedMILE(:,1),coverBounded(:,1),'Vartype','unequal');[~,~,~,stats32BoundedMILE]=ttest2(coverBoundedMILE(:,2),coverBounded(:,1),'Vartype','unequal');[~,~,~,stats33BoundedMILE]=ttest2(coverBoundedMILE(:,3),coverBounded(:,1),'Vartype','unequal');
tcoverdifffiBoundedMILE=[stats31BoundedMILE.tstat,stats32BoundedMILE.tstat,stats33BoundedMILE.tstat]
disp('Table S3.3');
disp('Sample means for the true EQL at the estimated optimal control-factor setting given by ASLLGP, ALL, AGPL models');
minEQLBounded=[minASLLEQLBounded,minALLEQLBounded,minAGPLEQLBounded];
MEQLBounded=mean(minEQLBounded)
disp('Standard errors for the true EQL at the estimated optimal control-factor setting given by ASLLGP, ALL, AGPL models');
STDEQLBounded=std(minEQLBounded)./sqrt(ndesign)
disp('Paired t-statistics for the true EQL at the estimated optimal control-factor setting given by ASLLGP, ALL, AGPL models');
[~,~,~,s1Bounded]=ttest(minEQLBounded(:,2)-minEQLBounded(:,1));[~,~,~,s2Bounded]=ttest(minEQLBounded(:,3)-minEQLBounded(:,1));
tEQLBounded=[0,s1Bounded.tstat,s2Bounded.tstat]
disp('Sample means for the true EQL at the estimated optimal control-factor setting given by shifted log loss GP, lognormal loss, and GP for loss models');
minEQLBounded_MILE=[minSLLEQLBounded_MILE,minLLEQLBounded_MILE,minGPLEQLBounded_MILE];
MEQLBounded_MILE=mean(minEQLBounded_MILE)
disp('Standard errors for the true EQL at the estimated optimal control-factor setting given by shifted log loss GP, lognormal loss, and GP for loss models');
STDEQLBounded_MILE=std(minEQLBounded_MILE)./sqrt(ndesign)
disp('Two sample t-statistics for the true EQL at the estimated optimal control-factor setting given by shifted log loss GP, lognormal loss, and GP for loss models');
[~,~,~,s0Bounded_MILE]=ttest2(minEQLBounded_MILE(:,1),minEQLBounded(:,1),'Vartype','unequal');[~,~,~,s1Bounded_MILE]=ttest2(minEQLBounded_MILE(:,2),minEQLBounded(:,1),'Vartype','unequal');[~,~,~,s2Bounded_MILE]=ttest2(minEQLBounded_MILE(:,3),minEQLBounded(:,1),'Vartype','unequal');
tEQLBounded_MILE=[s0Bounded_MILE.tstat,s1Bounded_MILE.tstat,s2Bounded_MILE.tstat]

tMAREdf_MILE=[stats11BoundedMILE.df,stats12BoundedMILE.df,stats13BoundedMILE.df];
tMRCILdf_MILE=[stats21BoundedMILE.df,stats22BoundedMILE.df,stats23BoundedMILE.df];
tARE95df_MILE=[stats41BoundedMILE.df,stats42BoundedMILE.df,stats43BoundedMILE.df];
tRCIL95df_MILE=[stats51BoundedMILE.df,stats52BoundedMILE.df,stats53BoundedMILE.df];
tcoverdf_MILE=[stats31BoundedMILE.df,stats32BoundedMILE.df,stats33BoundedMILE.df];
mindf_MILE=min([min(tMAREdf_MILE),min(tMRCILdf_MILE),min(tARE95df_MILE),min(tRCIL95df_MILE),min(tcoverdf_MILE)]);
maxdf_MILE=max([max(tMAREdf_MILE),max(tMRCILdf_MILE),max(tARE95df_MILE),max(tRCIL95df_MILE),max(tcoverdf_MILE)]);
tEQLdf_MILE=[s0Bounded_MILE.df,s1Bounded_MILE.df,s2Bounded_MILE.df];
