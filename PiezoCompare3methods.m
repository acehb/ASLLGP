clear;
load('Designs for the piezoelectric actuator example.mat');
load('400 test points for the piezoelectric actuator example.mat');
load('Piezo True EQL on grid points.mat');
b=0.0013;
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

trans_ASLL=zeros(ndesign,4);
MARE=zeros(ndesign,3);
MRCIL=zeros(ndesign,3);
cover=zeros(ndesign,3);
ARE95=zeros(ndesign,3);
RCIL95=zeros(ndesign,3);
temp1=repmat(ri_noise,1,length(ri_noise));%
points=[reshape(temp1,[],1),reshape(temp1',[],1)];%

[indextrue1,indextrue2]=find(qgrid==min(min(qgrid)));xcmin=[xc1(indextrue1,indextrue2),xc2(indextrue1,indextrue2)];%
xcminASLL=zeros(ndesign,2);
xcminALL=zeros(ndesign,2);xcminAGPL=zeros(ndesign,2);
minALLEQL=zeros(ndesign,1);minASLLEQL=zeros(ndesign,1);minAGPLEQL=zeros(ndesign,1);
for k=1:ndesign
    %fit 3 models
    xl=reshape(xlall(k,:),d,[])';x=xl(1:m,:);%
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

   expectedlossAGPL=zeros(nx,nx);
    parfor i=1:nx
        for j=1:nx
        [expectedlossAGPL(i,j),~,~]=AGPLexpectedloss(points,w_noise,xl,AGPLthetal,[xc1(i,j),xc2(i,j)],AGPLbetal,AGPLtao2l,AGPLirResl,AGPLirxl,...
            x,AGPLtheta,AGPLrou,AGPLbeta,AGPLtao2,AGPLirRes,AGPLirx);
        end
    end
    [indexAGPL1,indexAGPL2]=find(expectedlossAGPL==min(min(expectedlossAGPL)));xcminAGPL(k,:)=[xc1(indexAGPL1,indexAGPL2),xc2(indexAGPL1,indexAGPL2)];

     expectedlossALL=zeros(nx,nx);
      parfor i=1:nx
        for j=1:nx
        [ expectedlossALL(i,j),~,~]=ALLexpectedloss(points,w_noise,xl,ALLthetal,[xc1(i,j),xc2(i,j)],ALLbetal,ALLtao2l,ALLirResl,ALLirxl,...
            x,ALLtheta,ALLrou,ALLbeta,ALLtao2,ALLirRes,ALLirx);
        end
      end
    [indexALL1,indexALL2]=find(expectedlossALL==min(min(expectedlossALL)));xcminALL(k,:)=[xc1(indexALL1,indexALL2),xc2(indexALL1,indexALL2)];
end
parfor k=1:ndesign
    minALLEQL(k)=PiezoTrueEQL(xcminALL(k,:),points,w_noise);
    minASLLEQL(k)=PiezoTrueEQL(xcminASLL(k,:),points,w_noise);
    minAGPLEQL(k)=PiezoTrueEQL(xcminAGPL(k,:),points,w_noise);
end

minq=min(min(qgrid))
%max(max(qgrid)),mean(mean(qgrid))
expq=trapz(trapz(qgrid))*(0.005^2)
minEQL=[minASLLEQL,minALLEQL,minAGPLEQL];
MEQL=mean(minEQL)
STDEQL=std(minEQL)./sqrt(ndesign)
[~,~,~,s1]=ttest(minEQL(:,2)-minEQL(:,1));[~,~,~,s2]=ttest(minEQL(:,3)-minEQL(:,1));
tEQL=[0,s1.tstat,s2.tstat]

% five performance measures and the paired sample t-statistics
mMARE=mean(MARE)*100
stderrMARE=std(MARE)./sqrt(ndesign)*100
[~,~,~,s11]=ttest(MARE(:,2)-MARE(:,1));[~,~,~,s12]=ttest(MARE(:,3)-MARE(:,1));
tMARE=[0,s11.tstat,s12.tstat]

mMRCIL=mean(MRCIL)*100
stderrMRCIL=std(MRCIL)./sqrt(ndesign)*100
[~,~,~,s21]=ttest(MRCIL(:,2)-MRCIL(:,1));[~,~,~,s22]=ttest(MRCIL(:,3)-MRCIL(:,1));
tMRCIL=[0,s21.tstat,s22.tstat]

mARE95=mean(ARE95)*100
stderrARE95=std(ARE95)./sqrt(ndesign)*100
[~,~,~,s41]=ttest(ARE95(:,2)-ARE95(:,1));[~,~,~,s42]=ttest(ARE95(:,3)-ARE95(:,1));
tARE95=[0,s41.tstat,s42.tstat]

mRCIL95=mean(RCIL95)*100
stderrRCIL95=std(RCIL95)./sqrt(ndesign)*100
[~,~,~,s51]=ttest(RCIL95(:,2)-RCIL95(:,1));[~,~,~,s52]=ttest(RCIL95(:,3)-RCIL95(:,1));
tRCIL95=[0,s51.tstat,s52.tstat]

mcover=mean(cover)
stderrcover=std(cover)./sqrt(ndesign)
[~,~,~,s31]=ttest(cover(:,2)-cover(:,1));[~,~,~,s32]=ttest(cover(:,3)-cover(:,1));
tcover=[0,s31.tstat,s32.tstat]
%%
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
transBounded_ASLL=zeros(ndesign,4);
MAREBounded=zeros(ndesign,3);
MRCILBounded=zeros(ndesign,3);
coverBounded=zeros(ndesign,3);
ARE95Bounded=zeros(ndesign,3);
RCIL95Bounded=zeros(ndesign,3);

for k=1:ndesign
    %fit 3 models
    xl=reshape(xlall(k,:),d,[])';x=xl(1:m,:);%
    LhBounded=LhallBounded(:,k);LlBounded=LlallBounded(:,k);LhlBounded=LlallBounded(1:m,k);
    [ASLLthetalBounded,ASLLfvallBounded,ASLLbetalBounded,ASLLtao2lBounded,ASLLirReslBounded,ASLLirxlBounded,transparlBounded]=gpfitASLL1level(xl,LlBounded,b);
    ASLLhlBounded=log(transparlBounded(2)*LhlBounded+transparlBounded(1));
    [ASLLthetaBounded,ASLLfvalBounded,ASLLrouBounded,ASLLbetaBounded,ASLLtao2Bounded,ASLLirResBounded,ASLLirxBounded,transparBounded]=gpfitASLL2level(x,ASLLhlBounded,LhBounded,b);
     transBounded_ASLL(k,1:2)=transparlBounded;  transBounded_ASLL(k,3:4)=transparBounded;  
    [AGPLthetalBounded,AGPLfvallBounded,AGPLbetalBounded,AGPLtao2lBounded,AGPLirReslBounded,AGPLirxlBounded]=gpfit1level(xl,LlBounded);
    [AGPLthetaBounded,AGPLfvalBounded,AGPLrouBounded,AGPLbetaBounded,AGPLtao2Bounded,AGPLirResBounded,AGPLirxBounded]=gpfit2level(x,LhlBounded,LhBounded);
        
    ALLLhBounded=log(LhBounded);ALLLlBounded=log(LlBounded); ALLhlBounded=ALLLlBounded(1:m);
    [ALLthetalBounded,ALLfvallBounded,ALLbetalBounded,ALLtao2lBounded,ALLirReslBounded,ALLirxlBounded]=gpfit1level(xl,ALLLlBounded);
    [ALLthetaBounded,ALLfvalBounded,ALLrouBounded,ALLbetaBounded,ALLtao2Bounded,ALLirResBounded,ALLirxBounded]=gpfit2level(x,ALLhlBounded,ALLLhBounded);
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

    end

% five performance measures and the paired sample t-statistics
mMAREBounded=mean(MAREBounded)*100
stderrMAREBounded=std(MAREBounded)./sqrt(ndesign)*100
[~,~,~,s11Bounded]=ttest(MAREBounded(:,2)-MAREBounded(:,1));[~,~,~,s12Bounded]=ttest(MAREBounded(:,3)-MAREBounded(:,1));
tMAREBounded=[0,s11Bounded.tstat,s12Bounded.tstat]

mMRCILBounded=mean(MRCILBounded)*100
stderrMRCILBounded=std(MRCILBounded)./sqrt(ndesign)*100
[~,~,~,s21Bounded]=ttest(MRCILBounded(:,2)-MRCILBounded(:,1));[~,~,~,s22Bounded]=ttest(MRCILBounded(:,3)-MRCILBounded(:,1));
tMRCILBounded=[0,s21Bounded.tstat,s22Bounded.tstat]

mARE95Bounded=mean(ARE95Bounded)*100
stderrARE95Bounded=std(ARE95Bounded)./sqrt(ndesign)*100
[~,~,~,s41Bounded]=ttest(ARE95Bounded(:,2)-ARE95Bounded(:,1));[~,~,~,s42Bounded]=ttest(ARE95Bounded(:,3)-ARE95Bounded(:,1));
tARE95Bounded=[0,s41Bounded.tstat,s42Bounded.tstat]

mRCIL95Bounded=mean(RCIL95Bounded)*100
stderrRCIL95Bounded=std(RCIL95Bounded)./sqrt(ndesign)*100
[~,~,~,s51Bounded]=ttest(RCIL95Bounded(:,2)-RCIL95Bounded(:,1));[~,~,~,s52Bounded]=ttest(RCIL95Bounded(:,3)-RCIL95Bounded(:,1));
tRCIL95Bounded=[0,s51Bounded.tstat,s52Bounded.tstat]

mcoverBounded=mean(coverBounded)
stderrcoverBounded=std(coverBounded)./sqrt(ndesign)
[~,~,~,s31Bounded]=ttest(coverBounded(:,2)-coverBounded(:,1));[~,~,~,s32Bounded]=ttest(coverBounded(:,3)-coverBounded(:,1));
tcoverBounded=[0,s31Bounded.tstat,s32Bounded.tstat]
