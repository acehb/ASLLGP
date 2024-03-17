clear;
load('Designs for the piezoelectric actuator example.mat');
load('400 test points for the piezoelectric actuator example.mat');
load('Piezo True EQL on grid points.mat');
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

