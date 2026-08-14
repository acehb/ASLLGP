%% generate 400 test points for Example 2 and Example 3 and compute the true EQL(based on the unbounded and bounded loss functions) at those points
clear all;
load('Piezo True EQL on grid points.mat');
nt=400;
%b=round(maxl/100,3,'significant');
qtest=zeros(nt,1);
qBoundedtest=zeros(nt,1);
xctest=lhsdesign(nt,2);

ydtest_all=zeros(nt,Nx_e); 
ymtest_all=zeros(nt,Nx_e); 
HFtime=zeros(nt,Nx_e);
parfor i=1:nt
    
    Lent=100e-3+xctest(i,1)*100e-3;
    Heit=0.6e-3+xctest(i,2)*0.6e-3;
      d31t=1.5e-11+x_e(:,1)*1.5e-11;
    vt=90+x_e(:,2)*20;
    
      ydtest=zeros(Nx_e,1);
    ymtest=zeros(Nx_e,1);
    for k=1:Nx_e
        tic
        [ydtest(k),ymtest(k)]=PiezoelectricActuator(Lent,Heit,d31t(k),vt(k),2);
    HFtime(i,k)=toc;
    end
    ydtest_all(i,:)=ydtest';  
    ymtest_all(i,:)=ymtest';  
end
%save('400 test points for the piezoelectric actuator example yd ym.mat');
parfor i=1:nt
    losstesttemps=zeros(Nx_e,1);
    losstesttemp=zeros(Nx_e,1);
    losstesttempm=zeros(Nx_e,1);
    
    for k=1:Nx_e
        [losstesttemp(k),losstesttemps(k),losstesttempm(k)]=lossPiezo(ydtest_all(i,k),ymtest_all(i,k));
    end
    losstest_1=w_noise*losstesttemps;
    qtest(i)=losstest_1+losstesttempm(1);
end  
parfor i=1:nt
    losstesttemps=zeros(Nx_e,1);
    losstesttemp=zeros(Nx_e,1);
    losstesttempm=zeros(Nx_e,1);
 
    for k=1:Nx_e
        [losstesttemp(k),losstesttemps(k),losstesttempm(k)]=lossPiezo(ydtest_all(i,k),ymtest_all(i,k),b);
    end
    qBoundedtest(i)=w_noise*losstesttemp;
end
save('400 test points for the piezoelectric actuator example.mat', 'qtest', 'qBoundedtest', 'w_noise', 'nt', 'x_e', 'xctest', 'ri_noise');