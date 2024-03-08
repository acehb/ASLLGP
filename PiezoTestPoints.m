% generate 400 test points for Example 2 and compute the true EQL at those points
clear;
[ri_noise,ci_noise]=lgwt(6,0,1);
w_noise=zeros(1,length(ri_noise)*length(ri_noise));
k=1;
sigma_noise=1/6;
for i=1:length(ri_noise)
    for j=1:length(ri_noise)
        fx3 =normpdf(ri_noise(j),0.5,sigma_noise)/(normcdf(1,0.5,sigma_noise)-normcdf(0,0.5,sigma_noise));
        fx4 =normpdf(ri_noise(i),0.5,sigma_noise)/(normcdf(1,0.5,sigma_noise)-normcdf(0,0.5,sigma_noise));
        w_noise(k)=ci_noise(j)*ci_noise(i)*fx3*fx4;
        k=k+1;
    end
    
end
temp1=repmat(ri_noise,1,length(ri_noise));
x_e=[reshape(temp1,[],1),reshape(temp1',[],1)];
nt=400;Nx_e=k-1;
xctest=lhsdesign(nt,2);qtest=zeros(nt,1);
parfor i=1:nt
  [qtest(i)]=PiezoTrueEQL(xctest(i,:),x_e,w_noise);

end
save('400 test points for the piezoelectric actuator example.mat');
