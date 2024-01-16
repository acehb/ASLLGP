function [L,L_s,L_m]=lossShaft(ys,ym)
%Inputs ys:the maximum von Mises stress in the stepped shaft;ym: the mass of the small shaft.
%Outputs L:total quality cost; L_s:quality cost incurred by ys;L_m:quality cost incurred by ym.
 maxas=390;sigmas=(-170+maxas)/6;means=(170+maxas)/2;
 k_1=1.9;k_2=1;tm=9;%
L_m=k_2*(ym-tm).*(1-((ym-tm)<0));
L_s=k_1*((ys-means).*normcdf(ys,means,sigmas)+sigmas./sqrt(2.*pi).*exp(-((ys-means)./sigmas).^2./2));       
        L=L_s+L_m;
end
