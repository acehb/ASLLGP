function [L,L_s,L_m]=lossBar(ys,ym)
%Inputs ys:the maximum von Mises stress in the bar;ym: the mass of the small bar.
%Outputs L:total quality cost; L_s:quality cost incurred by ys;L_m:quality cost incurred by ym.
 
sigmas=32.53;means=273.76;k_1=1.5;
tm=8.25;k_2=1;

L_m=k_2*(ym-tm).*(1-((ym-tm)<0));
L_s=k_1*((ys-means).*normcdf(ys,means,sigmas)+sigmas./sqrt(2.*pi).*exp(-((ys-means)./sigmas).^2./2));       
        L=L_s+L_m;
end
