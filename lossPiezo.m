function [L,L_d,L_m]=lossPiezo(yd,ym,b)
%Inputs yd:tip deflection of the beam;ym: mass of the beam;%b:scaling parameter of the bounded loss(omit for unbounded loss function).
%Outputs L:total quality cost; L_d:quality cost incurred by the tip deflection;L_m:quality cost incurred by the mass of the beam.
tm=7.25e-5;td=-7.8e-5;
k_1=3300;k_2=1;
L_m=k_2*(ym-tm).*(1-((ym-tm)<0));
L_d=k_1*(yd-td).^2;
if(nargin==2)
    L=L_d+L_m;
else
    L0=L_d+L_m;L=L0/b/(L0/b+1);
end
