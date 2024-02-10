function [L,L_d,L_m]=lossPiezo(yd,ym)
%Inputs yd:tip deflection of the beam;ym: mass of the beam
%Outputs L:total quality cost; L_d:quality cost incurred by the tip deflection;L_m:quality cost incurred by the mass of the beam.
tm=7.5e-5;td=-8e-5;
k_1=3e3;k_2=1;
L_m=k_2*(ym-tm).*(1-((ym-tm)<0));
L_d=k_1*(yd-td).^2; 
L=L_d+L_m;
end
