function [ys,ym]=ShaftInterp(stress,mass,x,fidelity)
%Inputs stress: the maximum Von Mises stress outputs of the HF and LF simulators on the 11^2 grid {0,0.1,бн,1}^2; mass:the mass of the small shaft output of both simulators on the 11^2 grid {0,0.1,бн,1}^2;x:input setting at which interpolated outputs is desired;
%fidelity can be 1 or 2, fidelity=1:interpolator for the low-fidelity simulator,fidelity=2:interpolator for the high-fidelity simulator.
%Outputs ys:the interpolated maximum von mises stress;ym: the interpolated mass of the small shaft.
[xcgrid,xegrid] = meshgrid(0:0.1:1);
ngrid=11;
ymtemp=reshape(mass(:,1),[],ngrid);
r=0.2;
ymgrid=ymtemp(1,:)-((r-0.02*2).^2+(r-0.02*2)*0.02*4+pi*(0.02.^2))*0.2*7800;
if(fidelity==1) %LF
    ylgrid=reshape(stress(:,1),[],ngrid);
    ys=interp2(xcgrid,xegrid,ylgrid,x(:,1),x(:,2));
end
if(fidelity==2) %HF
yhgrid=reshape(stress(:,2),[],ngrid);
ys=interp2(xcgrid,xegrid,yhgrid,x(:,1),x(:,2));
end
ym=interp1(xcgrid(1,:),ymgrid,x(:,1)); 

end
