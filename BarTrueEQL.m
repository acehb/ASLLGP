function q=BarTrueEQL(x_c,stress,mass,w,x_e)
%Inputs xc:control factor at which estimation of the true expected loss is concerned; x_e:quadrature nodes; w_noise:quadrature weights times xi(x_e);stress: the maximum Von Mises stress outputs of the HF and LF simulators on the 11^2 grid {0,0.1,бн,1}^2; mass:the mass of the small bar output of both simulators on the 11^2 grid {0,0.1,бн,1}^2
%Outputs q:true EQL
Nx_e=size(x_e,1);
ystemp=zeros(Nx_e,1);ymtemp=zeros(Nx_e,1);
losstemps=zeros(Nx_e,1);losstemp=zeros(Nx_e,1);
losstempm=zeros(Nx_e,1);
for k=1:Nx_e
    [ystemp(k),ymtemp(k)]=BarInterp(stress,mass,[x_c,x_e(k)],2);
    [losstemp(k),losstemps(k),losstempm(k)]=lossBar(ystemp(k),ymtemp(k));
end
loss_1=w*losstemps;
q=loss_1+losstempm(1);
end
