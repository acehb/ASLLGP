function q=PiezoTrueEQL(x_c,x_e,w_noise)
%Inputs xc:control factor at which estimation of the true expected loss is concerned; x_e:quadrature nodes; w_noise:quadrature weights times xi(x_e)
%Outputs q:true EQL

Nx_e=size(x_e,1);
Lent=100e-3+x_c(1)*100e-3;
Heit=0.6e-3+x_c(2)*0.6e-3;
d31t=1.5e-11+x_e(:,1)*1.5e-11;
vt=90+x_e(:,2)*20;
ydtemp=zeros(Nx_e,1);ymtemp=zeros(Nx_e,1);
losstemps=zeros(Nx_e,1);losstemp=zeros(Nx_e,1);
losstempm=zeros(Nx_e,1);
for k=1:Nx_e
     [ydtemp(k),ymtemp(k)]=PiezoelectricActuator(Lent,Heit,d31t(k),vt(k),2);
    [losstemp(k),losstemps(k),losstempm(k)]=lossPiezo(ydtemp(k),ymtemp(k));
end
loss_1=w_noise*losstemps;
q=loss_1+losstempm(1);

end