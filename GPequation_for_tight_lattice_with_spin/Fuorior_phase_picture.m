clear
profile on
V0=0;
gamma_P=1/6;gamma_nR=2*gamma_P;R=0.01;P0=gamma_nR*gamma_P/R;
scalet=1;scale=1.52*scalet;
dT=0.001;g_nR=0.01;Omega_F=4;g_P=0.005;
nx=20;
F=zeros(2*nx,1);
J_matrix=2*ones(1,nx);
phi=0.3*pi*ones(1,nx);
Driven_T=8.3;
Step=5000;Cord=40;Omega_length=400;site=nx/2;

P0_up=P0*linspace(0.6,1.6,20);
P0_down=P0*linspace(0.6,1.6,20);
Amplitude_T_Sx=zeros(20,20);
Amplitude_T_Sy=zeros(20,20);
Amplitude_T_Sz=zeros(20,20);
for i=1:20
    for j=1:20
    [Vpo,Pump,Coupling,Psi0] = Calculate_potential(nx,phi,V0,P0_up(i),P0_down(j),0.5*pi,0);
    [Result_matrix_up,Result_matrix_down,~,Time] = Iterative_program_Multiple_gpu(dT,Step,Cord,...
        scale*Vpo,scale*J_matrix,scale*Coupling,scalet*Pump,scalet*Driven_T,scale*F,scale*Omega_F,scale*g_P,scale*g_nR,scalet*R,scalet*gamma_P,scalet*gamma_nR,Psi0);
    Sz=(abs(Result_matrix_up).^2-abs(Result_matrix_down).^2)./(abs(Result_matrix_up).^2+abs(Result_matrix_down).^2);
    Psi_H=(Result_matrix_up+Result_matrix_down)/sqrt(2);
    Psi_V=(Result_matrix_up-Result_matrix_down)/sqrt(2);
    Sx=(abs(Psi_H).^2-abs(Psi_V).^2)./(abs(Psi_H).^2+abs(Psi_V).^2);
    Psi_D=(exp(1j*0.25*pi)*Result_matrix_up+exp(-1j*0.25*pi)*Result_matrix_down)/sqrt(2);
    Psi_A=(exp(1j*0.25*pi)*Result_matrix_up-exp(-1j*0.25*pi)*Result_matrix_down)/sqrt(2);
    Sy=(abs(Psi_D).^2-abs(Psi_A).^2)./(abs(Psi_D).^2+abs(Psi_A).^2);
    %[~,Sphi,Stheta]=cart2sph(Sx,Sy,Sz);
    Omega0=2*pi/Driven_T;
    
    omega=linspace(0.2*Omega0,0.8*Omega0,Omega_length);
    Amplitude=zeros(4,Omega_length);
    for k=1:Omega_length
        Exptime=exp(-1j*Time*omega(k));
        Amplitude(1,k)=sum((Sx(site,:))*Exptime)/Step;
        Amplitude(2,k)=sum((Sy(site,:))*Exptime)/Step;
        Amplitude(3,k)=sum((Sz(site,:))*Exptime)/Step;
    end
    Amplitude_T_Sx(i,j)=max(abs(Amplitude(1,:)).^2);
    Amplitude_T_Sy(i,j)=max(abs(Amplitude(2,:)).^2);
    Amplitude_T_Sz(i,j)=max(abs(Amplitude(3,:)).^2);
    end
end


[tx,ty]=meshgrid(P0_down/P0,P0_up/P0);
figure(2)
subplot(2,2,1)
pcolor(tx,ty,Amplitude_T_Sx);
colormap(hot);
colorbar
caxis([0,0.2])
shading interp
xlabel('P_{down}/P0')
ylabel('P_{up}/P0')
title('S_x')
subplot(2,2,2)
pcolor(tx,ty,Amplitude_T_Sy);
colormap(hot);
colorbar
caxis([0,0.2])
shading interp
xlabel('P_{down}/P0')
ylabel('P_{up}/P0')
title('S_y')
subplot(2,2,3)
pcolor(tx,ty,Amplitude_T_Sz);
colormap(hot);
colorbar
caxis([0,0.2])
shading interp
xlabel('P_{down}/P0')
ylabel('P_{up}/P0')
title('S_z')



profile viewer