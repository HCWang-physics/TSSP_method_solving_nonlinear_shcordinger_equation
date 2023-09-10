clear
profile on
V0=0;
gamma_P=1/6;gamma_nR=2*gamma_P;R=0.01;P0=gamma_nR*gamma_P/R;
P0_up=1.2*P0;P0_down=1.0*P0;
scalet=1;scale=1.52*scalet;
dT=0.001;g_nR=0.01;Omega_F=4;g_P=0.005;
nx=30;
rng(1);
F=zeros(2*nx,1);
J_matrix=2*ones(1,nx);
phi=0.3*pi*ones(1,nx);
[Vpo,Pump,Coupling,Psi0] = Calculate_potential(nx,phi,V0,P0_up,P0_down);
site=nx/2;Step=20000;Cord=10;Omega_length=400;

Driven_T0=linspace(2,10,40);
Amplitude_T_Sx=zeros(400,40);
Amplitude_T_Sy=zeros(400,40);
Amplitude_T_Sz=zeros(400,40);
for i=1:40
    Driven_T=Driven_T0(i);
    [Result_matrix_up,Result_matrix_down,nR,Time] = Iterative_program_Multiple_gpu(dT,Step,Cord,...
        scale*Vpo,scale*J_matrix,scale*Coupling,scalet*Pump,scalet*Driven_T,scale*F,scale*Omega_F,scale*g_P,scale*g_nR,scalet*R,scalet*gamma_P,scalet*gamma_nR,Psi0);
    
    Sz=(abs(Result_matrix_up).^2-abs(Result_matrix_down).^2)./(abs(Result_matrix_up).^2+abs(Result_matrix_down).^2);
    Psi_H=(Result_matrix_up+Result_matrix_down)/sqrt(2);
    Psi_V=(Result_matrix_up-Result_matrix_down)/sqrt(2);
    Sx=(abs(Psi_H).^2-abs(Psi_V).^2)./(abs(Psi_H).^2+abs(Psi_V).^2);
    Psi_D=(exp(1j*0.25*pi)*Result_matrix_up+exp(-1j*0.25*pi)*Result_matrix_down)/sqrt(2);
    Psi_A=(exp(1j*0.25*pi)*Result_matrix_up-exp(-1j*0.25*pi)*Result_matrix_down)/sqrt(2);
    Sy=(abs(Psi_D).^2-abs(Psi_A).^2)./(abs(Psi_D).^2+abs(Psi_A).^2);
    %[~,Sphi,Stheta]=cart2sph(Sx,Sy,Sz);
    
    
    omega=linspace(0,2*pi/Driven_T,Omega_length);
    Amplitude=zeros(4,Omega_length);
    for j=1:Omega_length
        Exptime=exp(-1j*Time*omega(j));
        Amplitude(1,j)=sum((Sx(site,:))*Exptime)/Step;
        Amplitude(2,j)=sum((Sy(site,:))*Exptime)/Step;
        Amplitude(3,j)=sum((Sz(site,:))*Exptime)/Step;
    end
    Omega0=2*pi/Driven_T;
    Amplitude_T_Sx(:,i)=abs(Amplitude(1,:)).^2;
    Amplitude_T_Sy(:,i)=abs(Amplitude(2,:)).^2;
    Amplitude_T_Sz(:,i)=abs(Amplitude(3,:)).^2;
end

[tx,ty]=meshgrid(Driven_T0,omega/Omega0);
figure(1);
subplot(2,2,1)
pcolor(tx,ty,Amplitude_T_Sx);
colormap(hot);
colorbar
shading interp
xlabel('T0(ps)')
ylabel('\omega/\omega_{0}')
subplot(2,2,2)
pcolor(tx,ty,Amplitude_T_Sy);
colormap(hot);
colorbar
shading interp
xlabel('T0(ps)')
ylabel('\omega/\omega_{0}')
subplot(2,2,3)
pcolor(tx,ty,Amplitude_T_Sz);
colormap(hot);
colorbar
shading interp
xlabel('T0(ps)')
ylabel('\omega/\omega_{0}')



profile viewer