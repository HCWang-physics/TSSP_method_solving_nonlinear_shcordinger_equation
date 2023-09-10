clear
profile on
V0=0;
gamma_P=1/6;gamma_nR=2*gamma_P;R=0.01;P0=gamma_nR*gamma_P/R;
scalet=1;scale=1.52*scalet;rate=1;
dT=0.001;g_nR=0.01;Omega_F=4;g_P=0.005;
nx=20;
F=zeros(2*nx,1);
J_matrix=2*ones(1,nx);
phi=0.3*pi*ones(1,nx);
Number=20;
Driven_Time=linspace(2,8,Number);
Chaos_K=zeros(Number,1);
Step=5000;Cord=80;Omega_length=400;site=nx/2;

P0_up=P0*linspace(0.6,1.6,20);
P0_down=P0*linspace(0.6,1.6,20);
Amplitude_T_Sx=zeros(20,Omega_length);
Amplitude_T_Sy=zeros(20,Omega_length);
Amplitude_T_Sz=zeros(20,Omega_length);
for i=1:Number
    Driven_T=Driven_Time(i);
    [Vpo,Pump,Coupling,Psi0] = Calculate_potential(nx,phi,V0,2*P0,2*P0,0.5*pi,0);
    [Result_matrix_up,Result_matrix_down,~,Time] = Iterative_program_Multiple_gpu(dT,Step,Cord,...
        scale*Vpo,scale*J_matrix,scale*Coupling,scalet*Pump,scalet*Driven_T,scale*F,scale*Omega_F, ...
        scale*g_P,scale*g_nR,scalet*R,scalet*gamma_P,scalet*gamma_nR,rate,Psi0);
    Sz=(abs(Result_matrix_up).^2-abs(Result_matrix_down).^2)./(abs(Result_matrix_up).^2+abs(Result_matrix_down).^2);
    Psi_H=(Result_matrix_up+Result_matrix_down)/sqrt(2);
    Psi_V=(Result_matrix_up-Result_matrix_down)/sqrt(2);
    Sx=(abs(Psi_H).^2-abs(Psi_V).^2)./(abs(Psi_H).^2+abs(Psi_V).^2);
    Psi_D=(exp(1j*0.25*pi)*Result_matrix_up+exp(-1j*0.25*pi)*Result_matrix_down)/sqrt(2);
    Psi_A=(exp(1j*0.25*pi)*Result_matrix_up-exp(-1j*0.25*pi)*Result_matrix_down)/sqrt(2);
    Sy=(abs(Psi_D).^2-abs(Psi_A).^2)./(abs(Psi_D).^2+abs(Psi_A).^2);
    [~,Sphi,Stheta]=cart2sph(Sx,Sy,Sz);
    Omega0=2*pi/Driven_T;
    
    omega=linspace(0.1*Omega0,1*Omega0,Omega_length);
    Amplitude=zeros(4,Omega_length);
    for k=1:Omega_length
        Exptime=exp(-1j*Time*omega(k));
        Amplitude(1,k)=sum((Sx(site,:))*Exptime)/Step;
        Amplitude(2,k)=sum((Sy(site,:))*Exptime)/Step;
        Amplitude(3,k)=sum((Sz(site,:))*Exptime)/Step;
    end
    Amplitude_T_Sx(i,:)=abs(Amplitude(1,:)).^2;
    Amplitude_T_Sy(i,:)=abs(Amplitude(2,:)).^2;
    Amplitude_T_Sz(i,:)=abs(Amplitude(3,:)).^2;

    Phi=Sphi(site,:)+Stheta(site,:);
    c_c=linspace(pi/5,4*pi/5,20);al1=0.1;
    Kc_c=zeros(20,1);
    Ephi=mean(Phi);
    for j=1:20
    c=c_c(j);
    qc_t=zeros(Step,1);
    pc_t=zeros(Step,1);
    thetac_t=zeros(Step,1);
    qc=0;
    pc=0;
    thetac=0;
    for k=1:Step
        pc=pc+Phi(k)*cos(thetac);
        qc=qc+Phi(k)*sin(thetac);
        thetac=thetac+c+al1*Phi(k);
        qc_t(k,1)=qc;
        pc_t(k,1)=pc;
        thetac_t(k,1)=thetac;
    end
    N_Cut=Step/10;
    IntervalN=1:N_Cut;
    Mc=zeros(N_Cut,1);
    Dc=zeros(N_Cut,1);
    for k=1:N_Cut
        N_x=IntervalN(k);
        Mc(k,1)=(sum((pc_t(N_x+1:Step-N_x,1)-pc_t(1:Step-2*N_x,1)).^2)+sum((qc_t(N_x+1:Step-N_x,1)-qc_t(1:Step-2*N_x,1)).^2))/(Step-2*N_x);
        Dc(k,1)=Mc(k,1)-Ephi^2*(1-cos(N_x*c))/(1-cos(c));
    end
    Kc=(IntervalN-mean(IntervalN))*(Dc-mean(Dc))/N_Cut;
    Kc=Kc/sqrt(mean((IntervalN-mean(IntervalN)).^2));
    Kc=Kc/sqrt(mean((Dc-mean(Dc)).^2));
    Kc_c(j,1)=Kc;
    end
    Chaos_K(i,1)=mean(Kc_c);
end


[tx,ty]=meshgrid(omega/Omega0,Driven_Time);
figure(2)
subplot(2,2,1)
pcolor(tx,ty,Amplitude_T_Sx);
colormap(hot);
colorbar
shading interp
xlabel('\omega')
ylabel('DrivenT(ps)')
title('S_x')
subplot(2,2,2)
pcolor(tx,ty,Amplitude_T_Sy);
colormap(hot);
colorbar
shading interp
xlabel('\omega')
ylabel('DrivenT(ps)')
title('S_y')
subplot(2,2,3)
pcolor(tx,ty,Amplitude_T_Sz);
colormap(hot);
colorbar
shading interp
xlabel('\omega')
ylabel('DrivenT(ps)')
title('S_z')
subplot(2,2,4)
plot(Driven_Time,Chaos_K,'-or')
xlabel('DrivenT(ps)')
ylabel('K')



profile viewer