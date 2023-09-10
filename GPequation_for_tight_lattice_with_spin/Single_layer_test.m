clear
tic
V0=0;
gamma_P=1/4;gamma_nR=2*gamma_P;R=0.05;P0=gamma_nR*gamma_P/R;
P0_up=P0;P0_down=0*P0;P1=1*P0;P2=0*P0;
scalet=1;scale=1.52*scalet;
dT=0.005;g_P=0.005;g_nR=8*g_P;
kp=3;nx=10*kp;
F=zeros(2*nx,1);
J_matrix=-2*ones(1,nx-1);
phi=2*pi*(1:nx)/kp;x0=nx/2;
[Vpo,Pump,Pump_pulse,Coupling,Psi0] = Calculate_potential(nx,phi,V0,P0_up,P0_down,P1,P2,0.5*pi,0*pi,x0);

%%
Step=2000;Cord=10;
Kicked_period=2000;Kick_T1=500;Kick_T2=1500;
Driven_T=dT*Kicked_period;
kick=zeros(Kicked_period,2);
kick(Kick_T1,1)=1;
kick(Kick_T2,2)=1;
[Result_matrix_up,Result_matrix_down,nR,Time] = Iterative_program_Multiple_gpu(dT,Step,Cord,...
    scale*Vpo,scale*J_matrix,scale*Coupling,scalet*Pump,scalet*Pump_pulse,scale*g_P, ...
    scale*g_nR,scalet*R,scalet*gamma_P,scalet*gamma_nR,kick,Kicked_period,Psi0);

%%
Sz=(abs(Result_matrix_up).^2-abs(Result_matrix_down).^2)./(abs(Result_matrix_up).^2+abs(Result_matrix_down).^2);
Psi_H=(Result_matrix_up+Result_matrix_down)/sqrt(2);
Psi_V=(Result_matrix_up-Result_matrix_down)/sqrt(2);
Sx=(abs(Psi_H).^2-abs(Psi_V).^2)./(abs(Psi_H).^2+abs(Psi_V).^2);
Psi_D=(exp(1j*0.25*pi)*Result_matrix_up+exp(-1j*0.25*pi)*Result_matrix_down)/sqrt(2);
Psi_A=(exp(1j*0.25*pi)*Result_matrix_up-exp(-1j*0.25*pi)*Result_matrix_down)/sqrt(2);
Sy=(abs(Psi_D).^2-abs(Psi_A).^2)./(abs(Psi_D).^2+abs(Psi_A).^2);
[Sphi,Stheta,~]=cart2sph(Sx,Sy,Sz);
Stheta=pi/2-Stheta;
%%
site=1:nx;
[xt,yt]=meshgrid(Time,site);
figure(1);
subplot(2,2,1)
pcolor(xt,yt,abs(Result_matrix_up).^2);
colormap(hot);
colorbar
shading interp
xlabel('T(T0)')
ylabel('site')
subplot(2,2,2)
pcolor(xt,yt,abs(Result_matrix_down).^2+abs(Result_matrix_up).^2);
colormap(hot);
colorbar
xlabel('T(T0)')
ylabel('site')
shading interp
subplot(2,2,3)
pcolor(xt,yt,Sx);
colormap(hot);
colorbar
xlabel('T(T0)')
ylabel('site')
shading interp
subplot(2,2,4)
pcolor(xt,yt,Sz);
colormap(hot);
colorbar
xlabel('T(T0)')
ylabel('site')
shading interp

figure(2);
subplot(2,2,1)
plot(Time,sum(abs(Result_matrix_up).^2));
shading interp
xlabel('T(T0)')
ylabel('site')
subplot(2,2,2)
plot(Time,sum(abs(Result_matrix_down).^2));
xlabel('T(T0)')
ylabel('site')
shading interp
%%
site=nx/2;
figure(3)
sphere
axis equal
re=[0.75 0.75 0.75];
colormap(re)
shading interp
alpha(0.5)
hold on
plot3(Sx(site,:),Sy(site,:),Sz(site,:),'ok','LineWidth',1,'MarkerSize',2);
hold on
plot3(Sx(1,:),Sy(1,:),Sz(1,:),'or','LineWidth',1,'MarkerSize',2);

Omega0=2*pi/Driven_T;
Omega_length=400;
omega=linspace(0.4*Omega0,0.6*Omega0,Omega_length);
Amplitude=zeros(3,Omega_length);
for i=1:Omega_length
    Exptime=exp(-1j*Time*omega(i));
    Amplitude(1,i)=sum((Sx(site,:))*Exptime)/Step;
    Amplitude(2,i)=sum((Sy(site,:))*Exptime)/Step;
    Amplitude(3,i)=sum((Sz(site,:))*Exptime)/Step;
end
figure(4)
plot(omega/Omega0,abs(Amplitude(1,:)).^2,'-ok','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3);
hold on
plot(omega/Omega0,abs(Amplitude(2,:)).^2,'-hm','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',3);
hold on
plot(omega/Omega0,abs(Amplitude(3,:)).^2,'-db','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',3);
legend('Sx','Sy','Sz')
ylabel('Amplitude')
xlabel('\omega')

Phi=Sphi(site,:)+Stheta(site,:);
al1=1;
Kc_c=zeros(100,1);
Ephi=mean(Phi);
for j=1:100
c=3*pi*rand(1,1)/5+pi/5;
qc_t=zeros(Step,1);
pc_t=zeros(Step,1);
thetac_t=zeros(Step,1);
qc=0;
pc=0;
thetac=0;
for i=1:Step
    pc=pc+Phi(i)*cos(thetac);
    qc=qc+Phi(i)*sin(thetac);
    thetac=thetac+c+al1*Phi(i);
    qc_t(i,1)=qc;
    pc_t(i,1)=pc;
    thetac_t(i,1)=thetac;
end
N_Cut=Step/10;
IntervalN=1:N_Cut;
Mc=zeros(N_Cut,1);
Dc=zeros(N_Cut,1);
for i=1:N_Cut
    N_x=IntervalN(i);
    Mc(i,1)=(sum((pc_t(N_x+1:Step-N_x,1)-pc_t(1:Step-2*N_x,1)).^2)+sum((qc_t(N_x+1:Step-N_x,1)-qc_t(1:Step-2*N_x,1)).^2))/(Step-2*N_x);
    Dc(i,1)=Mc(i,1)-Ephi^2*(1-cos(N_x*c))/(1-cos(c));
end
Kc=(IntervalN-mean(IntervalN))*(Dc-mean(Dc))/N_Cut;
Kc=Kc/sqrt(mean((IntervalN-mean(IntervalN)).^2));
Kc=Kc/sqrt(mean((Dc-mean(Dc)).^2));
Kc_c(j,1)=Kc;
end
median(Kc_c)

figure(5)
subplot(4,1,1)
plot(qc_t,pc_t,'ok',MarkerSize=1);
subplot(4,1,2)
plot(Time,Phi,'-ok',MarkerSize=1);
ylabel('N')
subplot(4,1,3)
plot(IntervalN,Mc,'-or',MarkerSize=1);
hold on
plot(IntervalN,Dc,'--ok',MarkerSize=1);


toc