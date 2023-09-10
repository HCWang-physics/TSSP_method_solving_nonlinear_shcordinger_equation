clear
profile on
V0=0;
gamma_P=1/4;gamma_nR=2*gamma_P;R=0.05;P0=gamma_nR*gamma_P/R;
P0_up=2*P0;P0_down=2*P0;rate=1;
scalet=1;scale=1.52*scalet;
dT=0.001;g_nR=0.01;Omega_F=4;g_P=0.005;
nx=30;
F=zeros(2*nx,1);
J_matrix=2*ones(1,nx);
phi=0.3*pi*ones(1,nx);
Driven_T=4.21;
Step=5000;Cord=80;Omega_length=400;site=nx/2;

Theta0=linspace(0,0.5*pi,5);
phi0=linspace(-pi,pi,10);
x_site_L=zeros(Step,5*10);
y_site_L=zeros(Step,5*10);
x_site_R=zeros(Step,5*10);
y_site_R=zeros(Step,5*10);
for i=1:5
    for j=1:10
    index=i+10*(j-1);
    [Vpo,Pump,Coupling,Psi0] = Calculate_potential(nx,phi,V0,P0_up,P0_down,Theta0(i),phi0(j));
    [Result_matrix_up,Result_matrix_down,~,~] = Iterative_program_Multiple_gpu(dT,Step,Cord,...
        scale*Vpo,scale*J_matrix,scale*Coupling,scalet*Pump,scalet*Driven_T,scale*F,scale*Omega_F, ...
        scale*g_P,scale*g_nR,scalet*R,scalet*gamma_P,scalet*gamma_nR,rate,Psi0);
    
    
    Sz=(abs(Result_matrix_up).^2-abs(Result_matrix_down).^2)./(abs(Result_matrix_up).^2+abs(Result_matrix_down).^2);
    Psi_H=(Result_matrix_up+Result_matrix_down)/sqrt(2);
    Psi_V=(Result_matrix_up-Result_matrix_down)/sqrt(2);
    Sx=(abs(Psi_H).^2-abs(Psi_V).^2)./(abs(Psi_H).^2+abs(Psi_V).^2);
    Psi_D=(exp(1j*0.25*pi)*Result_matrix_up+exp(-1j*0.25*pi)*Result_matrix_down)/sqrt(2);
    Psi_A=(exp(1j*0.25*pi)*Result_matrix_up-exp(-1j*0.25*pi)*Result_matrix_down)/sqrt(2);
    Sy=(abs(Psi_D).^2-abs(Psi_A).^2)./(abs(Psi_D).^2+abs(Psi_A).^2);
    [Sphi,Stheta,~]=cart2sph(Sx,Sy,Sz);
    Stheta=pi/2-Stheta;
    x_site_L(:,index)=Stheta(site,:)/pi;
    y_site_L(:,index)=Sphi(site,:)/pi;

    [Vpo,Pump,Coupling,Psi0] = Calculate_potential(nx,phi,V0,P0_up,P0_down,Theta0(i)+0.5*pi,phi0(j));
    [Result_matrix_up,Result_matrix_down,nR,Time] = Iterative_program_Multiple_gpu(dT,Step,Cord,...
        scale*Vpo,scale*J_matrix,scale*Coupling,scalet*Pump,scalet*Driven_T,scale*F,scale*Omega_F, ...
        scale*g_P,scale*g_nR,scalet*R,scalet*gamma_P,scalet*gamma_nR,rate,Psi0);
    
    
    Sz=(abs(Result_matrix_up).^2-abs(Result_matrix_down).^2)./(abs(Result_matrix_up).^2+abs(Result_matrix_down).^2);
    Psi_H=(Result_matrix_up+Result_matrix_down)/sqrt(2);
    Psi_V=(Result_matrix_up-Result_matrix_down)/sqrt(2);
    Sx=(abs(Psi_H).^2-abs(Psi_V).^2)./(abs(Psi_H).^2+abs(Psi_V).^2);
    Psi_D=(exp(1j*0.25*pi)*Result_matrix_up+exp(-1j*0.25*pi)*Result_matrix_down)/sqrt(2);
    Psi_A=(exp(1j*0.25*pi)*Result_matrix_up-exp(-1j*0.25*pi)*Result_matrix_down)/sqrt(2);
    Sy=(abs(Psi_D).^2-abs(Psi_A).^2)./(abs(Psi_D).^2+abs(Psi_A).^2);
    [Sphi,Stheta,~]=cart2sph(Sx,Sy,Sz);
    Stheta=pi/2-Stheta;
    x_site_R(:,index)=Stheta(site,:)/pi;
    y_site_R(:,index)=Sphi(site,:)/pi;
    end
end


figure(1)
subplot(1,2,1)
plot(x_site_L,y_site_L,'.','MarkerSize',4)
hold on
plot(x_site_R,y_site_R,'.','MarkerSize',4)
xlabel('\theta(\pi)')
ylabel('\phi(\pi)')
xlim([0,1])
ylim([-1,1])
subplot(1,2,2)
plot(x_site_L,y_site_L,'.r','MarkerSize',4)
hold on
plot(x_site_R,y_site_R,'.k','MarkerSize',4)
xlabel('\theta(\pi)')
ylabel('\phi(\pi)')
xlim([0,1])
ylim([-1,1])

figure(2)
sphere
axis equal
re=[0.75 0.75 0.75];
colormap(re)
shading interp
alpha(0.5)
hold on
plot3(sin(x_site_R*pi).*cos(y_site_R*pi),sin(x_site_R*pi).*sin(y_site_R*pi),cos(x_site_R*pi),'o','MarkerSize',2);






profile viewer