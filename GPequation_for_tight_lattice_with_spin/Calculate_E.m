clear
tic
V0=0;
gamma_P=1/4;gamma_nR=2*gamma_P;R=0.05;P0=gamma_P;
P0_up=2*P0;P0_down=0*P0;P1=0*P0;P2=0*P0;
scalet=1;scale=1.52*scalet;
dT=0.005;g_P=0.005;g_nR=8*g_P;
kp=3;nx=10*kp;
F=zeros(2*nx,1);
J_matrix=-2*ones(1,nx-1);
phi=2*pi*(1:nx)/kp;x0=nx/2;

[Vpo,Pump,Pump_pulse,Coupling,~] = Calculate_potential(nx,phi,V0,P0_up,P0_down,P1,P2,0.5*pi,0*pi,x0);
K_matrix=diag(J_matrix,1)+diag(J_matrix,-1);
A1=zeros(2);A1(2,1)=1;A2=zeros(2);A2(1,2)=1;
H0=kron(eye(2),K_matrix)+kron(A2,K_matrix)+kron(A1,K_matrix')+1j*diag(Pump-gamma_P,0);
[VectorL,Evalue,VectorR]=eig(H0);
Evalue=diag(Evalue);
B=real(Evalue);
[Y,I] = sort(B);
Evalue1=Evalue(I);

site=1:2*nx;
figure(1)
subplot(2,1,1)
plot(site,real(Evalue1),'.r')
subplot(2,1,2)
plot(site,imag(Evalue1),'.r')
figure(2)
plot(site,abs(VectorR(:,I([15,46]))).^2,'-')