%%114514
tic
a0=4;
sitea=1.2*a0/4;
redus=0.8*a0/4; 
V0=40;
ny=5;
nx=5;
Boundary=0.25*a0;
Boundarx=0.25*a0;
Ly=ny*a0+2*Boundary;
Lx=nx*a0+2*Boundarx;

Xnumber=1100+1;Ynumber=1100+1;%discrete point number along x and y direction
dx=Lx/Xnumber
dy=Ly/Ynumber
F=zeros(Ynumber,Xnumber);
F0=0;
gamma_P=1/8;R=0.05;gamma_nR=2*gamma_P;P0=gamma_P*gamma_nR/R;
[SpaceX,SpaceY,Vpo,Pump,Pulse_P,PumpVortex1,PumpVortex2] = Calculate_potential_all(a0,nx,ny,Boundary,Boundarx,redus,sitea, ...
    Xnumber,Ynumber,V0,1.5*P0,3*P0);

%% 1919810
scalet=1;
scale=1.52;
Step=3000;
Cord=50;
dT=0.002;
m=25;
Em=76.3/(2*m);
Omega_F=18.7;
kx=0*pi/(5*a0);
g_P=0.008;
g_nR=0;
tau=40;
[Result_matrix4,nR4,Number4,Number_boundary4,m4,T_time4] = Iterative_program_TSSM_Incoherent(a0,nx,ny,Boundary,Boundarx,sitea,dT,Cord,Step,SpaceX,SpaceY,...
    scale*Vpo,scale*Em,scalet*Pump,PumpVortex1,PumpVortex2,scale*0.005,scalet*Pulse_P,scale*F0,kx,scale*Omega_F,scale*g_P, ...
    scale*g_nR,scalet*R,scalet*gamma_P,scalet*gamma_nR,scalet*tau);
% save 'm_LU-1_RB1_0p235.mat' m
%%


theta=linspace(0,2*pi,100);
x=zeros(100,4*nx*ny);
y=zeros(100,4*nx*ny);
for i=1:nx
    for j=1:ny
        site=i+(j-1)*nx;
        x(:,4*(site-1)+1)=a0/2+sitea+(i-1)*a0+Boundarx+redus*cos(theta);
        y(:,4*(site-1)+1)=a0/2+sitea+(j-1)*a0+Boundary+redus*sin(theta);
        x(:,4*(site-1)+2)=a0/2-sitea+(i-1)*a0+Boundarx+redus*cos(theta);
        y(:,4*(site-1)+2)=a0/2-sitea+(j-1)*a0+Boundary+redus*sin(theta);
        x(:,4*(site-1)+3)=a0/2+sitea+(i-1)*a0+Boundarx+redus*cos(theta);
        y(:,4*(site-1)+3)=a0/2-sitea+(j-1)*a0+Boundary+redus*sin(theta);
        x(:,4*(site-1)+4)=a0/2-sitea+(i-1)*a0+Boundarx+redus*cos(theta);
        y(:,4*(site-1)+4)=a0/2+sitea+(j-1)*a0+Boundary+redus*sin(theta);
    end
end
%%
figure(1);
subplot(2,2,1);
pcolor(SpaceX,SpaceY,abs(Result_matrix4(:,:)).^2);
hold on
plot(x,y,'--black');
colormap(hot);
colorbar
shading interp
xlabel('x(\mum)')
ylabel('y(\mum)')
subplot(2,2,2);
pcolor(SpaceX,SpaceY,angle(Result_matrix4(:,:)));
hold on
plot(x,y,'--black');
colorbar
colormap(hot);
shading interp
xlabel('x(\mum)')
ylabel('y(\mum)')
xlim([0,Lx])
subplot(2,2,3)
plot(T_time4,real(m4(1,:)),'-b','LineWidth',1.5);
hold on
plot(T_time4,real(m4(2,:)),'-r','LineWidth',1.5);
hold on
plot(T_time4,real(m4(3,:)),'-k','LineWidth',1.5);
hold on
plot(T_time4,real(m4(4,:)),'-m','LineWidth',1.5);
legend('LB','RB','LU','RU')
ylabel('m')
xlim([0,300])
xlabel('T(ps)')
ylim([-3,3])
subplot(2,2,4)
plot(T_time4,Number_boundary4(1,:),'-b','LineWidth',2);
hold on
plot(T_time4,Number_boundary4(2,:),'-r','LineWidth',2);
hold on
plot(T_time4,Number_boundary4(3,:),'-k','LineWidth',2);
hold on
plot(T_time4,Number_boundary4(4,:),'-m','LineWidth',2);
xlabel('T(ps)')
ylabel('Density intensity')
xlim([0,300])
ylim([0,2])
legend('LB','RB','LU','RU')
%%
figure(2);
subplot(2,2,1);
pcolor(SpaceX,SpaceY,angle(PumpVortex(:,:)));
clim([-pi pi])
colormap(hot);
colorbar
shading interp
xlabel('x(\mum)')
ylabel('y(\mum)')
subplot(2,2,2);
pcolor(SpaceX,SpaceY,Pulse_P);
colorbar
colormap(hot);
shading interp
xlabel('x(\mu)')
ylabel('y(\mu)')
xlim([0,Lx])
subplot(2,2,3);
ylim([0,Ly])
pcolor(SpaceX,SpaceY,Vpo);
colorbar
colormap(hot);
shading interp
xlabel('x(\mu)')
ylabel('y(\mu)')
xlim([0,Lx])
subplot(2,2,4);
ylim([0,Ly])
pcolor(SpaceX,SpaceY,abs(PumpVortex(:,:)));
colorbar
colormap(hot);
shading interp
xlabel('x(\mu)')
ylabel('y(\mu)')
xlim([0,Lx])
%%
figure(3);
plot(T_time4,Number4(1,:),'-b','LineWidth',1.5);
%%
toc