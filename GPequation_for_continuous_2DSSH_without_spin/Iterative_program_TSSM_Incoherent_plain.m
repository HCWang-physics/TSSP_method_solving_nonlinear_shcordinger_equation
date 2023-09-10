function [Result_matrix,nR,Number,Number_boundary,m,T_time] = Iterative_program_TSSM_Incoherent_plain(a0,nx,ny,Boundary,Boundarx,sitea,dT,Cord,Step,Space_X,Space_Y,Vpo,Em,Pump,PumpVortex,Vortex0,Pulse_P,F0,kx,Omega_F,g_P,g_nR,R,gamma_P,gamma_nR,tau)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
[Y_number,X_number]=size(Space_X);
Middle_X=floor(X_number/2);
Middle_Y=floor(Y_number/2);
dx=Space_X(1,2)-Space_X(1,1);
dy=Space_Y(2,1)-Space_Y(1,1);
Number=zeros(1,Step);
T_time=zeros(1,Step);
m=zeros(4,Step);
Number_boundary=zeros(2,Step);
g=gpuDevice(1);
reset(g);
% Y_max=155;X_max=155;
% Y_0=45;X_0=45;
Y_max=floor(0.75*(Y_number-1)/5.5)+1;X_max=floor(0.75*(X_number-1)/5.5)+1;
Y_0=floor(0.25*(Y_number-1)/5.5);X_0=floor(0.25*(X_number-1)/5.5);

Fourier_matrix_X=zeros(X_number,X_number);
TimeRevolution_X=zeros(Y_number,X_number);
Fourier_matrix_X_R=zeros(X_number,X_number);
Fourier_matrix_Y=zeros(Y_number,Y_number);
TimeRevolution_Y=zeros(Y_number,X_number);
Fourier_matrix_Y_R=zeros(Y_number,Y_number);
nR=zeros(Y_number,X_number);
% Noise=zeros(Y_number,X_number);
Lx=dx*X_number;
Ly=dy*Y_number;
% A=X_number*Y_number;
T=0;
% excitedsite=1;
PumpVortexconj=conj(PumpVortex);

Random1=randn(Y_number,X_number)+1j*randn(Y_number,X_number);
% Random1=exp(2j*Space_X-(Space_X-Lx/2).^2-(Space_Y-Ly/2).^2);
Random1(:,1)=0;
Random1(:,X_number)=0;
Random1(1,:)=0;
Random1(Y_number,:)=0;
Random1=Random1-sum(sum(Random1))/(Y_number*X_number);
Model1=sum(sum(abs(Random1).^2*dx*dy));
Result_matrix=Random1/sqrt(Model1);
for i=1:X_number
    mu_x=2*(i-Middle_X-1)*pi/Lx;
	Fourier_matrix_X(1:X_number,i)=exp(-1i*mu_x*Space_X(1,1:X_number))/X_number;
	Fourier_matrix_X_R(i,1:X_number)=exp(1i*mu_x*Space_X(1,1:X_number));
end
for i=1:Y_number
    mu_y=2*(i-Middle_Y-1)*pi/Ly;
	Fourier_matrix_Y(i,1:Y_number)=exp(-1i*mu_y*Space_Y(1:Y_number,1))/Y_number;
	Fourier_matrix_Y_R(1:Y_number,i)=exp(1i*mu_y*Space_Y(1:Y_number,1));
end
for i=1:Y_number
    for j=1:X_number
        mu_x=2*(j-Middle_X-1)*pi/Lx;
        mu_y=2*(i-Middle_Y-1)*pi/Ly;
        TimeRevolution_X(i,j)=exp(-1j*(Em*mu_x^2)*dT/2);
        TimeRevolution_Y(i,j)=exp(-1j*(Em*mu_y^2)*dT/2);
    end
end
Fourier_matrix_X_gpu=gpuArray(Fourier_matrix_X);
TimeRevolution_X_gpu=gpuArray(TimeRevolution_X);
Fourier_matrix_X_R_gpu=gpuArray(Fourier_matrix_X_R);
Fourier_matrix_Y_gpu=gpuArray(Fourier_matrix_Y);
TimeRevolution_Y_gpu=gpuArray(TimeRevolution_Y);
Fourier_matrix_Y_R_gpu=gpuArray(Fourier_matrix_Y_R);

T0=0;
T1=70;
%Main_iterative_program
Middle_Matrix=gpuArray(single(Result_matrix));
nR_gpu=gpuArray(single(nR));
Vpo_gpu=gpuArray(single(Vpo));
Pulse_P_gpu=gpuArray(single(Pulse_P));
Pump_gpu=gpuArray(single(Pump));
PumpVortex_gpu=gpuArray(single(PumpVortex));
PumpVortexconj_gpu=gpuArray(single(PumpVortexconj));
Space_X_gpu=gpuArray(single(Space_X));
Space_Y_gpu=gpuArray(single(Space_Y));
for k=1:Step
    for z=1:Cord
%     Middle_Matrix=Middle_Matrix+(0.5*dT*exp(-((T-T0)/tau)^2))*PumpVortex_gpu.*exp(-1j*Omega_F*T);

    Middle_Matrix=Middle_Matrix*Fourier_matrix_X_gpu.*TimeRevolution_X_gpu*Fourier_matrix_X_R_gpu;%First step:calculate the coupling between up and down and the kinetic term--1
    Middle_Matrix=Fourier_matrix_Y_gpu*Middle_Matrix;%First step:calculate the coupling between up and down and the kinetic term--2
    Middle_Matrix=TimeRevolution_Y_gpu.*Middle_Matrix;
    Middle_Matrix=Fourier_matrix_Y_R_gpu*Middle_Matrix;

    W=0.5*(-R*nR_gpu+gamma_P);
    V=Vpo_gpu+g_nR*nR_gpu;
    Psi=abs(Middle_Matrix).^2;
    middle_w=dT*W;
    Middle_Matrix(W~=0)=exp(-1j*(0.5*V(W~=0)*dT+g_P*Psi(W~=0).*(1-exp(-middle_w(W~=0)))./(2*W(W~=0))) ...
        -0.5*middle_w(W~=0)).*Middle_Matrix(W~=0);
    Middle_Matrix(W==0)=exp(-0.5j*(V(W==0)+g_P*Psi(W==0))*dT).*Middle_Matrix(W==0);

    Middle_Matrix(:,1)=0;Middle_Matrix(:,X_number)=0;Middle_Matrix(1,:)=0;Middle_Matrix(Y_number,:)=0;

%     Middle_Matrix(:,1)=0;Middle_Matrix(:,X_number)=0;Middle_Matrix(1,:)=0;Middle_Matrix(Y_number,:)=0;
%     if T<T1-tau
%         Pulse=Pulse_P_gpu*exp(-((T-T0)/tau)^2)+Pump_gpu;
%         C_vortex=Vortex0*sum(sum(PumpVortexconj_gpu.*Middle_Matrix))*dx*dy;
% %         Middle_Matrix=Middle_Matrix+C_vortex*PumpVortex_gpu*dT;
%         Middle_Matrix=Middle_Matrix+C_vortex*Pulse.*PumpVortex_gpu*dT;
%     else
%         Pulse=0.61*Pulse_P_gpu*exp(-((T-T1)/tau)^2)+Pump_gpu;
%         C_vortex=Vortex0*sum(sum(PumpVortex_gpu.*Middle_Matrix))*dx*dy;
% %         Middle_Matrix=Middle_Matrix+C_vortex*PumpVortexconj_gpu*dT;
%         Middle_Matrix=Middle_Matrix+C_vortex*Pulse.*PumpVortexconj_gpu*dT;
%     end

    Pulse=Pulse_P_gpu+Pump_gpu;
%     Pulse=Pulse_P_gpu*exp(-((T-T0)/tau)^2)+Pump_gpu;
%     C_vortex=(Vortex0*(P1*exp(-((T-T0)/tau)^2)+P0))*sum(sum(PumpVortexconj_gpu.*Middle_Matrix))*dx*dy*dT;
%     Middle_Matrix=Middle_Matrix+C_vortex*PumpVortex_gpu;

%     Middle_Matrix1=Middle_Matrix1+F0*dT*exp(-1j*Omega_F*T-((T-T0)/tau)^2)*(PumpVortex(:,:,1)+PumpVortex(:,:,2));
%     Middle_Matrix1=Middle_Matrix1+3*F0*dT*exp(-1j*Omega_F*T-((T-T1)/tau)^2)*(PumpVortexconj(:,:,1));

%     Middle_Matrix(:,1)=0;Middle_Matrix(:,X_number)=0;Middle_Matrix(1,:)=0;Middle_Matrix(Y_number,:)=0;
    Middle_Matrix(:,1)=0;Middle_Matrix(:,X_number)=0;Middle_Matrix(1,:)=0;Middle_Matrix(Y_number,:)=0;

    W=0.5*(-R*nR_gpu+gamma_P);
    V=Vpo_gpu+g_nR*nR_gpu;
    Psi=abs(Middle_Matrix).^2;
    middle_w=dT*W;
    Middle_Matrix(W~=0)=exp(-1j*(0.5*V(W~=0)*dT+g_P*Psi(W~=0).*(1-exp(-middle_w(W~=0)))./(2*W(W~=0))) ...
        -0.5*middle_w(W~=0)).*Middle_Matrix(W~=0);
    Middle_Matrix(W==0)=exp(-0.5j*(V(W==0)+g_P*Psi(W==0))*dT).*Middle_Matrix(W==0);
    
    Middle_Matrix=Fourier_matrix_Y_gpu*Middle_Matrix;
    Middle_Matrix=TimeRevolution_Y_gpu.*Middle_Matrix;
    Middle_Matrix=Fourier_matrix_Y_R_gpu*Middle_Matrix;
    Middle_Matrix=Middle_Matrix*Fourier_matrix_X_gpu.*TimeRevolution_X_gpu*Fourier_matrix_X_R_gpu;

    nR_gpu=exp(-(gamma_nR+R*(abs(Middle_Matrix).^2))*dT).*nR_gpu+dT*Pulse;
%     nR=exp(-(gamma_nR+R*(abs(Middle_Matrix2).^2))*dT).*nR+dT*Pump;
    T=T+dT;
%     Middle_Matrix=Middle_Matrix+(0.5*dT*exp(-((T-T0)/tau)^2))*PumpVortex_gpu.*exp(-1j*Omega_F*T);
    %Noise=sqrt(dT.*(R*nR+gamma_P)./(dx*dy)).*(randn(Y_number,X_number)+1j*randn(Y_number,X_number));
    %Noise(:,1)=0;
    %Noise(:,X_number)=0;
    %Noise(1,:)=0;
    %Noise(Y_number,:)=0;
    %Middle_Matrix2=Middle_Matrix2+Noise;
    end
    Number(1,k)=sum(sum(abs(Middle_Matrix).^2*dx*dy));
    T_time(1,k)=T;

    xsite=0;ysite=0;
    Psi=Middle_Matrix(Y_0:Y_max,X_0:X_max);
    [FX,FY] = gradient(Psi);
    PartialPsi=((Space_X_gpu(Y_0:Y_max,X_0:X_max)-a0/2+sitea-xsite*a0-Boundarx).*FY./dy-(Space_Y_gpu(Y_0:Y_max,X_0:X_max)-a0/2+sitea-ysite*a0-Boundary).*FX./dx);
	Lz=-1j*sum(sum(conj(Psi).*PartialPsi))*dx*dy;
    Number1=sum(sum(abs(Psi).^2*dx*dy));
    m(1,k)=Lz/Number1;
    Number_boundary(1,k)=Number1;

    xsite=nx-1;ysite=0;
    Psi=Middle_Matrix(Y_0:Y_max,X_number-X_max:X_number-X_0);
    [FX,FY] = gradient(Psi);
    PartialPsi=((Space_X_gpu(Y_0:Y_max,X_number-X_max:X_number-X_0)-a0/2-sitea-xsite*a0-Boundarx).*FY./dy-(Space_Y_gpu(Y_0:Y_max,X_number-X_max:X_number-X_0)-a0/2+sitea-ysite*a0-Boundary).*FX./dx);
	Lz=-1j*sum(sum(conj(Psi).*PartialPsi))*dx*dy;
    Number2=sum(sum(abs(Psi).^2*dx*dy));
    m(2,k)=Lz/Number2;
    Number_boundary(2,k)=Number2;

    xsite=0;ysite=ny-1;
    Psi=Middle_Matrix(Y_number-Y_max:Y_number-Y_0,X_0:X_max);
    [FX,FY] = gradient(Psi);
    PartialPsi=((Space_X_gpu(Y_number-Y_max:Y_number-Y_0,X_0:X_max)-a0/2+sitea-xsite*a0-Boundarx).*FY./dy-(Space_Y_gpu(Y_number-Y_max:Y_number-Y_0,X_0:X_max)-a0/2-sitea-ysite*a0-Boundary).*FX./dx);
	Lz=-1j*sum(sum(conj(Psi).*PartialPsi))*dx*dy;
    Number3=sum(sum(abs(Psi).^2*dx*dy));
    m(3,k)=Lz/Number3;
    Number_boundary(3,k)=Number3;

    xsite=nx-1;ysite=ny-1;
    Psi=Middle_Matrix(Y_number-Y_max:Y_number-Y_0,X_number-X_max:X_number-X_0);
    [FX,FY] = gradient(Psi);
    PartialPsi=((Space_X_gpu(Y_number-Y_max:Y_number-Y_0,X_number-X_max:X_number-X_0)-a0/2-sitea-xsite*a0-Boundarx).*FY./dy-(Space_Y_gpu(Y_number-Y_max:Y_number-Y_0,X_number-X_max:X_number-X_0)-a0/2-sitea-ysite*a0-Boundary).*FX./dx);
	Lz=-1j*sum(sum(conj(Psi).*PartialPsi))*dx*dy;
    Number4=sum(sum(abs(Psi).^2*dx*dy));
    m(4,k)=Lz/Number4;
    Number_boundary(4,k)=Number4;
%     xsite=0;ysite=0;
%     Psi=Middle_Matrix(25:75,25:75);
%     [FX,FY] = gradient(Psi);
%     PartialPsi=((Space_X_gpu(25:75,25:75)-a0/2+sitea-xsite*a0-Boundarx).*FY./dy-(Space_Y_gpu(25:75,25:75)-a0/2+sitea-ysite*a0-Boundary).*FX./dx);
% 	Lz=-1j*sum(sum(conj(Psi).*PartialPsi))*dx*dy;
%     Number1=sum(sum(abs(Psi).^2*dx*dy));
%     m(1,k)=Lz/Number1;
%     Number_boundary(1,k)=Number1;
% 
%     xsite=nx-1;ysite=0;
%     Psi=Middle_Matrix(25:75,X_number-75:X_number-25);
%     [FX,FY] = gradient(Psi);
%     PartialPsi=((Space_X_gpu(25:75,X_number-75:X_number-25)-a0/2-sitea-xsite*a0-Boundarx).*FY./dy-(Space_Y_gpu(25:75,X_number-75:X_number-25)-a0/2+sitea-ysite*a0-Boundary).*FX./dx);
% 	Lz=-1j*sum(sum(conj(Psi).*PartialPsi))*dx*dy;
%     Number2=sum(sum(abs(Psi).^2*dx*dy));
%     m(2,k)=Lz/Number2;
%     Number_boundary(2,k)=Number2;
% 
%     xsite=0;ysite=ny-1;
%     Psi=Middle_Matrix(Y_number-75:Y_number-25,25:75);
%     [FX,FY] = gradient(Psi);
%     PartialPsi=((Space_X_gpu(Y_number-75:Y_number-25,25:75)-a0/2+sitea-xsite*a0-Boundarx).*FY./dy-(Space_Y_gpu(Y_number-75:Y_number-25,25:75)-a0/2-sitea-ysite*a0-Boundary).*FX./dx);
% 	Lz=-1j*sum(sum(conj(Psi).*PartialPsi))*dx*dy;
%     Number3=sum(sum(abs(Psi).^2*dx*dy));
%     m(3,k)=Lz/Number3;
%     Number_boundary(3,k)=Number3;
% 
%     xsite=nx-1;ysite=ny-1;
%     Psi=Middle_Matrix(Y_number-75:Y_number-25,X_number-75:X_number-25);
%     [FX,FY] = gradient(Psi);
%     PartialPsi=((Space_X_gpu(Y_number-75:Y_number-25,X_number-75:X_number-25)-a0/2-sitea-xsite*a0-Boundarx).*FY./dy-(Space_Y_gpu(Y_number-75:Y_number-25,X_number-75:X_number-25)-a0/2-sitea-ysite*a0-Boundary).*FX./dx);
% 	Lz=-1j*sum(sum(conj(Psi).*PartialPsi))*dx*dy;
%     Number4=sum(sum(abs(Psi).^2*dx*dy));
%     m(4,k)=Lz/Number4;
%     Number_boundary(4,k)=Number4;
end
Result_matrix=gather(Middle_Matrix);
end

