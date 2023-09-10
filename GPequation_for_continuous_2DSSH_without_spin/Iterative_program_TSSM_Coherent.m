function [Result_matrix] = Iterative_program_TSSM_Coherent(dT,Cord,Step,Space_X,Space_Y,Vpo,Em,Pump,F,Omega_F,kx,g_P,alpha,g_Pump,R,gamma_P,tau)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
rotation=0;
[Y_number,X_number]=size(Space_X);
dx=Space_X(1,2)-Space_X(1,1);
dy=Space_Y(2,1)-Space_Y(1,1);
Middle_X=fix(X_number/2);
Middle_Y=fix(Y_number/2);
Result_matrix=zeros(Y_number,X_number,Step+1);
Fourier_matrix_X=zeros(X_number,X_number);
TimeRevolution_X=zeros(Y_number,X_number);
Fourier_matrix_X_R=zeros(X_number,X_number);
Fourier_matrix_Y=zeros(Y_number,Y_number);
TimeRevolution_Y=zeros(Y_number,X_number);
Fourier_matrix_Y_R=zeros(Y_number,Y_number);
Noise=zeros(Y_number,X_number);
Lx=dx*X_number;
Ly=dy*Y_number;
A=X_number*Y_number;
T=0;

Random1=randn(Y_number,X_number)+1j*randn(Y_number,X_number);
Random1(:,1)=0;
Random1(:,X_number)=0;
Random1(1,:)=0;
Random1(Y_number,:)=0;
Model1=sum(sum(abs(Random1).^2*dx*dy));
Result_matrix(:,:,1)=Random1/sqrt(Model1);
for i=1:X_number
    mu_x=(i)*pi/Lx;
	Fourier_matrix_X(1:X_number,i)=2*sin(mu_x*Space_X(1,1:X_number))/X_number;
	Fourier_matrix_X_R(i,1:X_number)=sin(mu_x*Space_X(1,1:X_number));
end
for i=1:Y_number
    mu_y=(i)*pi/Ly;
	Fourier_matrix_Y(i,1:Y_number)=2*sin(mu_y*Space_Y(1:Y_number,1))/Y_number;
	Fourier_matrix_Y_R(1:Y_number,i)=sin(mu_y*Space_Y(1:Y_number,1));
end
for i=1:Y_number
    for j=1:X_number
        mu_x=(j)*pi/Lx;
        mu_y=(i)*pi/Ly;
        TimeRevolution_X(i,j)=exp(-1j*(2*Em*mu_x^2+2*rotation*Space_Y(i,j))*dT/4);
        TimeRevolution_Y(i,j)=exp(-1j*(2*Em*mu_y^2+2*rotation*Space_X(i,j))*dT/4);
    end
end
Middle_Matrix=zeros(Y_number,X_number);
Middle_Matrix1=zeros(Y_number,X_number);
Middle_Matrix2=zeros(Y_number,X_number);

%Main_iterative_program
Middle_Matrix2=Result_matrix(:,:,1);
for k=1:Step
    for z=1:Cord
    Middle_Matrix=Middle_Matrix2*Fourier_matrix_X.*TimeRevolution_X*Fourier_matrix_X_R;%First step:calculate the the kinetic term--1
    Middle_Matrix=Fourier_matrix_Y*Middle_Matrix;
    Middle_Matrix=TimeRevolution_Y.*Middle_Matrix;
    Middle_Matrix=Fourier_matrix_Y_R*Middle_Matrix;
    
    for i=1:Y_number  %Second step:calculate the potential, interaction, gain and loss--2
        for j=1:X_number
            V=Vpo(i,j)+g_Pump*Pump(i,j);
            Psi=abs(Middle_Matrix(i,j))^2;
            W=-R*Pump(i,j)/2+gamma_P/2+alpha*Psi;
            if  W==0;
                Middle_Matrix1(i,j)=exp(-1j*(V+g_P*Psi)*dT)*Middle_Matrix(i,j);
            else
                Middle_Matrix1(i,j)=exp(-1j*(V*dT+g_P*Psi*(1-exp(-2*dT*W))/(2*W))-dT*W)*Middle_Matrix(i,j);
            end
        end
    end
    %Middle_Matrix1=Middle_Matrix1+dT.*cos(kx*Space_X).*F(:,:).*exp(-1j*Omega_F*T);%-(T/tau)^2
    if T<200
        Middle_Matrix1=Middle_Matrix1+(5*(T-0)/200).*dT.*cos(kx*Space_X).*F(:,:).*exp(-1j*Omega_F*T);
    else
        Middle_Matrix1=Middle_Matrix1+(5*(400-T)/200).*dT.*cos(kx*Space_X).*F(:,:).*exp(-1j*Omega_F*T);
    end
    
    Middle_Matrix2=Fourier_matrix_Y*Middle_Matrix1;%Third step:calculate the kinetic term--3
    Middle_Matrix2=TimeRevolution_Y.*Middle_Matrix2;
    Middle_Matrix2=Fourier_matrix_Y_R*Middle_Matrix2;
    Middle_Matrix2=Middle_Matrix2*Fourier_matrix_X.*TimeRevolution_X*Fourier_matrix_X_R;
	T=T+dT;
    end
    Result_matrix(:,:,k+1)=Middle_Matrix2;
end
end

