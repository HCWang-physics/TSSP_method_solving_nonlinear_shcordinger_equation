function [Result_matrix_up,Result_matrix_down,nR,Time] = Iterative_program_Multiple_gpu(dT,Step,Cord,Vpo,J_matrix,Coupling,Pump,Pump_pulse,g_P,g_nR,R,gamma_P,gamma_nR,kick,Kicked_period,Initial_Psi)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
[X_number,~]=size(Vpo);
nR=zeros(2*X_number,1);
Result_matrix_up=zeros(X_number,Step);
Result_matrix_down=zeros(X_number,Step);
Kinetic_matrix=zeros(X_number,X_number);
Coupling_M=zeros(X_number,X_number);
Time=zeros(Step,1);
for i=1:X_number-1
    Kinetic_matrix(i,i+1)=J_matrix(1,i);
    Kinetic_matrix(i+1,i)=J_matrix(1,i);
    Coupling_M(i,i)=Coupling(1,i);
end
Coupling_M(X_number,X_number)=Coupling(1,X_number);
K_matrix=kron(eye(2,2),Kinetic_matrix);
A=zeros(2,2);
A(1,2)=1;
Coupling_matrix=kron(A,Coupling_M);
T=0;
Mixed_Matrix=(Coupling_matrix+Coupling_matrix')+K_matrix;
Mixed_Matrix=-1j*dT*Mixed_Matrix;

%Main_iterative_program
clear Kinetic_matrix

iteration_site=1;
Middle_Matrix=Initial_Psi;

for k=1:Step
    for ls=1:Cord
    kicked_site=mod(iteration_site,Kicked_period)+1;
    PumpT=dT*Pump+kick(kicked_site,1)*Pump_pulse(:,1)+kick(kicked_site,2)*Pump_pulse(:,2);
    nR=exp(-(gamma_nR+R*(abs(Middle_Matrix).^2))*dT).*nR+PumpT;
    for i=1:X_number  %Second step:calculate the potential, interaction, gain and loss
        V=Vpo(i,1)+g_nR*nR(i,1);
        W=0.5*(-R*nR(i,1)+gamma_P);
        Psi_0=Middle_Matrix(i,1);
        Psi=abs(Psi_0)^2;
        if  W~=0
            Middle_w=dT*W;
            Middle_Matrix(i,1)=exp(-1j*(0.5*V*dT+g_P*Psi*(1-exp(-Middle_w))/(2*W))-0.5*Middle_w)*Psi_0;
        else
            Middle_Matrix(i,1)=exp(-0.5j*(V+g_P*Psi)*dT)*Psi_0;
        end

        V=Vpo(i,1)+g_nR*nR(i+X_number,1);
        W=0.5*(-R*nR(i+X_number,1)+gamma_P);
        Psi_0=Middle_Matrix(i+X_number,1);
        Psi=abs(Psi_0)^2;
        if  W~=0
            Middle_w=dT*W;
            Middle_Matrix(i+X_number,1)=exp(-1j*(0.5*V*dT+g_P*Psi*(1-exp(-Middle_w))/(2*W))-0.5*Middle_w)*Psi_0;
        else
            Middle_Matrix(i+X_number,1)=exp(-0.5j*(V+g_P*Psi)*dT)*Psi_0;
        end
    end

    Middle_Matrix=expm(Mixed_Matrix)*Middle_Matrix;

    for i=1:X_number  %Second step:calculate the potential, interaction, gain and loss
        V=Vpo(i,1)+g_nR*nR(i,1);
        W=0.5*(-R*nR(i,1)+gamma_P);
        Psi_0=Middle_Matrix(i,1);
        Psi=abs(Psi_0)^2;
        if  W~=0
            Middle_w=dT*W;
            Middle_Matrix(i,1)=exp(-1j*(0.5*V*dT+g_P*Psi*(1-exp(-Middle_w))/(2*W))-0.5*Middle_w)*Psi_0;
        else
            Middle_Matrix(i,1)=exp(-0.5j*(V+g_P*Psi)*dT)*Psi_0;
        end

        V=Vpo(i,1)+g_nR*nR(i+X_number,1);
        W=0.5*(-R*nR(i+X_number,1)+gamma_P);
        Psi_0=Middle_Matrix(i+X_number,1);
        Psi=abs(Psi_0)^2;
        if  W~=0
            Middle_w=dT*W;
            Middle_Matrix(i+X_number,1)=exp(-1j*(0.5*V*dT+g_P*Psi*(1-exp(-Middle_w))/(2*W))-0.5*Middle_w)*Psi_0;
        else
            Middle_Matrix(i+X_number,1)=exp(-0.5j*(V+g_P*Psi)*dT)*Psi_0;
        end
    end
    iteration_site=iteration_site+1;
    T=T+dT;
    end
    Time(k,1)=T;
    Result_matrix_up(:,k)=Middle_Matrix(1:X_number,1);
    Result_matrix_down(:,k)=Middle_Matrix(1+X_number:2*X_number,1);
end
end

