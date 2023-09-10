clear
profile on
V0=0;
gamma_P=1/4;gamma_nR=2*gamma_P;R=0.05;P0=gamma_nR*gamma_P/R;
P0_up=1.1*P0;P0_down=1.1*P0;
scalet=1;scale=1.52*scalet;
dT=0.001;g_nR=0.01;Omega_F=4;g_P=0.004;
nx=30;
F=zeros(2*nx,1);
J_matrix=2*ones(1,nx);
phi=0.2*pi*ones(1,nx);
Driven_T=8.2;
Step=10000;Cord=10;

P0_up=P0*linspace(0.6,1.6,20);
P0_down=P0*linspace(0.6,1.6,20);
Fidelity=zeros(20,20);
for i=1:20
    for j=1:20
    [Vpo,Pump,Coupling,Psi0] = Calculate_potential(nx,phi,V0,P0_up(i),P0_down(j),0.1*pi,0);
    [Result_matrix_up,Result_matrix_down,~,Time] = Iterative_program_Multiple_gpu(dT,Step,Cord,...
        scale*Vpo,scale*J_matrix,scale*Coupling,scalet*Pump,scalet*Driven_T,scale*F,scale*Omega_F,scale*g_P,scale*g_nR,scalet*R,scalet*gamma_P,scalet*gamma_nR,Psi0);
    Result_matrix=cat(1,Result_matrix_up,Result_matrix_down);
    rho=sum(abs(Result_matrix).^2);
    F_time=Psi0'*Result_matrix./sqrt(rho);
    Fidelity(i,j)=mean(abs(F_time).^2);
    end
end

figure(1)
plot(Time,abs(F_time).^2);

[tx,ty]=meshgrid(P0_down/P0,P0_up/P0);
figure(2)
pcolor(tx,ty,Fidelity);
colormap(hot);
colorbar
shading interp
xlabel('P_{down}')
ylabel('P_{up}')



profile viewer