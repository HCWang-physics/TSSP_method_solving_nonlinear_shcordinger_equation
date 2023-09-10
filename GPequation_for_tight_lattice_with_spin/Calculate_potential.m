function [Ver,Pump,Pump_pulse,Coupling,Psi] = Calculate_potential(nx,phi,V0,P0_up,P0_down,P1,P2,theta0,phi0,site0)
Ver=V0*randn(nx,1);
Pump=zeros(2*nx,1);
Pump_pulse=zeros(2*nx,2);
Coupling=1*exp(1j*phi);

Pump(1:nx,1)=P0_up;
Pump(1+nx:2*nx,1)=P0_down;
Pump_pulse(1:nx,1)=P1;
Pump_pulse(1+nx:2*nx,2)=P2;

site=zeros(nx,1);
site(:,1)=1:nx;
Random0=rand(nx,1);
Random1=2*pi*rand(nx,1);
%Random2=2*pi*rand(nx,1);
Random_up=cos(theta0/2).*exp(-1j*(phi0+Random1)-(site-site0).^2/(0.25*nx^2));
Random_down=sin(theta0/2).*exp(-1j*Random1-(site-site0).^2/(0.25*nx^2));
Random_up=Random0.*Random_up;
Random_down=Random0.*Random_down;


Psi(1:nx,1)=Random_up;
Psi(1+nx:2*nx,1)=Random_down;
Model1=sum(abs(Psi(:,1)).^2);
Psi(:,1)=Psi(:,1)./sqrt(Model1);
end

