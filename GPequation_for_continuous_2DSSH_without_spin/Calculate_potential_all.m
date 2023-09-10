function [SpaceX,SpaceY,Vpo,P,Vortex_Pulse,Vortex_Pump1,Vortex_Pump2] = Calculate_potential_all(a0,nx,ny,Boundary,Boundarx,redus,sitea,Xnumber,Ynumber,V0,P0,P1)
%UNTITLED3 此处提供此函数的摘要
%   此处提供详细说明

Vpo=V0*ones(Ynumber,Xnumber);
P=zeros(Ynumber,Xnumber);
Vortex_Pump1=zeros(Ynumber,Xnumber);Vortex_Pump2=zeros(Ynumber,Xnumber);
Vortex_Pulse=zeros(Ynumber,Xnumber);
SpaceX=zeros(Ynumber,Xnumber);
SpaceY=zeros(Ynumber,Xnumber);

for i=1:Xnumber
    for j=1:Ynumber
        x=(nx*a0+2*Boundarx)*(i-1)/(Xnumber);
        y=(ny*a0+2*Boundary)*(j-1)/(Ynumber);
        SpaceX(j,i)=x;
        SpaceY(j,i)=y;
%         for ysite=0:ny-1
%             for xsite=0:nx-1
%                 if ((x-a0/2-sitea-xsite*a0-Boundarx)^2+(y-a0/2-sitea-ysite*a0-Boundary)^2)<((1*redus)^2)
%                    Vpo(j,i)=0;
%                 elseif ((x-a0/2+sitea-xsite*a0-Boundarx)^2+(y-a0/2+sitea-ysite*a0-Boundary)^2)<((1*redus)^2)
%                    Vpo(j,i)=0;
%                 elseif ((x-a0/2-sitea-xsite*a0-Boundarx)^2+(y-a0/2+sitea-ysite*a0-Boundary)^2)<((1*redus)^2)
%                    Vpo(j,i)=0;
%                 elseif ((x-a0/2+sitea-xsite*a0-Boundarx)^2+(y-a0/2-sitea-ysite*a0-Boundary)^2)<((1*redus)^2)
%                    Vpo(j,i)=0;
%                 end
% %                 Vpo(j,i)=Vpo(j,i)-P0*exp(-((x-a0/2-sitea-xsite*a0-Boundarx)^2+(y-a0/2-sitea-ysite*a0-Boundary)^2)/((0.7*redus)^2));
% %                 Vpo(j,i)=Vpo(j,i)-P0*exp(-((x-a0/2+sitea-xsite*a0-Boundarx)^2+(y-a0/2+sitea-ysite*a0-Boundary)^2)/((0.7*redus)^2));
% %                 Vpo(j,i)=Vpo(j,i)-P0*exp(-((x-a0/2-sitea-xsite*a0-Boundarx)^2+(y-a0/2+sitea-ysite*a0-Boundary)^2)/((0.9*redus)^2));
% %                 Vpo(j,i)=Vpo(j,i)-P0*exp(-((x-a0/2+sitea-xsite*a0-Boundarx)^2+(y-a0/2-sitea-ysite*a0-Boundary)^2)/((1*redus)^2));
% %                 P(j,i)=P(j,i)+P0*exp(-((x-a0/2+sitea-xsite*a0-Boundarx-redus*cos(theta)/2)^2+ ...
% %                     (y-a0/2-sitea-ysite*a0-Boundary-redus*sin(theta)/2)^2)/((0.5*redus)^2));
% %                 P(j,i)=P(j,i)+0*gamma*exp(-((x-a0/2-sitea-xsite*a0-Boundarx)^2+(y-a0/2-sitea-ysite*a0-Boundary)^2)/((1*redus)^2));
% %                 P(j,i)=P(j,i)+0*gamma*exp(-((x-a0/2+sitea-xsite*a0-Boundarx)^2+(y-a0/2+sitea-ysite*a0-Boundary)^2)/((1*redus)^2));
% %                 P(j,i)=P(j,i)+0*gamma*exp(-((x-a0/2-sitea-xsite*a0-Boundarx)^2+(y-a0/2+sitea-ysite*a0-Boundary)^2)/((1*redus)^2));
% %                 P(j,i)=P(j,i)+1*gamma*exp(-((x-a0/2+sitea-xsite*a0-Boundarx)^2+(y-a0/2-sitea-ysite*a0-Boundary)^2)/((1*redus)^2));
%             end
%         end
    end
end

Width=1*(1*redus-0.4*redus);
center_r=(1*redus+0.4*redus)/2;
redus1=1*redus;
for ysite=0:ny-1
    for xsite=0:nx-1
        [phi,r]=cart2pol(SpaceX-a0/2+sitea-xsite*a0-Boundarx,SpaceY-a0/2-sitea-ysite*a0-Boundary);
        %P=P+P0*exp(-(r.^2)/(redus1^2));%abs(exp(-((r-center_r)./Width).^2)).^2;
        Vpo=Vpo-V0*exp(-(r.^20)/(redus^20));

        [phi,r]=cart2pol(SpaceX-a0/2-sitea-xsite*a0-Boundarx,SpaceY-a0/2-sitea-ysite*a0-Boundary);
        %P=P+P0*exp(-(r.^2)/(redus1^2));
        Vpo=Vpo-V0*exp(-(r.^20)/(redus^20));

        [phi,r]=cart2pol(SpaceX-a0/2+sitea-xsite*a0-Boundarx,SpaceY-a0/2+sitea-ysite*a0-Boundary);
        %P=P+P0*exp(-(r.^2)/(redus1^2));
        Vpo=Vpo-V0*exp(-(r.^20)/(redus^20));

        [phi,r]=cart2pol(SpaceX-a0/2-sitea-xsite*a0-Boundarx,SpaceY-a0/2+sitea-ysite*a0-Boundary);
        %P=P+P0*exp(-(r.^2)/(redus1^2));
        Vpo=Vpo-V0*exp(-(r.^20)/(redus^20));
    end
end
% for ysite=1:ny-2
%     for xsite=1:nx-2
%         [phi,r]=cart2pol(SpaceX-a0/2+sitea-xsite*a0-Boundarx,SpaceY-a0/2-sitea-ysite*a0-Boundary);
%         P=P-P0*exp(-(r.^2)/(redus1^2));%abs(exp(-((r-center_r)./Width).^2)).^2;
% 
%         [phi,r]=cart2pol(SpaceX-a0/2-sitea-xsite*a0-Boundarx,SpaceY-a0/2-sitea-ysite*a0-Boundary);
%         P=P-P0*exp(-(r.^2)/(redus1^2));
% 
%         [phi,r]=cart2pol(SpaceX-a0/2+sitea-xsite*a0-Boundarx,SpaceY-a0/2+sitea-ysite*a0-Boundary);
%         P=P-P0*exp(-(r.^2)/(redus1^2));
% 
%         [phi,r]=cart2pol(SpaceX-a0/2-sitea-xsite*a0-Boundarx,SpaceY-a0/2+sitea-ysite*a0-Boundary);
%         P=P-P0*exp(-(r.^2)/(redus1^2));
%     end
% end
% xsite=0;ysite=0;
% [phi,r]=cart2pol(SpaceX-a0/2+sitea-xsite*a0-Boundarx,SpaceY-a0/2+sitea-ysite*a0-Boundary);
% P=P+P0*abs(exp(-((r-center_r)./Width).^2)).^2;
% xsite=nx-1;ysite=0;
% [phi,r]=cart2pol(SpaceX-a0/2-sitea-xsite*a0-Boundarx,SpaceY-a0/2+sitea-ysite*a0-Boundary);
% P=P+P0*abs(exp(-((r-center_r)./Width).^2)).^2;
% xsite=nx-1;ysite=ny-1;
% [phi,r]=cart2pol(SpaceX-a0/2-sitea-xsite*a0-Boundarx,SpaceY-a0/2-sitea-ysite*a0-Boundary);
% P=P+P0*abs(exp(-((r-center_r)./Width).^2)).^2;
%LU
xsite=0;ysite=ny-1;
[phi,r]=cart2pol(SpaceX-a0/2+sitea-xsite*a0-Boundarx,SpaceY-a0/2-sitea-ysite*a0-Boundary);
Vortex_Pump1(:,:)=Vortex_Pump1(:,:)+exp(1j*phi).*exp(-((r-center_r)./(Width)).^2);
Vortex_Pump2(:,:)=Vortex_Pump2(:,:)+exp(1j*phi).*exp(-((r-center_r)./(Width)).^2);
Vortex_Pulse(:,:)=Vortex_Pulse(:,:)+P1*exp(-((r-center_r)./Width).^2).^2;%(sin(phi-3*pi/4).^2).*
P=P+P0*exp(-((r-center_r)./Width).^2).^2;
%RU
% xsite=nx-1;ysite=ny-1;
% [phi,r]=cart2pol(SpaceX-a0/2-sitea-xsite*a0-Boundarx,SpaceY-a0/2-sitea-ysite*a0-Boundary);
% % Vortex_Pump(:,:)=Vortex_Pump(:,:)+0.5*exp(mru*1j*phi).*exp(-((r-center_r)./(Width)).^2);
% Vortex_Pulse(:,:)=Vortex_Pulse(:,:)+P1*(sin(phi-3*pi/4).^2).*exp(-((r-center_r)./Width).^2).^2;%(sin(phi-3*pi/4).^2).*
%LB
% xsite=0;ysite=0;
% [phi,r]=cart2pol(SpaceX-a0/2+sitea-xsite*a0-Boundarx,SpaceY-a0/2+sitea-ysite*a0-Boundary);
% Vortex_Pump(:,:)=Vortex_Pump(:,:)+0.5*exp(mlb*1j*phi).*exp(-((r-center_r)./(Width)).^2);
% Vortex_Pulse(:,:)=Vortex_Pulse(:,:)+P1*(sin(phi-0*pi/4).^2).*exp(-((r-center_r)./Width).^2).^2;
%RB
xsite=nx-1;ysite=0;
[phi,r]=cart2pol(SpaceX-a0/2-sitea-xsite*a0-Boundarx,SpaceY-a0/2+sitea-ysite*a0-Boundary);
Vortex_Pump1(:,:)=Vortex_Pump1(:,:)+exp(1j*phi).*exp(-((r-center_r)./(Width)).^2);
Vortex_Pump2(:,:)=Vortex_Pump2(:,:)+exp(-1j*phi).*exp(-((r-center_r)./(Width)).^2);
Vortex_Pulse(:,:)=Vortex_Pulse(:,:)+P1*exp(-((r-center_r)./Width).^2).^2;%(sin(phi-1*pi/4).^2).*
P=P+P0*exp(-((r-center_r)./Width).^2).^2;

% xsite=floor(nx/2);ysite=ny-1;
% [phi,r]=cart2pol(SpaceX-a0/2-sitea-xsite*a0-Boundarx,SpaceY-a0/2-sitea-ysite*a0-Boundary);
% Vortex_Pump(:,:)=Vortex_Pump(:,:)+P1*exp(-1j*phi).*exp(-((r-center_r)./(Width)).^2);
% Vortex_Pulse=Vortex_Pulse+P1*exp(-((r-center_r)./Width).^2).^2;
% xsite=nx-1;ysite=0;
% [phi,r]=cart2pol(SpaceX-a0/2-sitea-xsite*a0-Boundarx,SpaceY-a0/2+sitea-ysite*a0-Boundary);
% Vortex_Pump(:,:,2)=exp(1i*phi).*exp(-((r-center_r)./(Width)).^2);
% Vortex_Pulse=Vortex_Pulse+P1*exp(-((r-center_r)./Width).^2).^2;
%Vpo=Vpo-5*SpaceX/Lx-5*SpaceY/Ly;

end