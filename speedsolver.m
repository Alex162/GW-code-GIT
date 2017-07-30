clc
clear all
close all
tic
M=0.5;
L=2;
sigma=-3;
%Reggie wheeler parameters and potential eqn:
VRW= @(r) (1-2*M./r).*((L.*(L+1))./(r.^2)+sigma.*(2*M./(r.^3)));

endtime=200;
dt=0.005;
dx=0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rint=100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%computation parameters
p=(dt^2)/(dx^2);

rgrid=2*M:dx:400*M;
%spatial domain

VRWpot=VRW(rgrid);

t=[0,dt,2*dt];
%initalising time vector




g(1:length(rgrid))=0 ;
%initial condition of the derivative of h (ie. dh/dt(at t= 0) = g(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ha= @(x) 1*exp(-1*(x-15).^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=ha(rgrid);%initial condition
% LHBC=0;
% RHBC=0;

%solving for the seccond time point using initial condition of derivative
%(seccond order accuracy)
for j=2:(rgrid(end)-2*M)/dx
     h(2,j)=1/2*p*h(1,j+1)+(1-p+(dt^2)/2*VRWpot(j))*h(1,j)+1/2*p*h(1,j-1)+dt*g(j);
   % h(2,j)=dt*g(j)+h(1,j);
    % h(2,1)=LHBC;
    % h(2,end)=RHBC;
    h(2,1)=h(2,1) + dt/dx*(h(2,2)-h(2,1));
    h(2,end)=h(1,end)-dt/dx*(h(1,end)-h(1,end-1));
    %applying radiative BC
end
vectint=h(:,rint);
%solving the equation for the bulk (seccond order accuracy)
for i=2:endtime/dt
    t(i+1)=t(i)+dt;
    for j=2:(rgrid(end)-2*M)/dx
        h(3,j)=p*h(2,j+1) + p*h(2,j-1) - h(2,j)*(2*p-2+dt^2*VRWpot(j))-h(1,j);
      
    end
    % h(i+1,1)=LHBC;
    % h(i+1,end)=RHBC;
    h(3,1)=h(3,1) + dt/dx*(h(2,2)-h(2,1));
    h(3,end)=h(2,end)-dt/dx*(h(2,end)-h(2,end-1));
%     h(i+1,end)=(1/dt^2-1/dt)^-1 * (h(i,end)*(2/(dt^2)-1/dt+VRWpot(end)+1/(dx)^2+1/dx)...
%         + h(i,end-1)*(-2/(dx^2)-1/dx) + h(i,end-2)/(dx^2) - h(i-1,end)/(dt^2));

%     if mod(i,30)==0
%         %generates real time picture of the system
%       drawnow
%       
%         plot(rgrid,h(3,:))
%         axis([rgrid(1),rgrid(end), -0.015, 0.015])
%         disp(t(i+1))
%     end
    h(1,:)=h(2,:);
    h(2,:)=h(3,:);
    vectint(i+1)=h(2, (rint-1)/dx);
    
end
figure(2)
plot(t,vectint)
figure(3)
semilogy(t,abs(vectint))
axis([0,endtime,10^-20,1])

title(strcat(char(ha),'100 test station'))
toc
%plots time history at a given point.