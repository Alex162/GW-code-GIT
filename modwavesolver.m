clc
clear all
close all
tic
M=0.5;
L=2;
sigma=-3;
%Reggie wheeler parameters and potential eqn:
VRW= @(r) (1-2*M./r).*((L.*(L+1))./(r.^2)+sigma.*(2*M./(r.^3)));

endtime=100;
dt=0.005;
dx=0.0;
%computation parameters
p=(dt^2)/(dx^2);

rgrid=2*M:dx:40*M;
%spatial domain

VRWpot=VRW(rgrid);

t=[0,dt,2*dt];
%initalising time vector




g(1:length(rgrid))=0 ;
%initial condition of the derivative of h (ie. dh/dt(at t= 0) = g(x)

h=exp(-(rgrid-15).^2);%initial condition
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

%solving the equation for the bulk (seccond order accuracy)
for i=2:endtime/dt
    t(i+1)=t(i)+dt;
    for j=2:(rgrid(end)-2*M)/dx
        h(i+1,j)=p*h(i,j+1) + p*h(i,j-1) - h(i,j)*(2*p-2+dt^2*VRWpot(j))-h(i-1,j);
      
    end
    % h(i+1,1)=LHBC;
    % h(i+1,end)=RHBC;
    h(i+1,1)=h(i,1) + dt/dx*(h(i,2)-h(i,1));
    h(i+1,end)=h(i,end)-dt/dx*(h(i,end)-h(i,end-1));
%     h(i+1,end)=(1/dt^2-1/dt)^-1 * (h(i,end)*(2/(dt^2)-1/dt+VRWpot(end)+1/(dx)^2+1/dx)...
%         + h(i,end-1)*(-2/(dx^2)-1/dx) + h(i,end-2)/(dx^2) - h(i-1,end)/(dt^2));
%     if mod(i,5)==0
%         %generates real time picture of the system
%       drawnow
%       
%         plot(rgrid,h(i+1,:))
%         axis([rgrid(1),rgrid(end), -1.5, 1.5])
%         disp(t(i+1))
%     end
end
figure(2)
plot(t,h(:,(17-1)/dx))
figure(3)
semilogy(t,abs(h(:,(17-1)/dx)))
% plots the time history at 17. the presence of the minus one is to account for no storage of information within the BH
axis([0,endtime,10^-10,1])
%plots time history at a given point.
toc