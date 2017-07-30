clc
clear all
close all

tic
dt=0.001;
dx=0.002;
p=(dt^2)/(dx^2);

endtime=1;

x=0:dx:10;
g(1:length(x))=0 ;
%initial condition of the derivative of h (ie. dh/dt(at t= 0) = g(x)

h=exp(-400*(x-0.3).^2)+1;
% LHBC=0;
% RHBC=0;




for j=2:x(end)/dx
    h(2,j)=1/2*p*h(1,j+1)+(1-p)*h(1,j)+1/2*p*h(1,j-1)+dt*g(j);
    % h(2,1)=LHBC;
    % h(2,end)=RHBC;
    h(2,1)=h(2,1) + dt/dx*(h(2,2)-h(2,1));
    h(2,end)=h(1,end)-dt/dx*(h(1,end)-h(1,end-1));
end

for i=2:endtime/dt
    for j = 2:x(end)/dx
        h(i+1,j)=p*h(i,j+1) + p*h(i,j-1) - h(i,j)*(2*p-2)-h(i-1,j);
      
    end
    % h(i+1,1)=LHBC;
    % h(i+1,end)=RHBC;
    h(i+1,1)=h(i,1) + dt/dx*(h(i,2)-h(i,1));
    h(i+1,end)=h(i,end)-dt/dx*(h(i,end)-h(i,end-1));
%     if mod(i,1)==0
%       drawnow
%       
%         plot(x,h(i,:))
%         axis([0,1, -1,1])
%         disp(i*dt)
%     end
end
% plot(0:dt:endtime, h(:,200))
toc
