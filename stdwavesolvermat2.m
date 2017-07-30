clc
clear all
close all

tic
dt=0.001;
dx=0.002;
p=(dt^2)/(dx^2);

endtime=1;

x=0:dx:20;
g(1:length(x))=0 ;
%initial condition of the derivative of h (ie. dh/dt(at t= 0) = g(x)

h=exp(-400*(x-0.3).^2)+1;
% LHBC=0;
% RHBC=0;

d0=2*(1-p)*ones(1,length(x));
d1=p*ones(1,length(x)-1);
D=diag(d0,0)+diag(d1,1) + diag(d1,-1);

B=zeros(1,length(x));
for j=2:x(end)/dx
    h(2,j)=1/2*p*h(1,j+1)+(1-p)*h(1,j)+1/2*p*h(1,j-1)+dt*g(j);
    % h(2,1)=LHBC;
    % h(2,end)=RHBC;
    h(2,1)=h(2,1) + dt/dx*(h(2,2)-h(2,1));
    h(2,end)=h(1,end)-dt/dx*(h(1,end)-h(1,end-1));
end

for i=2:endtime/dt
%     for j = 2:x(end)/dx
%         h(i+1,j)=p*h(i,j+1) + p*h(i,j-1) - h(i,j)*(2*p-2)-h(i-1,j);
%       
%     end
%      h(i+1,1)=0;
%      h(i+1,end)=0;
    B(1)=h(2,1)*p;
    B(end)=h(2,end)*p;
    h(3,:)=D*h(2,:)'-h(1,:)' + B';
    h(3,1)=h(2,1) + dt/dx*(h(2,2)-h(2,1));
    h(3,end)=h(2,end)-dt/dx*(h(2,end)-h(2,end-1));
%     if mod(i,1)==0
%       drawnow
%       
%         plot(x,h(2,:))
%         axis([0,1, 0,2])
%         disp(i*dt)
%     end
h(1,:)=h(2,:);
h(2,:)=h(3,:);
end
toc