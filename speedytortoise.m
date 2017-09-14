clc
clear all
close all
tic
M=0.5;
L=2;
sigma=-3;
%Reggie wheeler parameters and potential eqn:
VRW= @(r) (1-2*M./r).*((L.*(L+1))./(r.^2)+sigma.*(2*M./(r.^3)));
tortoiser= @(x) 2*M*lambertw(exp(x./((2*M))-1))+2*M;

endtime=1200;
dt=0.005;
dx=0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rint=400;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%computation parameters

p=(dt^2)/(dx^2);

rgrid=-40*M:dx:1800*M;
tortgrid=tortoiser(rgrid);
%spatial domain

VRWpot=VRW(tortgrid);
%VRWpot=VRWpot+0.1*exp(-1*(rgrid-5).^2);
t=[0,dt,2*dt];
%initalising time vector

rinttort=tortoiser(rint);
rintindex=find(tortgrid==rinttort);


g(1:length(rgrid))=0;
%initial condition of the derivative of h (ie. dh/dt(at t= 0) = g(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sourcecentre=400;

ThicknessPara=1;
Amp=1;
ha= @(x) Amp*exp(-ThicknessPara*(x-sourcecentre).^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=ha(tortgrid);%initial condition
% LHBC=0;
% RHBC=0;

%solving for the seccond time point using initial condition of derivative
%(seccond order accuracy)
for j=2:length(tortgrid)-2
     h(2,j)=1/2*p*h(1,j+1)+(1-p+(dt^2)/2*VRWpot(j))*h(1,j)+1/2*p*h(1,j-1)+dt*g(j);
   % h(2,j)=dt*g(j)+h(1,j);
    % h(2,1)=LHBC;
    % h(2,end)=RHBC;
    h(2,1)=h(2,1) + dt/dx*(h(2,2)-h(2,1));
    h(2,end)=h(1,end)-dt/dx*(h(1,end)-h(1,end-1));
    %applying radiative BC
end
vectint=h(:,rintindex);
%solving the equation for the bulk (seccond order accuracy)
for i=2:endtime/dt
    t(i+1)=t(i)+dt;
    for j=2:length(tortgrid)-1
        h(3,j)=p*h(2,j+1) + p*h(2,j-1) - h(2,j)*(2*p-2+dt^2*VRWpot(j))-h(1,j);
      
    end
    % h(i+1,1)=LHBC;
    % h(i+1,end)=RHBC;
    h(3,1)=h(2,1) + dt/dx*(h(2,2)-h(2,1));
    %h(3,1)=2/3*(2*h(2,1)-1/2*h(1,1)+dt/dx*(-3/2*h(2,1)+2*h(2,2)-1/2*h(2,3)));
    h(3,end)=h(2,end)-dt/dx*(h(2,end)-h(2,end-1));
%     h(i+1,end)=(1/dt^2-1/dt)^-1 * (h(i,end)*(2/(dt^2)-1/dt+VRWpot(end)+1/(dx)^2+1/dx)...
%         + h(i,end-1)*(-2/(dx^2)-1/dx) + h(i,end-2)/(dx^2) - h(i-1,end)/(dt^2));
% 
       if mod(i,200)==0
%       % generates real time picture of the system
%       drawnow
%       
%         plot(tortgrid,h(3,:))
%         xlabel('x')
%         ylabel('h')
%        
%         axis([tortgrid(1)-20,tortgrid(end)+20, -1, 1])
            disp(t(i+1))
       end
    h(1,:)=h(2,:);
    h(2,:)=h(3,:);
    vectint(i+1)=h(2, rintindex);
 
end
figure(2)
plot(t,vectint)
figure(3)
semilogy(t(1:end),abs(vectint))
axis([0,endtime,10^-20,1])
annotation('textbox','String',[.0,.0,.3,.3],'String',strcat('Amp=',num2str(Amp),...
    ' ThicknessPara=',num2str(ThicknessPara),' sourcecentre=',num2str(sourcecentre)),'FitBoxToText','on')
title(strcat(char(ha),'100 test station VRW+0.1*exp(-10*(x-5).^2)'))
toc

str=strcat('VRWT_',num2str(rgrid(end)),'_',num2str(endtime),'_',num2str(rint),...
    '_',num2str(sourcecentre),'_',num2str(ThicknessPara),'_',num2str(Amp))
save(str)

%syntax: Potential, domain size, endtime, test station, source location, thickness,
%amplitude

%plots time history at a given point.