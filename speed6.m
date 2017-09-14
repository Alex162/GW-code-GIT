clc
clear all
close all
tic
M=0.5;
L=2;
sigma=-3;
%Reggie wheeler parameters and potential eqn:
VRW= @(r) (1-2*M./r).*((L.*(L+1))./(r.^2)+sigma.*(2*M./(r.^3)));

endtime=400;
dt=0.00005;
dx=0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rint=100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%computation parameters
p=(dt^2)/(dx^2);

rgrid=2*M:dx:600*M;
%spatial domain

VRWpot=VRW(rgrid);
%VRWpot=VRWpot+0.1*exp(-1*(rgrid-5).^2);
t=[0,dt,2*dt];
%initalising time vector




g(1:length(rgrid))=0 ;
%initial condition of the derivative of h (ie. dh/dt(at t= 0) = g(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sourcecentre=80;

ThicknessPara=1;
Amp=1;
ha= @(x) Amp*exp(-ThicknessPara*(x-sourcecentre).^2);
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
for i=2:6
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

%     if mod(i,30)==0
%         %generates real time picture of the system
%       drawnow
%       
%         plot(rgrid,h(3,:))
%         axis([rgrid(1),rgrid(end), -0.015, 0.015])
%         disp(t(i+1))
%     end

end
    vectint(:)=h(:, (rint-1)/dx);
    

%solving the equation for the bulk (seccond order accuracy)
for i=7:endtime/dt
    t(i+1)=t(i)+dt;
    for j=4:(rgrid(end)-2*M)/dx -3
        h(8,j)=-90/469 * ( (-223/(10)+49*dt^2/(18*dx^2)-dt^2*VRWpot(j))*h(7,j)...
            +879/20*h(6,j) -949/18*h(5,j) +41*h(4,j) -201/10*h(3,j) +1019/180*h(2,j)... 
            -7/10*h(1,j)    +p*(-1/90*h(7,j-3) +3/20*h(7,j-2) -3/2*h(7,j-1)...
            -3/2*h(7,j+1) +3/20*h(7,j+2)-1/90*h(7,j+3)));
      
    end
        for j=2:3
        h(8,j)=-90/469 * ( (-223/(10)+49*dt^2/(18*dx^2)-dt^2*VRWpot(j))*h(7,j)...
            +879/20*h(6,j) -949/18*h(5,j) +41*h(4,j) -201/10*h(3,j) +1019/180*h(2,j)... 
            -7/10*h(1,j)    +p*(223/10*h(7,j+1) -879/20*h(7,j+2) + 949/18*h(7,j+3)...
            -41*h(7,j+4) +201/10*h(7,j+5)-1019/180*h(7,j+6)));
        end
        
      for j=(rgrid(end)-2*M)/dx-2:(rgrid(end)-2*M)/dx-1
        -90/469 * ( (-223/(10)+49*p/18-dt^2*VRWpot(j))*h(7,j)...
            +879/20*h(6,j) -949/18*h(5,j) +41*h(4,j) -201/10*h(3,j) +1019/180*h(2,j)... 
            -7/10*h(1,j)    +p*(223/10*h(7,j-1) -879/20*h(7,j-2) + 949/18*h(7,j-3)...
            -41*h(7,j-4) +201/10*h(7,j-5)-1019/180*h(7,j-6)));
        end
    % h(i+1,1)=LHBC;
    % h(i+1,end)=RHBC;
    h(8,1)=h(7,1) + dt/dx*(h(7,2)-h(7,1));
    h(8,end)=h(7,end)-dt/dx*(h(7,end)-h(7,end-1));
%     h(i+1,end)=(1/dt^2-1/dt)^-1 * (h(i,end)*(2/(dt^2)-1/dt+VRWpot(end)+1/(dx)^2+1/dx)...
%         + h(i,end-1)*(-2/(dx^2)-1/dx) + h(i,end-2)/(dx^2) - h(i-1,end)/(dt^2));
if mod(i,1000)==0
        %generates real time picture of the system
      drawnow
      
        plot(rgrid,h(3,:))
        axis([rgrid(1),rgrid(end), -1, 1])
        disp(t(i+1))
end
    h(1,:)=h(2,:);
    h(2,:)=h(3,:);
    h(3,:)=h(4,:);
    h(4,:)=h(5,:);
    h(5,:)=h(6,:);
    h(6,:)=h(7,:);
    h(7,:)=h(8,:);
   
    
    
    
    vectint(i+1)=h(7, (rint-1)/dx);
    
end
figure(2)
plot(t,vectint)
figure(3)
semilogy(t,abs(vectint))
axis([0,endtime,10^-20,1])
annotation('textbox','String',[.0,.0,.3,.3],'String',strcat('Amp=',num2str(Amp),...
    ' ThicknessPara=',num2str(ThicknessPara),' sourcecentre=',num2str(sourcecentre)),'FitBoxToText','on')
title(strcat(char(ha),'100 test station VRW+0.1*exp(-10*(x-5).^2)'))
toc

str=strcat('VRWBp1w_',num2str(rgrid(end)),'_',num2str(endtime),'_',num2str(rint),...
    '_',num2str(sourcecentre),'_',num2str(ThicknessPara),'_',num2str(Amp))
save(str)

%syntax: Potential, domain size, endtime, test station, source location, thickness,
%amplitude

%plots time history at a given point.