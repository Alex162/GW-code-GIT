clc
clear all
close all
tic
1
M=0.5;
L=2;
sigma=-3;
%Reggie wheeler parameters and potential eqn:
VRW= @(r) (1-2*M./r).*((L.*(L+1))./(r.^2)+sigma.*(2*M./(r.^3)));

endtime=700;
dt=0.01;
dx=0.02;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rint=100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%computation parameters
p=(dt^2)/(dx^2);

rgrid=2*M:dx:900*M;
disp('rgridsize')
disp(size(rgrid))
%spatial domain

VRWpot=VRW(rgrid);
%VRWpot=VRWpot+0.1*exp(-1*(rgrid-5).^2);
t=[0,dt,2*dt];
%initalising time vector



for i=1:length(rgrid)-3
    lhsplus1(i)=(1-p)/12*(1+(VRWpot(i+1)*dt^2)/12);
    lhsminus1(i)=(1-p)/12*(1+(VRWpot(i)*dt^2)/12);
    
    nmainplus1(i)=(1-p)/12*(2+(VRWpot(i+1)*dt^2)/6) + p*(1-(VRWpot(i+1)*dt^2)/12) ...
    -(VRWpot(i+1)*(1-p)*dt^2)/12;
    nmainminus1(i)=(1-p)/12*(2+(VRWpot(i)*dt^2)/6) + p*(1-(VRWpot(i)*dt^2)/12) ...
    -(VRWpot(i)*(1-p)*dt^2)/12;
    
    nprevplus1(i)=-(1-p)/12*(1+(VRWpot(i+1)*(dt^2))/12);
    nprevminus1(i)=-(1-p)/12*(1+(VRWpot(i)*(dt^2))/12);
end

for i=1:length(rgrid)-2
    lhsmid(i)=1 + (VRWpot(i)*dt^2)/12 + (1-p)/12*(-2-(VRWpot(i)*dt^2)/6);
    nmainmid(i)=2 + (VRWpot(i)*dt^2)/6 - (1-p)/12*(4+(VRWpot(i)*dt^2)/3) -2*p...
        + (VRWpot(i)*p*dt^2)/6 - (dt^2)*(VRWpot(i) - ((1-p)/6)*VRWpot(i));
    nprevmid(i)= -1 - (VRWpot(i)*dt^2)/12 + (1-p)/12*(2+(VRWpot(i)*dt^2)/6);
end

disp('computing LHS')
LHS=diag(lhsmid,0)+diag(lhsplus1,1)+diag(lhsminus1,-1);
disp('computing NMAIN')
NMAIN=diag(nmainmid,0)+diag(nmainplus1,1)+diag(nmainminus1,-1);
disp('computing NPREV')
NPREV=diag(nprevmid,0)+diag(nprevplus1,1)+diag(nprevminus1,-1);
size(LHS)
size(NMAIN)
size(NPREV)
BCVEC=zeros(1,length(rgrid)-2);
BCVECprev=zeros(1,length(rgrid)-2);
BCVECnext=zeros(1,length(rgrid)-2);
disp('invertingLHS')
pause(0.5)
invLHS=inv(LHS);
1
toc




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
vectint=h(:,(rint-1)/dx);
%solving the equation for the bulk (seccond order accuracy)
for i=2:endtime/dt
    t(i+1)=t(i)+dt;
    h(3,1)=h(2,1) + dt/dx*(h(2,2)-h(2,1));
    %h(3,1)=2/3*(2*h(2,1)-1/2*h(1,1)+dt/dx*(-3/2*h(2,1)+2*h(2,2)-1/2*h(2,3)));
    h(3,end)=h(2,end) - dt/dx*(h(2,end)-h(2,end-1));
    
    BCVEC(1)=h(2,1)*nmainminus1(1);

    BCVEC(end)=h(2,end)*nmainplus1(end);

    BCVECprev(1)=h(1,1)*nprevminus1(1);

    BCVECprev(end)=h(1,end)*nprevplus1(end);
    
    BCVECnext(1)=h(3,1)*lhsminus1(1);
    
    BCVECnext(end)=h(3,end)*lhsplus1(end);
    
    h(3,2:end-1)=(NMAIN*h(2,2:end-1)' + NPREV*h(1,2:end-1)' + BCVEC' + BCVECprev'-BCVECnext');
    h(3,2:end-1)=invLHS*h(3,2:end-1)';

% 
%     if mod(i,10)==0
%       drawnow
%       
%         plot(rgrid,h(2,:))
%         %axis([rgrid(1)-1,rgrid(end)+1, -1, 1])
%         disp(i*dt)
%     end

h(1,:)=h(2,:);
h(2,:)=h(3,:);
   
    

    vectint(i+1)=h(2, (rint-1)/dx);
    
end
% figure(2)
% plot(t,vectint)
% figure(3)
% semilogy(t,abs(vectint))
% axis([0,endtime,10^-20,1])
% annotation('textbox','String',[.0,.0,.3,.3],'String',strcat('Amp=',num2str(Amp),...
%     ' ThicknessPara=',num2str(ThicknessPara),' sourcecentre=',num2str(sourcecentre)),'FitBoxToText','on')
% title(strcat(char(ha),'100 test station VRW+0.1*exp(-10*(x-5).^2)'))
toc

str=strcat('VRWx_xm2tm2_',num2str(rgrid(end)),'_',num2str(endtime),'_',num2str(rint),...
    '_',num2str(sourcecentre),'_',num2str(ThicknessPara),'_',num2str(Amp))
save(str)

%syntax: Potential, domain size, endtime, test station, source location, thickness,
%amplitude

%plots time history at a given point.