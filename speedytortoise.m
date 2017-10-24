clc
clear all
close all
tic
%IMPORTANT NOTE: IF FOR ANY REASON YOU CHANGE M, YOU MUST RERUN R_TORTOISE
%FUNCTION BEFORE USE
load('r_of_x.mat')
rofx=Expression1;%r(x)
%MAKE SURE THIS ACTUALLY SUITS X DOMAIN- SEE MATHEMATICA FILE
%load('m_terp.mat') %misner data
M=1/2;
L=2;
sigma=-3;
%Reggie wheeler parameters and potential eqn:
VRW= @(r) (1-2*M./r).*((L.*(L+1))./(r.^2)+sigma.*(2*M./(r.^3)));



%VRW=@(r) (1-2.*M./r).*(1./((1+3.*M./(2.*r)).^2).*(9.*M.^3./(2.*r.^5)...
%    - (3.*M./(r.^3)).*(1-3.*M./r)) + 6./(r.^(2).*(1+3.*M./(2.*r))));
endtime=1200;
dt=0.005;
dx=0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xint=200;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%computation parameters

p=(dt^2)/(dx^2);

rgrid=-2400*M:dx:3000*M;




%spatial domain

VRWpot=VRW(rofx);
%VRWpot=VRWpot+0.1*exp(-1*(rgrid-5).^2);
t=[0,dt,2*dt];
%initalising time vector


xintindex=find(rgrid==xint);


g(1:length(rgrid))=0;
%initial condition of the derivative of h (ie. dh/dt(at t= 0) = g(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sourcecentre=80;

 ThicknessPara=1;
 Amp=0.001;
% 
% 
% 
% 
 ha= @(x) Amp*exp(-ThicknessPara*(x-sourcecentre).^2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 h=ha(rgrid);%initial condition

%h(1,:)=3.8/0.08*ppval(misnerinterp,rgrid);
% LHBC=0;
% RHBC=0;

%solving for the seccond time point using initial condition of derivative
%(seccond order accuracy)
for j=2:length(rgrid)-2
     h(2,j)=1/2*p*h(1,j+1)+(1-p+(dt^2)/2*VRWpot(j))*h(1,j)+1/2*p*h(1,j-1)+dt*g(j);
   % h(2,j)=dt*g(j)+h(1,j);
    % h(2,1)=LHBC;
    % h(2,end)=RHBC;
    h(2,1)=h(1,1) + dt/dx*(h(1,2)-h(1,1));
    h(2,end)=h(1,end)-dt/dx*(h(1,end)-h(1,end-1));
    %applying radiative BC
end
vectint=h(:,xintindex);
%solving the equation for the bulk (seccond order accuracy)
for i=2:endtime/dt
    t(i+1)=t(i)+dt;
    for j=2:length(rgrid)-1
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
      % generates real time picture of the system
%       drawnow
% %       
%          plot(rgrid,h(3,:))
%          xlabel('x')
%          ylabel('h')
% %        
%          axis([-100,200, -1, 4])
           disp(t(i+1))
       end
    h(1,:)=h(2,:);
    h(2,:)=h(3,:);
    vectint(i+1)=h(2, xintindex);
 
end
figure(2)
plot(t,vectint)
figure(3)
semilogy(t(1:end),abs(vectint))
axis([0,endtime,10^-20,1])

title(strcat('100 test station VRW+0.1*exp(-10*(x-5).^2)'))
toc

str=strcat('VRW_TEST',num2str(rgrid(1)),'_',num2str(rgrid(end)),...
    '_',num2str(endtime),'_',num2str(xint),'_',num2str(sourcecentre),'_',num2str(ThicknessPara),'_',num2str(Amp))
save(str)
load('gong.mat')
soundsc(y)

%syntax: Potential, left bound, right bound, endtime, test station,
% source location, thickness, amplitude

%plots time history at a given point.