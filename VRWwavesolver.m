clc
clear all
close all
tic

M=0.5
L=2
sigma=-3

dr=0.001
dt=0.002

p=(dt^2)/(dr^2)

endtime=100
rgrid=2*M:dr:40*M;


VRW= @(r) (1-2*M./r).*((L.*(L+1))./(r.^2)+sigma.*(2*M./(r.^3)))

VRWpot=VRW(rgrid)

h=zeros(2,length(rgrid))

for i=2:endtime/dt
    for j = 2:(rgrid(end)-2*M)/dr
        h(i+1,j)=p*h(i,j+1) + p*h(i,j-1) - h(i,j)*(2*p-2-VRWpot(j))-h(i-1,j);
      
    end
      drawnow
        plot(rgrid,h(1,:))
end

 toc