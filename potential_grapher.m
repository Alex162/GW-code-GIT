clc
clear all
close all

% G=6.67*10^-11
% c=3*10^8
% 
% geoconstmass=G/(c^2)
% Msungeo=2*10^30*geoconstmass
% M=30*Msungeo
M=0.5
L=2
sigma=-3


VRW= @(r) (1-2*M./r).*((L.*(L+1))./(r.^2)+sigma.*(2*M./(r.^3)))

%tort= @(r) 1+2*M*log(r./(2*M)-1)

r= @(x) 2*M*lambertw(exp(x./((2*M))-1))+2*M


x=[-10:0.01:20];


 plot(x,VRW(r(x)))


figure(2)

hold on

x=2*M:0.01:1000*M
VRWx=VRW(x)
plot(x,VRW(x))
title('VRW vs radius')
figure(3)
semilogy(x,VRWx)
title('VRW vs radius')

