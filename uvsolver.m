clc
clear all
close all
tic
M=0.5
L=2
sigma=-3
load('r_of_x.mat')%from mathematica
rofx=Expression1;
VRW= @(r) (1-2*M./r).*((L.*(L+1))./(r.^2)+sigma.*(2*M./(r.^3)));

%r= @(x) 2*M*lambertw(exp(x./((2*M))-1))+2*M;
dx=.1
uvstep=sqrt(2*dx)
x=0:dx:600;
x2=-1000:dx:1000;

sourcecentre=100;

ThicknessPara=1;
Amp=1;

ha= @(x) Amp*exp(-ThicknessPara*(x-sourcecentre).^2);
h=zeros(length(x));
u(1,:)=ha(x);
v(1,:)=zeros(length(x),1);
h(1,:)=u;
h(:,1)=v';
xint=200;
xintindex=find(x==xint)
plotvect(1)=v(1)
potindex=find(x2==xint)
for i=2:length(u)-1
    
    t(i+1)=i*dx;
    for j=2:length(v)-1
    h(i,j)=h(i,j-1) + h(i-1,j) - h(i-1,j-1) - ((uvstep)^2)/8*VRW(rofx(potindex+j-i))*(h(i,j-1)+h(i-1,j));
    end
   
end
% figure(1)
% contourf(h)
% colorbar


figure(3)
plotything=diag(h);
plot(t,plotything)

figure(4)
semilogy(2*t,abs(plotything))
axis([500,700,10^-10,1])
toc
% 
% figure(5)
% loglog(2*t,abs(plotything))
% axis([1,2*t(end),10^-20,1])
% for p=1:length(h)
% if mod(p,5)==0
%     figure(5)
% drawnow
%     plot(x2,VRW(rofx))
%     hold on
%     plot(x(1:end),h(p,:))
%     hold off
% end
% end


toc
