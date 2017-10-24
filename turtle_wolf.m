clc
clear all
close all
tic
load('VRWmatrix.mat')
VRWmat=Expression1;
load('usize.mat')
usize=Expression1;
load('step.mat')
step=Expression1;
load('test_station.mat')
samplelocation=Expression1;
load('end.mat')
endterm=Expression1;
sourcecentre=20;
ThicknessPara=1;
Amp=1;

ha= @(x) Amp*exp(-ThicknessPara*(x-sourcecentre).^2);
h=zeros(usize);

u(1,:)=ha(linspace(0,usize*step,usize));
v(1,:)=zeros(usize,1);
h(1,:)=u;
h(:,1)=v';

plot(linspace(0,usize*step,usize),u(1,:))

for i=2:usize-1
    
    t(i+1)=i*step;
    for j=2:length(v)-1
    h(i,j)=h(i,j-1) + h(i-1,j) - h(i-1,j-1) - ((step)^2)/8*VRWmat(i,j)*(h(i,j-1)+h(i-1,j));
    end
   
end

figure(3)
plotything=diag(h);
plot(t,plotything)
diag(VRWmat)
figure(4)
semilogy(t,abs(plotything))
axis([0,t(end),10^-10,10])
toc
str=strcat('TWexpdist_',num2str(endterm),'_',num2str(samplelocation))
save(str)

