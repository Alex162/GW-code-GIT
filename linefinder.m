clc
clear all
close all

load('VRWT_900_1200_400_80_1_1')
%syntax:Potential, possibly mesh change indicator, domain size, endtime, test station, source location, thickness,
%amplitude


permissibleerror=10^-7

logh=log10(vectint);
t=t;
plot(t,logh)
axis([0,endtime,-20,1])
[xgin,ygin]=ginput(3);
tic
xgin=sort(xgin);
for i=1:3
    A=xgin(i)-t;
    Asorted=sort(abs(A));


    indexy=find(A==Asorted(i) | A==-Asorted(i));
    ygin(i)=logh(indexy);
end
xgin;
ygin;

xleft=xgin(1)
ylogleft=ygin(1)

xmid=xgin(2)
ylogmid=ygin(2)

xright=xgin(3)
ylogright=ygin(3)
residual=inf
residualmin=inf


for n=0:0.1:endtime
xlogleft=log10(xleft-n);
xlogmid=log10(xmid-n);
xlogright=log10(xright-n);

p=polyfit([xlogleft,xlogright],[ylogleft,ylogright],1);
fitfunct= p(1)*xlogmid + p(2);

residual=abs(fitfunct-ylogmid);

if residual < permissibleerror
    nmin=n;
    residualmin=residual;
    plindex=p(1);
    disp('Sequence broken, permissible error reached')
   break
end
if residual<residualmin;
    nmin=n;
    residualmin=residual;
    plindex=p(1);
end


end

toc
semilogx(t-nmin,logh)
hold on
semilogx(xleft-nmin,ylogleft,'*')
semilogx(xmid-nmin,ylogmid,'*')
semilogx(xright-nmin,ylogright,'*')
axis([10,endtime,-20,1])

disp('minimum residual is ')
disp(residualmin)
disp('subtraction term is ')
disp(nmin)
disp('power law index is ')
disp(plindex)