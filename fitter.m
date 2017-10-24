clc
clear all
close all
tic
% 
% x = 70:0.1:100;
% y = (x-7).^-7
% 
% fo = fitoptions('Method','NonlinearLeastSquares',...
%                'Lower',[-Inf,-Inf,-8],...
%                'Upper',[Inf,Inf,-6],...
%                'StartPoint',[1 1 1]);
% ft = fittype('a*(x-b)^-n','independent',{'x'},...
%     'coefficients',{'a','b','n'},'options',fo);
% 
% myfit=fit(x',y',ft)
% plot(myfit,x,y)
sourcecentre=0;
load('VZerT_misnerXimp_-600_1500_1200_200_')
xint=rint
%syntax:Potential, possibly mesh change indicator, domain size, endtime, test station, source location, thickness,
%amplitude
tbracket(1)=xint+sourcecentre+150;

tbracket(2)=tbracket(1)+400;



minimizermin=200000000000;

% logh=log10(abs(vectint));
% t=t;
% plot(t,logh)
% axis([0,endtime,-20,1])
% [xgin,ygin]=ginput(2);
% tic
% xgin=sort(xgin);
 for i=1:2
     A=tbracket(i)-t;
     Asorted=sort(abs(A));
 
 
     index(i)=find(A==Asorted(i) | A==-Asorted(i));
     
 end

 figure(1)
 semilogy(t,abs(vectint))
 hold on
 semilogy(t(index(1)),abs(vectint(index(1))),'*')
 semilogy(t(index(2)),abs(vectint(index(2))),'*')
 axis([0,endtime,10^-20,1])
 drawnow
for n=0:0.1:xint+sourcecentre
%for n=360
tclipped=t(index(1):index(2))-n;

ptail=vectint(index(1):index(2));






[fitfunct,S]=polyfit(log10(abs(tclipped)),log10(abs(ptail))',1);
[fit,delta]=polyval(fitfunct,log10(abs(tclipped)),S);
minimizer=norm(delta);

if abs(minimizer)<abs(minimizermin)
    minimizermin=minimizer;
    fitmin=fit;
    deltamin=delta;
    nmin=n;
    plawindexmin=fitfunct(1);
end
end
disp('Power law index is ')
disp(plawindexmin)
disp('Shift factor is')
disp(nmin)
disp('minimizer is')
disp(minimizermin)
figure(2)
plot(log10(tclipped),fitmin)
hold on
plot(log10(tclipped),log10(abs(ptail)))
legend('fit','data')
toc