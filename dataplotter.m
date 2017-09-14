clc
clear 
close all
% 
load('VRWBCspesh_500_800_100_80_1_1')
i100s=vectint;
t100s=t;
%syntax:Potential, possibly mesh change indicator, domain size, endtime, test station, source location, thickness,
%amplitude


load('VRWBCspesh_1200_1600_400_380_1_1')
i400s=vectint;
t400s=t;
%syntax:Potential, possibly mesh change indicator, domain size, endtime, test station, source location, thickness,
%amplitude


load('VRW_2000_2400_800_780_1_1')
i800=vectint;
t800=t;
%syntax:Potential, possibly mesh change indicator, domain size, endtime, test station, source location, thickness,
%amplitude
figure(2)



load('VRW_3600_4000_1600_1580_1_1')
i1600=vectint;
t1600=t;
%syntax:Potential, possibly mesh change indicator, domain size, endtime, test station, source location, thickness,
%amplitude


figure(3)

%loglog(t100s,abs(i100s))

% loglog(t400s-79-600, abs(i400s))
% loglog(t800-1400-80, abs(i800))
% loglog(t1600-3000-80, abs(i1600))
 loglog(t100s ,10^7.5*(t100s).^-7)
hold on
%loglog(t100s-110,abs(i100s))

 loglog(t100s,10^7.5*(t100s-20).^-7)
 loglog(t100s,10^7.5*(t100s-50).^-7)
axis([1,800,10^-20,1])
title('(t).^-7, (t-20).^-7, (t-50).^-7 vs t')
xlabel('t')
ylabel('h')
legend('(t).^-7','(t-20).^-7','(t-50).^-7')
%annotation('textbox','String',[.0,.0,.3,.3],'String',strcat('Amp=',num2str(Amp),...
% title(strcat(char(ha),'100 test station 0.1Bw'))
% hold on
% 
% load('VRW_900_800_100_400_1_1')
% sc400=vectint;
% tsc400=t;
% load('VRW_400_400_100_120_1_1')
% sc120=vectint;
% tsc120=t;
% load('VRW_900_800_100_80_1_1')
% sc80=vectint;
% tsc80=t;
% load('VRW_700_600_100_15_1_1')
% sc15=vectint;
% tsc15=t;
%
% 
% 
% %syntax:Potential, possibly mesh change indicator, domain size, endtime, test station, source location, thickness,
% %amplitude
% figure(1)
% 
% loglog(tsc400-385,abs(sc400))
% hold on
% loglog(tsc120-105,abs(sc120))
% loglog(tsc80-65,abs(sc80))
% loglog(tsc15,abs(sc15))
% loglog(t,4*10.^9*t.^-7)
% 
% axis([20,tsc400(end),10^-12,1])
% % annotation('textbox','String',[.0,.0,.3,.3],'String',strcat('Amp=',num2str(Amp),...
% %    ' ThicknessPara=',num2str(ThicknessPara),' sourcecentre=',num2str(sourcecentre)),'FitBoxToText','on')
% title('Source location variation , 100 test station, standard Gaussian')
% 
% legend('Source centre = 400', 'Source centre = 120', 'Source centre = 80',...
%     'Source centre = 15', 't.^-7')

% figure(2)
% for i=1:500/dt
% dif1(i)=(ts200(i+100/dt)-ts100(i))/ts200(i+100/dt);
% dif2(i)=(ts300(i+200/dt)-ts200(i+100/dt))/ts300(i+200/dt);
% dif3(i)=(ts400(i+300/dt)-ts300(i+200/dt))/ts400(i+300/dt);
% dif4(i)=(ts800(i+700/dt)-ts400(i+300/dt))/ts800(i+700/dt);
% 
% 
% end
% semilogy((t1(1:end-300/dt-1)),abs((dif1)))
% hold on
% semilogy(t1(1:end-300/dt-1),abs(dif2))
% semilogy(t1(1:end-300/dt-1),abs(dif3))
% semilogy(t1(1:end-300/dt-1),abs(dif4))
% title('modulus Relative difference for varying test station, source location: 80')
% legend('(ts200-ts100)/ts200','(ts300-ts200)/ts300','(ts400-ts300)/ts400','(ts800-ts400)/ts800')
% axis([0,t(end)+1,10^-15,10^-2])

% title('Long duration tail test, log-lin')
% figure(2)
% loglog(t,abs(vectint))
% axis([300,endtime,10^-20,1])
% hold on
% plot(10^-1*(t.^-14))
% plot(10^-3*(t.^-7))
% legend('simulation data','t.^-14','t.^-7')
% 
% grad=(log10(abs(vectint(900/dt)))-log10(abs(vectint(700/dt))))...
%     /(log10(t(900/dt))-log10(t(700/dt)))
% title('VRW__900_1400_400_80_1_1')

% plot(log10(t),log10(abs(vectint)))
% title('Long duration tail test, log-log')
%  axis([2.4,log10(endtime),-20,1])
%  hold on
% %  plot(log10(t),-7*log10(t)+11)
%  hold on
% %figure(2)
% load('VRW_900_800_100_80_1_1.mat')
% %syntax:Potential, domain size, endtime, test station, source location, thickness,
% %amplitude
% 
% loglog(t,abs(vectint))
% axis([0,endtime,10^-20,1])
% annotation('textbox','String',[.0,.0,.3,.3],'String',strcat('Amp=',num2str(Amp),...
%     ' ThicknessPara=',num2str(ThicknessPara),' sourcecentre=',num2str(sourcecentre)),'FitBoxToText','on')
% title(strcat(char(ha),'100 test station 0.2Bw'))
% % 
% % 
% 
% load('VRW_300_400_100_80_1_1.mat')
% %syntax:Potential, domain size, endtime, test station, source location, thickness,
% %amplitude
% hold on
% %figure(3)
% semilogy(t,abs(vectint),'-.k')
% axis([0,endtime,10^-20,1])
% annotation('textbox','String',[.0,.0,.3,.3],'String',strcat('Amp=',num2str(Amp),...
%    ' ThicknessPara=',num2str(ThicknessPara),' sourcecentre=',num2str(sourcecentre)),'FitBoxToText','on')
% title('100 test station skinny bumps 0.4, 0.4w, no 2nd peak')
% % 
%  legend('Additional amp peak 2: 0.4','Additional amp peak 2: 0.4w','No seccond peak')