clc
clear 
close all

load('VRW_TEST-1200_1500_1200_200_80_1_1')
t200=t;
h200=vectint;
load('VRW_TEST-1200_1500_1400_400_80_1_1')
t400=t;
h400=vectint;

load('VRW_TEST-1200_1500_1800_800_80_1_1')
t800=t;
h800=vectint;

load('VRW_TEST-1200_1500_1200_200_80_1_0p001')
t200p=t;
h200p=vectint;

load('VRW_TEST-1200_1500_1800_800_80_1_0p001')
t800p=t;
h800p=vectint;

figure(2)
semilogy(t200,abs(h200))
hold on
semilogy(t400-200,abs(h400))
semilogy(t800-600,abs(h800))
% plot(t200p,1000*log(abs(h200p)))
% plot(t800p-600,1000*log(abs(h800p)))
axis([250,450,10^-10,1])
ylabel('|h|')
xlabel('t')

legend('source centre=200','source centre=400','source centre=800')
% load('TWscalar_600_100')
% h100=plotything;
% t100=t;
% load('TWscalar_800_200')
% h200=plotything;
% t200=t;
% load('TWscalar_1000_300')
% h300=plotything;
% t300=t;
% load('TWscalar_1200_400')
% h400=plotything;
% t400=t;
% load('TWscalar_1400_500')
% h500=plotything;
% t500=t;
% load('TWscalar_1600_600')
% h600=plotything;
% t600=t;
% % load('TWscalar_1800_700')
% % h700=plotything;
% % t700=t;
% load('TWscalar_2000_800')
% h800=plotything;
% t800=t;


% load('TW_2200_900')
% h900=plotything;
% t900=t;
% load('TW_2400_1000')
% h1000=plotything;
% t1000=t;
% load('TW_3600_1600')
% h1600=plotything;
% t1600=t;
%  plot(t100,log(abs(h100)));
%  hold on
%  plot(t200-200,log(abs(h200)))
%  plot(t300-400,log(abs(h300)))
%  plot(t400-600,log(abs(h400)))
%  plot(t500-800,log(abs(h500)))
%  plot(t600-1000,log(abs(h600)))
% % semilogy(t700-1200,abs(h700))
%  plot(t800-1400,log(abs(h800)))
% % semilogy(t900-1600,abs(h900))
% % semilogy(t1000-1800,abs(h1000))
% % semilogy(t1600-3000,abs(h1600))
%  title('Over lay of varying test station with SCALAR location')
%  legend('100','200','300','400','500','600','800')
%  xlabel('t')
%  ylabel('h')
% axis([200,300,-10,1])

%annotation('textbox','String',[.0,.0,.3,.3],'String',strcat('Amp=',num2str(Amp),...
% title(strcat(char(ha),'100 test station 0.1Bw'))
% hold on
% 
% load('VRW_900_r3800_100_400_1_1')
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