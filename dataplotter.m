clc
clear all
close all


load('VRW_400_400_100_120_1_p001.mat')
%syntax:Potential, possibly mesh change indicator, domain size, endtime, test station, source location, thickness,
%amplitude
figure(1)
semilogy(t,abs(vectint))
axis([0,endtime,10^-20,1])
%annotation('textbox','String',[.0,.0,.3,.3],'String',strcat('Amp=',num2str(Amp),...
 %   ' ThicknessPara=',num2str(ThicknessPara),' sourcecentre=',num2str(sourcecentre)),'FitBoxToText','on')
title(strcat(char(ha),'100 test station'))
hold on
load('VRW_400_400_100_120_10_p001.mat')
%syntax:Potential, domain size, endtime, test station, source location, thickness,
%amplitude

semilogy(t,abs(vectint))
axis([0,endtime,10^-20,1])
%annotation('textbox','String',[.0,.0,.3,.3],'String',strcat('Amp=',num2str(Amp),...
%    ' ThicknessPara=',num2str(ThicknessPara),' sourcecentre=',num2str(sourcecentre)),'FitBoxToText','on')
title(strcat(char(ha),'100 test station'))
% 
% 

load('VRW_Mp5_400_400_100_120_100_p001')
%syntax:Potential, domain size, endtime, test station, source location, thickness,
%amplitude

semilogy(t,abs(vectint))
axis([0,endtime,10^-20,1])
%annotation('textbox','String',[.0,.0,.3,.3],'String',strcat('Amp=',num2str(Amp),...
%    ' ThicknessPara=',num2str(ThicknessPara),' sourcecentre=',num2str(sourcecentre)),'FitBoxToText','on')
title(strcat(char(ha),'100 test station'))

legend('GWF: 1','GWF: 10', 'GWF: 100')