%clc
clear all
close all


x=-4:.25:4;

y=exp(-10*(x.^2));

plot(x,y,'.')

trapz(x,y)
% 
% x=0:0.1:100;
% y=(x-20).^-7;
% semilogy(x,y)
% hold on
% y2=(x-50).^-7;
% semilogy(x-30,y2);