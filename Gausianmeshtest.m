%clc
clear all
close all


x=-4:0.005:4;

y=exp(-100*(x.^2));

plot(x,y,'.')

trapz(x,y)

