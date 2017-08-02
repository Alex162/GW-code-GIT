%clc
clear all
close all


x=-4:0.01:4;

y=10.0003*exp(-100*(x.^2));

plot(x,y,'.')

trapz(x,y)

