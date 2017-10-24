clc
clear all
close all

M=0.5
L=2
sigma=-3

r=0:0.01:100;

tort= @(r) r+2*M*log(r./(2*M)-1);

x=tort(r);

plot(r,x,'.')