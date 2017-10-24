function pp=r_tortoise
tic
M=0.5;
x=-200:0.01:600;
r= @(x) 2*M*lambertw(exp(x./((2*M))-1))+2*M;
rx=r(x);
pp=pchip(x,rx);
save('pp.mat','pp')
toc
end
