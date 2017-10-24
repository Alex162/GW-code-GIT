clc
clear all
close all
%initial data plotter vonmisner
load('pp.mat')
x=-35:0.01:2000;
d=-100:0.01:200;
rfunct= @(pp,x) ppval(pp,x);
r=rfunct(pp,x);

M=0.5;

R=(r.^(1/2) + sqrt(r-2*M)).^2;

g=(M^3)./(8*R.^3.*(1+M./(2.*R)));


ddrfunct=(M^3*(-2*M*(-3 + M^2)*sqrt(r) + (-3+7*M^2)*r.^(3/2) - 3*M*r.^(5/2)...
    -  3* M * sqrt(-2*M + r) + 3*M^3*sqrt(-2*M + r) + r.*sqrt(-2*M + r)...
    +4*M^2 .*r .*sqrt(-2*M + r) - 3*M.* r.^2 .*sqrt(-2*M + r)))./ ...
    (8*sqrt(1 - (2*M)./r).*(-2*M + r).^(3/2).*(sqrt(r) + sqrt(-2*M + r)).^6 ...
    .*(1 + M*(-M + r + sqrt(r).*sqrt(-2*M + r))).^2);

Q1=2.*r.*(1-2.*M./r).^2 .* (g./(1-2*M./r) - 1./sqrt(1-2.*M./r).*ddrfunct) +6.*r.*g;

psi=sqrt(4*pi/5)*Q1./(1+3*M./(2*r));

misnerinterp=pchip(x,psi);
save('m_terp','misnerinterp')
plot(x,psi,'o')
hold on
plot(d,ppval(misnerinterp,d))