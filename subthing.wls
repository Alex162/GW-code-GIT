r = 2*M*ProductLog[Exp[x/(2*M) - 1]] + 2*M
M = 0.5
Plot[r, {x, -1200, 1500}]
rx = Table[r, {x, -1000, 1000, 0.1}];
Export["/Users/alexkemp/Desktop/Grav wave \
lit/GW-code-GIT/r_of_x.mat", rx]