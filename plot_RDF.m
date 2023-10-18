clear
close all
M = readmatrix('RDF_Q3_1.csv');
r = M(:,1);
gr = M(:,2);

grfit = fit(r, gr, 'smoothingspline');

plot(grfit, r, gr)
yline(1, '-.r')
title("Radial distribution function (Fig 2 replica)")
xlabel('Distance [r_c]')
ylabel('g(r)')
xlim([0 1])
ylim([0 1.5])