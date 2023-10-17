clear
close all
M1 = readmatrix('RDF_Q3_1.csv');
r1 = M1(:,1);
gr1 = M1(:,2);

gr1fit = fit(r1, gr1, 'smoothingspline');

plot(gr1fit, r1, gr1)
yline(1, '-.r')
xlabel('Distance [r_c]')
ylabel('g(r)')
xlim([0 1])
ylim([0 1.5])