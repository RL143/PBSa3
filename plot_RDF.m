clear
close all
M1 = readmatrix('RDF.csv');
r1 = M1(:,1);
gr1 = M1(:,2);

% M2 = readmatrix('RDF1.csv');
% r2 = M2(:,1);
% gr2 = M2(:,2);

gr1fit = fit(r1, gr1, 'smoothingspline');


plot(gr1fit, r1, gr1)
yline(1, '-.r')
xlabel('Distance [r_c]')
ylabel('g(r)')
% xlim([0 1])
% ylim([0 1.5])
% plot(r1, gr1, r2, gr2)
% xlabel('Distance [r_c]')
% ylabel('g(r)')
% legend('nonconservative', 'conservative')