clear
close all

M1 = readmatrix('RDF_con.csv');
r1 = M1(:,1);
gr1 = M1(:,2);

M2 = readmatrix('RDF_noncon.csv');
r2 = M2(:,1);
gr2 = M2(:,2);

plot(r1, gr1, r2, gr2)
xlabel('Distance [r_c]')
ylabel('g(r)')
legend('conservative', 'nonconservative')
xlim([0 2.5])
title("Radial distribution conservative and nonconservative forces")