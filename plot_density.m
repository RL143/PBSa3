clear all, close all, clc

% loading the density data
density = readmatrix("density_3000.csv");

x_coordinate = density(:,1);
data_1 = density(:,2);
data_2 = density(:,3);
data_3 = density(:,4);
data_4 = density(:,5);

figure("Name", "Hopefully a density plot")

plot(x_coordinate, data_1)

hold on
plot(x_coordinate, data_2)
plot(x_coordinate, data_3)
xlabel("x-coordinate")
ylabel("densities")
title("Density profile")
legend("\rho_A", "\rho_B", "\rho_A + \rho_B")

hold off
