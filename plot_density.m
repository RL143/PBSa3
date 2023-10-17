clear all, close all, clc

density = readmatrix("density.csv");

density(:,6) = [];

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
legend("Da", "DB", "DA + DB")

hold off



density_3000 = readmatrix("density_3000.csv");

density_3000(:,6) = [];

x_coordinate = density_3000(:,1);
data_1_3 = density_3000(:,2);
data_2_3 = density_3000(:,3);
data_3_3 = density_3000(:,4);
data_4_3 = density_3000(:,5);


figure("Name", "3000")
plot(x_coordinate, data_1_3)

hold on
plot(x_coordinate, data_2_3)
plot(x_coordinate, data_3_3)
legend("Da", "DB", "DA + DB")

hold off

