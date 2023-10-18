clear all, close all, clc

%% generating x-axis

x = [8, 10, 12, 15, 18, 20, 25];
 
%% 8

density_8 = readmatrix("density_8.csv");
density_8(:,6) = [];
chi_8 = density_8(:,5);
chi(1,1) = mean(chi_8);

%% 10

density_10 = readmatrix("density_10.csv");
density_10(:,6) = [];
chi_10 = density_10(:,5);
chi(2,1) = mean(chi_10);

%% 12

density_12 = readmatrix("density_12.csv");
density_12(:,6) = [];
chi_12 = density_12(:,5);
chi(3,1) = mean(chi_12);

%% 15

density_15 = readmatrix("density_15.csv");
density_15(:,6) = [];
chi_15 = density_15(:,5);
chi(4,1) = mean(chi_15);

%% 18

density_18 = readmatrix("density_18.csv");
density_18(:,6) = [];
chi_18 = density_18(:,5);
chi(5,1) = mean(chi_18);

%% 20

density_20 = readmatrix("density_20.csv");
density_20(:,6) = [];
chi_20 = density_20(:,5);
chi_20(8:22) = [];
chi(6,1) = mean(chi_20);

%% 25

density_25 = readmatrix("density_25.csv");
density_25(:,6) = [];
chi_25 = density_25(:,5);
chi_25(14:20) = [];
chi(7,1) = mean(chi_25);

%% Fitting a linear line through the origin
X = x';
Y = chi;
slope = X\Y; %Find the linear regression

% Calculate endpoints of the line
x_endpoints = [0, 25];
y_endpoints = slope * x_endpoints;

% Plotting
figure("Name", "Chi")
plot(x, chi, "o")
hold on
plot(x_endpoints, y_endpoints, 'r--')  % Plotting the linear fit
xlabel("a_{AB} - a_{AA}")
ylabel("\chi-parameter")
xlim([0 25])
ylim([0 6.5])
title("Relation between excess repulsion and effective \chi-parameter")
legend("Data", "Linear Fit")
