close all
clear

% Load data from CSV file
data = csvread('Energy_Q2_1.csv', 1, 0);

% Extract columns
%Step = data(:,1);
Time = data(:,1);
Epot = data(:,2);
Ekin = data(:,3);
Etot = data(:,4);

% Plot energies over time
figure;
plot(Time, Epot,  'DisplayName', 'Epot');%'r-', 'LineWidth', 1.5,
hold on;
plot(Time, Ekin,  'DisplayName', 'Ekin');%'b-', 'LineWidth', 1.5,
plot(Time, Etot,  'DisplayName', 'Etot');%'g-', 'LineWidth', 1.5,
xlabel('Time [S]');
ylabel('Energy [kT]');
title('Energies over Time');
legend('Location', 'Best');
grid on;
hold off;
