close all
clear

% Load data from CSV file
data = csvread('Energy.csv', 1, 0);

% Extract columns
%Step = data(:,1);
Time = data(:,1);
Epot = data(:,2);
Ekin = data(:,3);
Etot = data(:,4);

% Plot energies over time
figure;
plot(Time, Epot,  'DisplayName', 'Epot');
hold on;
plot(Time, Ekin,  'DisplayName', 'Ekin');
plot(Time, Etot,  'DisplayName', 'Etot');
xlabel('Time [S]');
ylabel('Energy [kT]');
title('Energies over Time');
legend('Location', 'Best');
grid on;
hold off;
