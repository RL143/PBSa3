% Load data from CSV file
data = csvread('Energy_time.csv', 1, 0);

% Extract columns
Step = data(:,1);
Time = data(:,2);
Epot = data(:,3);
Ekin = data(:,4);
Etot = data(:,5);

% Plot energies over time
figure;
plot(Time, Epot, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Epot');
hold on;
plot(Time, Ekin, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Ekin');
plot(Time, Etot, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Etot');
xlabel('Time');
ylabel('Energy');
title('Energies over Time');
legend('Location', 'Best');
grid on;
hold off;
