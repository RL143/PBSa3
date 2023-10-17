clear
close all
M = readmatrix('Energy_Q2_3.csv');
t = M(:,1);
E_pot = M(:,2);
E_kin = M(:,3);
E_tot = M(:,4);

plot(t, E_pot, t, E_kin, t, E_tot)
legend('E_{pot}', 'E_{kin}', 'E_{tot}')
xlabel('Time [s]')
ylabel('Energy [kT]')
avg = mean(E_kin(25:end));
text(3,8500, ['Average E_{kin}: ', num2str(avg), '  kT'])