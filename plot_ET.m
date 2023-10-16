clear
close all
M1 = readmatrix('Energy_time.csv');
t1 = M1(:,2);
pot1 = M1(:,3);
kin1 = M1(:,4);
tot1 = M1(:,5);

plot(t1, pot1, t1, kin1, t1, tot1)
legend('E_{pot}', 'E_{kin}', 'E_{tot}')
xlim([0 40])
xlabel('Time [s]')
ylabel('Energy [kT]')
average = mean(kin1(3500:end));
text(20,100, ['Average E_{kin}: ', num2str(average), '  kT'])