clc
clear

num_bins = 100;
num_part = 400;
steps = 10000*(1-(3/4));
binsize = sqrt(2*(10*1)/1)/num_bins;
M = readmatrix('Histogram.csv');
v = M(2:end,1)*binsize;
data_sim = M(2:end,2)/(num_part*steps*binsize);
MB_data = zeros(num_bins-1,1);

for i = 1:num_bins-1
    MB_data(i) = (1/(2*pi()))^(3/2)*4*pi()*v(i)^2*exp(-v(i)^2 /2);
end

plot(v, data_sim, v, MB_data)
legend("Numerical", "Maxwell-Boltzmann",'interpreter', 'latex')
xlabel('Velocity [v/s]','interpreter', 'latex')
ylabel('Probability/$\Delta$ ', 'interpreter', 'latex')
