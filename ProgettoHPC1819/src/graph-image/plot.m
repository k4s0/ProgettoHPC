figure
load data/serial_out -ascii;
load data/omp_final -ascii;

serial_time = sum(serial_out, 2) / size(serial_out, 2);
omp_times = sum(omp_final, 2) / size(omp_final, 2);

speedup_omp = serial_time ./ omp_times;
p = size(omp_final, 1);

for i=1:p
    efficiency_omp(i) = speedup_omp(i) / i;
end
fig = figure;
plot(speedup_omp), xlabel('num threads'), ylabel('speedup'), xlim([0 p]), xticks(0:1:p), title('OMP SPEEDUP');
saveas(fig,'speedup_omp.png');

plot(efficiency_omp), xlabel('num threads'), ylabel('efficiency'), axis([0 p 0 1]), xticks(0:1:p), title('Strong Scaling Efficiency OMP');
saveas(fig,'efficiency_omp.png');

close all