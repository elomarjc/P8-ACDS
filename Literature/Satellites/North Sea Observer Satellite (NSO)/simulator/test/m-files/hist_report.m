%% Creates the histogram plot used in the report
% need data in z_angle_all

figure(15)
HIST(rad2deg(z_angle_all),20000);
ocur1 = HIST(rad2deg(z_angle_all),20000);
title('Histogram for NSO Robust Control')
ylabel('Occurrences [.]')
xlabel('Error angle [\circ]');
axis([0 0.005 0 max(ocur1)*1.1])

print(figure(15),'-depsc2', 'main_histogram_report.eps')

