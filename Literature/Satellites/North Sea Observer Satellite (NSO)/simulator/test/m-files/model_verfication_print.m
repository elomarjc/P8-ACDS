%Plots the output and runs the model verification

%%Prints 0 = off, 1 = on
printEPS = 1;

sim('test/model_verification.mdl');

figure(1);
subplot(3,1,1);
title('Model Error');
plot(model_veri_err.time,model_veri_err.signals(1,1).values)
ylabel('Error quaternion [.]');
axis([0 1600 -1 1]);
grid;
subplot(3,1,2);
plot(model_veri_err.time,model_veri_err.signals(1,2).values)
ylabel('Angular rate [\circ/s]');
axis([0 1600 -4e-16 2e-16]);
grid;
subplot(3,1,3);
plot(model_veri_err.time,model_veri_err.signals(1,3).values)
ylabel('Attitude error [\circ]');
xlabel('Simulation time [s]');
axis([0 1600 0 400]);
grid;

if(printEPS)
print('-depsc','-r300','model_veri_err.eps')
end

figure(2);
title('Model Error');
plot(model_veri_err.time,model_veri_err.signals(1,3).values)
xlabel('Simulation time [s]');
ylabel('Attitude error [\circ]');
axis([300 307 0 0.5]);
grid;

if(printEPS)
print('-depsc','-r300','model_veri_att_err.eps')
end
