w_max = 2200/30*pi;
tau = motor_new.Bw*w_max;
P1= tau*w_max;

i = tau/motor_new.Kt
V = motor_new.Ra*i+motor_new.Kt*w_max;
P2 = i*V