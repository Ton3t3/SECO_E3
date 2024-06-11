clear, clc, close all
%%
directory = "realOutputs/test-MOTOR3POS.txt";
name = "test_MOTOR3POS";

load(directory);

eval(sprintf('t1 = %s(:,1);',name));
eval(sprintf('y = %s(:,2);',name));
% ·-·-·-·-·-·-·-·-· %
load("simu.mat")
% ·-·-·-·-·-·-·-·-· %
K = 2652.28/23; % CON REDUCTORA
p = 64.986;
kp=3.9208;
tau_d1 = -0.0169;
tau_d2 = 0.1437;
tau_i = 0.2510;
tau_d = tau_d1+tau_d2;
num = [K*kp*tau_d K*kp K*kp/tau_i];
den = [1 p+K*kp*tau_d1 K*kp K*kp/tau_i];
H = tf(num,den);
[y2,t2] = step(H, 2);
% ·-·-·-·-·-·-·-·-· %
hold on
plot(t1,y, t, signals, t2, y2)
title("Comparación posición modelo real vs Simulink vs Matlab")
legend("Motor Real","Motor Simulink", "Motor Matlab")
xlabel("t(s)")
ylabel("y(t)")

