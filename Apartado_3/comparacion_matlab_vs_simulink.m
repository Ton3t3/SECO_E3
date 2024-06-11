clf
ficheros = ["params/beta_0.7_buena.mat", "params/beta_0.6_buena.mat", "params/beta_0.5_buena.mat", "params/beta_crit_buena.mat"];
K = 2652.28/23;    % Con reductora
p = 64.986;

% init
load(ficheros(1));

% Señal de control

% Define the parameters
Ts = 5e-3;
Kp = kp;
Ki = kp/tau_i;
Kd1 = kp*tau_d1;
Kd2 = kp*tau_d2;

Ki_dis = kp*Ts/tau_i;
Kd1_dis = kp*tau_d1/Ts;
Kd2_dis = kp*tau_d2/Ts;




%% Datos simulink
% Acceder a los datos exportados desde el Scope
time = out.tout;
time = time(5:end);
time(1)=0;
signals = out.pos_0_7.signals.values;

%% Datos matlab
%COEF_AMORTIGUAMIENTO = 0.7
load("params/beta_0.7_buena.mat")
amortig = 0.7;

num = [p+K*kp*tau_d1 K*kp K*kp/tau_i];
den = [1 p+K*kp*tau_d1 K*kp K*kp/tau_i];
H = tf(num, den);
[y,t] = step(H,time);
figure(5)
% subplot(3,2,5)
signals = signals(5:end);
plot(t,y,t,signals)
yline(stepinfo(H).Peak,"b-.")
xline(stepinfo(H).TransientTime,"m-.")
xline(stepinfo(H).RiseTime, "r-.")
yline(0.98,"g:")
yline(1.02,"g:")
yline(1,"-")
hold off
legend("Modelo matlab continuo", "Modelo simulink discreto Ts=5ms")
xlabel("t(s)")
ylabel("y(t)")
title("Step for \zeta = "+amortig)

[ceros_dpid,poles_dpid,~] = tf2zp(num,den);
% subplot(3,2,6)
figure(2)
zplane(ceros_dpid, poles_dpid);
title("Pole-Zero Plot for \zeta = "+amortig)

save("figuras/simulink_0_7_bueno.mat", "t", "signals")
