%% Definición de planta
% Definir los parámetros y el rango de x
K = 2652.28;
p = 64.986;
Kp = linspace(0.3, 5, 100); % Rango de x y número de puntos
tau_d2 = 0.01;

%% Reg. trans variando Kp
% amort, wn
wn = sqrt(Kp*K);
amort = (p+wn.^2*tau_d)./(2*wn);

% Auxiliares
wd = wn.*real(sqrt(1-amort.^2));
phi0 = real(acos(amort));
v = 0.02;

% Calcular M(x) para cada valor de x usando operaciones vectorizadas
Mp = (exp(-(amort ./ real(sqrt(1 - amort.^2))) * pi));
tr = (pi-phi0)./wd;
tp = pi./wd;
ts = log(1./(v*real(sqrt(1-amort.^2))))./(amort.*wn);

%% Plots
figure(1)

% Subplot 1: Mp vs Kp
subplot(2, 2, 1);
plot(Kp, Mp, 'LineWidth', 2);
xlabel('Kp');
ylabel('M_p(Kp)');
title('M_p respecto a Kp');
grid on;

% Subplot 2: tr vs Kp
subplot(2, 2, 2);
plot(Kp, tr, 'LineWidth', 2);
xlabel('Kp');
ylabel('t_r(Kp)');
title('t_r respecto a Kp');
grid on;

% Subplot 3: tp vs Kp
subplot(2, 2, 3);
plot(Kp, tp, 'LineWidth', 2);
xlabel('Kp');
ylabel('t_p(Kp)');
title('t_p respecto a Kp');
grid on;

% Subplot 4: ts vs Kp
subplot(2, 2, 4);
plot(Kp, ts, 'LineWidth', 2);
xlabel('Kp');
ylabel('t_s(Kp)');
title('t_s respecto a Kp');
grid on;

%% Reg. permanente varianso Kp
errPermRampaDprePID = tau_i*(p-K*Kp*tau_d2)./(K*Kp);
errPermRampaPIDparD = tau_i*(p+K*Kp*tau_d2)./(K*Kp);

figure(2)
% Error con Kp D|PID
subplot(2,1,1)
plot(Kp, errPermRampaDprePID);
title('Error permanente de un controlador D|PID para una entrada en rampa variando Kp');
xlabel('Ganancia proporcional (Kp)');
ylabel('Error permanente para controlador D|PID');

% Error con Kp PID-D
subplot(2,1,2)
plot(Kp, errPermRampaPIDparD);
title('Error permanente de un controlador PID-D para una entrada en rampa variando Kp');
xlabel('Ganancia proporcional (Kp)');
ylabel('Error permanente para controlador PID-D');

%% Reg. trans variando tau_d2
% Definiciones
tau_d2 = linspace(0.001, 0.1, 100); % Rango de x y número de puntos
Kp = 0.8;

% amort, wn
wn = sqrt(Kp*K);
amort = (p+wn.^2.*tau_d2)./(2*wn);

% Auxiliares
wd = wn.*real(sqrt(1-amort.^2));
phi0 = real(acos(amort));
v = 0.02;

% Calcular M(x) para cada valor de x usando operaciones vectorizadas
Mp = (exp(-(amort ./ real(sqrt(1 - amort.^2))) * pi));
tr = (pi-phi0)./wd;
tp = pi./wd;
ts = log(1./(v*real(sqrt(1-amort.^2))))./(amort.*wn);

%% Plots
figure(3)

% Subplot 1: Mp vs Kp
subplot(2, 2, 1);
plot(tau_d2, Mp, 'LineWidth', 2);
xlabel('tau_d2');
ylabel('M_p(tau_d2)');
title('M_p respecto a tau_d2');
grid on;

% Subplot 2: tr vs Kp
subplot(2, 2, 2);
plot(tau_d2, tr, 'LineWidth', 2);
xlabel('tau_d2');
ylabel('t_r(tau_d2)');
title('t_r respecto a tau_d2');
grid on;

% Subplot 3: tp vs Kp
subplot(2, 2, 3);
plot(tau_d2, tp, 'LineWidth', 2);
xlabel('tau_d2');
ylabel('t_p(tau_d2)');
title('t_p respecto a tau_d2');
grid on;

% Subplot 4: ts vs Kp
subplot(2, 2, 4);
plot(tau_d2, ts, 'LineWidth', 2);
xlabel('tau_d2');
ylabel('t_s(tau_d2)');
title('t_s respecto a tau_d2');
grid on;


%% Reg. permanente variando tau_d2
errPermRampaDprePID = tau_i*(p-K*Kp*tau_d2)./(K*Kp);
errPermRampaPIDparD = tau_i*(p+K*Kp*tau_d2)./(K*Kp);

figure(4)
% Error con Kp P-D
subplot(2,1,1)
plot(tau_d2, errPermRampaDprePID);
title('Error permanente de un controlador D|PID para una entrada en rampa variando tau_d2');
xlabel('tau_d2');
ylabel('Error permanente para controlador D|PID');

% Error con Kp PD
subplot(2,1,2)
plot(tau_d2, errPermRampaPIDparD, "b*");
title('Error permanente de un controlador PID-D para una entrada en rampa variando tau_d2');
xlabel('tau_d2');
ylabel('Error permanente para controlador PID-D');
ylim([-0.05 0.05])
 
