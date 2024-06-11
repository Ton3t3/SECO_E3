
%% Definición de planta
% Definir los parámetros y el rango de x
K = 2652.28;
p = 64.986;
Kp = linspace(0.3, 5, 100); % Rango de x y número de puntos
tau_d = 0.01;

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
errPermRampaPparD = (p+Kp*K*tau_d)./(Kp*K);
errPermRampaPD = p./(Kp*K);

figure(2)
% Error con Kp P-D
subplot(2,1,1)
plot(Kp, errPermRampaPparD);
title('Error permanente de un controlador P-D para una entrada en rampa variando Kp');
xlabel('Ganancia proporcional (Kp)');
ylabel('Error permanente para controlador P-D');

% Error con Kp PD
subplot(2,1,2)
plot(Kp, errPermRampaPD);
title('Error permanente de un controlador PD para una entrada en rampa variando Kp');
xlabel('Ganancia proporcional (Kp)');
ylabel('Error permanente para controlador PD');

%% Reg. trans variando tau_d
% Definiciones
tau_d = linspace(0.001, 0.1, 100); % Rango de x y número de puntos
Kp = 0.8;

% amort, wn
wn = sqrt(Kp*K);
amort = (p+wn.^2.*tau_d)./(2*wn);

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
plot(tau_d, Mp, 'LineWidth', 2);
xlabel('tau_d');
ylabel('M_p(tau_d)');
title('M_p respecto a tau_d');
grid on;

% Subplot 2: tr vs Kp
subplot(2, 2, 2);
plot(tau_d, tr, 'LineWidth', 2);
xlabel('tau_d');
ylabel('t_r(tau_d)');
title('t_r respecto a tau_d');
grid on;

% Subplot 3: tp vs Kp
subplot(2, 2, 3);
plot(tau_d, tp, 'LineWidth', 2);
xlabel('tau_d');
ylabel('t_p(tau_d)');
title('t_p respecto a tau_d');
grid on;

% Subplot 4: ts vs Kp
subplot(2, 2, 4);
plot(tau_d, ts, 'LineWidth', 2);
xlabel('tau_d');
ylabel('t_s(tau_d)');
title('t_s respecto a tau_d');
grid on;


%% Reg. permanente variando tau_d
errPermRampaPparD = (p+Kp*K*tau_d)./(Kp*K);
errPermRampaPD = p./(Kp*K);

figure(4)
% Error con Kp P-D
subplot(2,1,1)
plot(tau_d, errPermRampaPparD);
title('Error permanente de un controlador P-D para una entrada en rampa variando tau_d');
xlabel('tau_d');
ylabel('Error permanente para controlador P-D');

% Error con Kp PD
subplot(2,1,2)
plot(tau_d, errPermRampaPD, "b*");
title('Error permanente de un controlador PD para una entrada en rampa variando tau_d');
xlabel('tau_d');
ylabel('Error permanente para controlador PD');
ylim([-0.05 0.05])
 
