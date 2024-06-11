
%% Definición de planta
% Definir los parámetros y el rango de x
K = 2652.28;
p = 64.986;
Kp = linspace(0.1, 5, 100); % Rango de x y número de puntos

%% Reg. trans
% amort, wn
wn = sqrt(Kp*K);
amort = p ./ (2 .* wn);

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


%% Reg. permanente
errPermRampa = p/K./Kp;

t = 1:1000+1000;
u_trapecio = 1 .* ones(2000,1);
u_trapecio(1:1000) = 1/1000 * t(1:1000);

[y_trian, ~, ~] = lsim(H, u_trapecio, t/1000);

figure(2)

% Respuesta a trapecio
subplot(2,1,1)
plot(t, u_trapecio, t, y_trian);
title('Respuestas del sistema a la entrada trapecio');
xlabel('Tiempo');
ylabel('Amplitud');

% Error con Kp
subplot(2,1,2)
plot(Kp, errPermRampa);
title('Error permanente para una entrada en rampa');
xlabel('Ganancia proporcional (Kp)');
ylabel('Error permanente');





