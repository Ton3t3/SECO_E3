
%% Controlador PID
clear, clc, close all

%% Definición de planta
K = 2652.28;
p = 64.986;
planta = tf(K, [1 p 0]);

%% VARIANDO TAU_D, MANTENIENDO TAU_I, KP
% VALORES ESTABLES  
Kp = 0.8;
tau_d = [0.0010 0.0115 0.05 0.5 1 2 5];
tau_i = 0.3;


figure('Name','\tau_d variation PID','NumberTitle','off');
for i = 1:length(tau_d)
    num = Kp*K*tau_d(i) * [1 1/tau_d(i) 1/tau_d(i)/tau_i];
    denom = [1 p+K*Kp*tau_d(i) K*Kp K*Kp/tau_i];
    
    H = tf(num, denom);
    [y,t] = step(H,2);
    subplot(length(tau_d),2,[1,3,5]);
    hold on
    plot(t,y, 'DisplayName', sprintf('kp=%.3f, tau_d=%.2f tau_i=%.2f', Kp, tau_d(i), tau_i));
    [ceros(:,i), poles(:,i),gain(:,i)] = tf2zp(num, denom);

    tr(i) = stepinfo(H).RiseTime;
    ts(i) = stepinfo(H).SettlingTime;
    Mp(i) = stepinfo(H).Peak;
    tp(i) = stepinfo(H).PeakTime;

end

legend show;
yline(1,'--')
yline(1.02,':','Color','m')
yline(0.98,':','Color','m')
hold off     
for i = 1:length(tau_d)
    subplot(length(tau_d),2,i*2);
    zplane(ceros(:,i),poles(:,i));
    title(sprintf("Pole-Zero Plot for tau_d=%.2f", tau_d(i)))
end

% Reg trans experimental
figure

% Subplot 1: Mp vs tau_d
subplot(2, 2, 1);
plot(tau_d, Mp, 'LineWidth', 2);
xlabel('tau_d');
ylabel('Kp(tau_d)');
title('M_p respecto a tau_d');
grid on;

% Subplot 2: tr vs tau_d
subplot(2, 2, 2);
plot(tau_d, tr, 'LineWidth', 2);
xlabel('tau_d');
ylabel('t_r(tau_d)');
title('t_r respecto a tau_d');
grid on;

% Subplot 3: tp vs tau_d
subplot(2, 2, 3);
plot(tau_d, tp, 'LineWidth', 2);
xlabel('tau_d');
ylabel('t_p(tau_d)');
title('t_p respecto a tau_d');
grid on;

% Subplot 4: ts vs tau_d
subplot(2, 2, 4);
plot(tau_d, ts, 'LineWidth', 2);
xlabel('tau_d');
ylabel('t_s(tau_d)');
title('t_s respecto a tau_d');
grid on;

% errPermRampa2 = p.*tau_i/K/Kp;






%% VARIANDO TAU_I, MANTENIENDO TAU_D, KP
% VALORES ESTABLES  
Kp = 0.8;
tau_d = 0.1;
tau_i = [0.0115 0.05 0.5 1 2 5];

figure('Name','\tau_i variation PID','NumberTitle','off');
for i = 1:length(tau_i)
    num = Kp*K*tau_d * [1 1/tau_d 1/tau_d/tau_i(i)];
    denom = [1 p+K*Kp*tau_d K*Kp K*Kp/tau_i(i)];
    
    H = tf(num, denom);
    [y,t] = step(H,2);
    subplot(length(tau_i),2,[1,3,5]);
    hold on
    plot(t,y, 'DisplayName', sprintf('kp=%.3f, tau_d=%.2f tau_i=%.2f', Kp, tau_d, tau_i(i)));
    [ceros(:,i), poles(:,i),gain(:,i)] = tf2zp(num, denom);

    tr(i) = stepinfo(H).RiseTime;
    ts(i) = stepinfo(H).SettlingTime;
    Mp(i) = stepinfo(H).Peak;
    tp(i) = stepinfo(H).PeakTime;

end

legend show;
yline(1,'--')
yline(1.02,':','Color','m')
yline(0.98,':','Color','m')
hold off     
for i = 1:length(tau_i)
    subplot(length(tau_i),2,i*2);
    zplane(ceros(:,i),poles(:,i));
    title(sprintf("Pole-Zero Plot for tau_i=%.2f", tau_i(i)))
end

% Reg trans experimental
figure

% Subplot 1: Mp vs tau_i
subplot(2, 2, 1);
plot(tau_i, Mp, 'LineWidth', 2);
xlabel('tau_i');
ylabel('Kp(tau_i)');
title('M_p respecto a tau_i');
grid on;

% Subplot 2: tr vs tau_i
subplot(2, 2, 2);
plot(tau_i, tr, 'LineWidth', 2);
xlabel('tau_i');
ylabel('t_r(tau_i)');
title('t_r respecto a tau_i');
grid on;

% Subplot 3: tp vs tau_i
subplot(2, 2, 3);
plot(tau_i, tp, 'LineWidth', 2);
xlabel('tau_i');
ylabel('t_p(tau_i)');
title('t_p respecto a tau_i');
grid on;

% Subplot 4: ts vs tau_d
subplot(2, 2, 4);
plot(tau_i, ts, 'LineWidth', 2);
xlabel('tau_i');
ylabel('t_s(tau_i)');
title('t_s respecto a tau_i');
grid on;


errPermRampa2 = p.*tau_i/K/Kp;


%% VARIANDO KP, MANTENIENDO TAU_D, TAU_I
% VALORES ESTABLES  
Kp = [0.05 0.1 0.5];
tau_d = 0.3;
tau_i = 0.3;

figure('Name','\Kp variation PID','NumberTitle','off');
for i = 1:length(Kp)
    num = Kp(i)*K*tau_d * [1 1/tau_d 1/tau_d/tau_i];
    denom = [1 p+K*Kp(i)*tau_d K*Kp(i) K*Kp(i)/tau_i];
    
    H = tf(num, denom);
    [y,t] = step(H,8);
    subplot(length(Kp),2,[1,3,5]);
    hold on
    plot(t,y, 'DisplayName', sprintf('kp=%.3f, tau_d=%.2f tau_i=%.2f', Kp(i), tau_d, tau_i));
    [ceros(:,i), poles(:,i),gain(:,i)] = tf2zp(num, denom);

    tr(i) = stepinfo(H).RiseTime;
    ts(i) = stepinfo(H).SettlingTime;
    Mp(i) = stepinfo(H).Peak;
    tp(i) = stepinfo(H).PeakTime;

end

legend show;
yline(1,'--')
yline(1.02,':','Color','m')
yline(0.98,':','Color','m')
hold off     
for i = 1:length(Kp)
    subplot(length(Kp),2,i*2);
    zplane(ceros(:,i),poles(:,i));
    title(sprintf("Pole-Zero Plot for Kp=%.2f", Kp(i)))
end

% Reg trans experimental
figure

% Subplot 1: Mp vs Kp
subplot(2, 2, 1);
plot(Kp, Mp, 'LineWidth', 2);
xlabel('Kp');
ylabel('Kp(Kp)');
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

% Subplot 4: ts vs tau_d
subplot(2, 2, 4);
plot(Kp, ts, 'LineWidth', 2);
xlabel('Kp');
ylabel('t_s(Kp)');
title('t_s respecto a Kp');
grid on;

errPermRampa2 = p*tau_i/K./Kp;

%% Controlador PI-D

%% VARIANDO TAU_D, MANTENIENDO TAU_I, KP

Kp = 0.8;
tau_d = [0.0010 0.05 0.1];
tau_i = 0.3;

Mp = zeros(1, length(tau_d)+1);
tp = zeros(1, length(tau_d)+1);
tr = zeros(1, length(tau_d)+1);
ts = zeros(1, length(tau_d)+1);
ceros = zeros(1,length(tau_d)+1); % 1 cero
poles = zeros(3,length(tau_d)+1); % 3 polos
gain = zeros(1,length(tau_d)+1);
figure('Name','\tau_d variation PI-D','NumberTitle','off');
for i = 1:length(tau_d)
    num = K*Kp*[1 1/tau_i];
    denom = [1 p+K*Kp*tau_d(i) K*Kp K*Kp/tau_i];
    
    H = tf(num, denom);
    [y,t] = step(H,2);
    subplot(length(tau_d),2,[1,3,5]);
    hold on
    plot(t,y, 'DisplayName', sprintf('kp=%.3f, tau_d=%.2f, tau_i=%.2f', Kp, tau_d(i), tau_i));
    [ceros(:,i), poles(:,i),gain(:,i)] = tf2zp(num, denom);

    tr(i) = stepinfo(H).RiseTime;
    ts(i) = stepinfo(H).SettlingTime;
    Mp(i) = stepinfo(H).Peak;
    tp(i) = stepinfo(H).PeakTime;

end

legend show;
yline(1,'--')
yline(1.02,':','Color','m')
yline(0.98,':','Color','m')
hold off     
for i = 1:length(tau_d)
    subplot(length(tau_d),2,i*2);
    zplane(ceros(:,i),poles(:,i));
    title(sprintf("Pole-Zero Plot for tau_d=%.2f", tau_d(i)))
end

errPermRampa2 = p*tau_i/K/Kp;

%% VARIANDO TAU_I, MANTENIENDO TAU_D, KP
% VALORES ESTABLES  
Kp = 0.8;
tau_d = 0.3;
tau_i = [0.0115 0.05 0.5];

Mp = zeros(1, length(tau_i)+1);
tp = zeros(1, length(tau_i)+1);
tr = zeros(1, length(tau_i)+1);
ts = zeros(1, length(tau_i)+1);
ceros = zeros(2,length(tau_i)+1); % 2 ceros
poles = zeros(3,length(tau_i)+1); % 3 polos
gain = zeros(1,length(tau_i)+1);
figure('Name','\tau_i variation PI-D','NumberTitle','off');
for i = 1:length(tau_i)
    num = K*Kp*[1 1/tau_i(i)];
    denom = [1 p+K*Kp*tau_d K*Kp K*Kp/tau_i(i)];
    
    H = tf(num, denom);
    [y,t] = step(H,5);
    subplot(length(tau_i),2,[1,3,5]);
    hold on
    plot(t,y, 'DisplayName', sprintf('kp=%.3f, tau_d=%.2f tau_i=%.2f', Kp, tau_d, tau_i(i)));
    [~, poles(:,i),gain(:,i)] = tf2zp(num, denom);

    tr(i) = stepinfo(H).RiseTime;
    ts(i) = stepinfo(H).SettlingTime;
    Mp(i) = stepinfo(H).Peak;
    tp(i) = stepinfo(H).PeakTime;

end

legend show;
yline(1,'--')
yline(1.02,':','Color','m')
yline(0.98,':','Color','m')
hold off     
for i = 1:length(tau_i)
    subplot(length(tau_i),2,i*2);
    zplane(ceros,poles(:,i));
    title(sprintf("Pole-Zero Plot for tau_i=%.2f", tau_i(i)))
end

errPermRampa2 = p.*tau_i/K/Kp;


%% VARIANDO KP, MANTENIENDO TAU_D, TAU_I
% VALORES ESTABLES  
Kp = [0.1 2 10];
tau_d = 0.3;
tau_i = 0.5;

Mp = zeros(1, length(Kp)+1);
tp = zeros(1, length(Kp)+1);
tr = zeros(1, length(Kp)+1);
ts = zeros(1, length(Kp)+1);
ceros = zeros(1,length(Kp)+1); % 1 cero
poles = zeros(3,length(Kp)+1); % 3 polos
gain = zeros(1,length(Kp)+1);
figure('Name','\Kp variation PI-D','NumberTitle','off');
for i = 1:length(Kp)
    num = K*Kp(i)*[1 1/tau_i];
    denom = [1 p+K*Kp(i)*tau_d K*Kp(i) K*Kp(i)/tau_i];
    
    H = tf(num, denom);
    [y,t] = step(H,5);
    subplot(length(Kp),2,[1,3,5]);
    hold on
    plot(t,y, 'DisplayName', sprintf('kp=%.3f, tau_d=%.2f tau_i=%.2f', Kp(i), tau_d, tau_i));
    [ceros(:,i), poles(:,i),gain(:,i)] = tf2zp(num, denom);

    tr(i) = stepinfo(H).RiseTime;
    ts(i) = stepinfo(H).SettlingTime;
    Mp(i) = stepinfo(H).Peak;
    tp(i) = stepinfo(H).PeakTime;

end

legend show;
yline(1,'--')
yline(1.02,':','Color','m')
yline(0.98,':','Color','m')
hold off     
for i = 1:length(Kp)
    subplot(length(Kp),2,i*2);
    zplane(ceros(:,i),poles(:,i));
    title(sprintf("Pole-Zero Plot for Kp=%.2f", Kp(i)))
end

errPermRampa2 = p*tau_i/K./Kp;

