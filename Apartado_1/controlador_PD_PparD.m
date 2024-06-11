%% Controlador PD

%% 
clear, clc, clf
%% Controlador PD

% Definición de planta
K = 2652.28;
p = 64.986;
planta = tf(K, [1 p 0]);

% b2 = 0.5;
% amort = [0.3 0.707 1 2];
% wn = p/b2./amort;
% Kp = wn.^2 / sqrt(K);

%% VARIANDO TAU_D, MANTENIENDO KP
% VALORES ESTABLES  
Kp = 2;
tau_d = [0.0010 0.0127911276793297 0.1 0.3 1]; % diferentes valores de tau_d para análisis
Mp = zeros(1, length(tau_d)+1);
tp = zeros(1, length(tau_d)+1);
tr = zeros(1, length(tau_d)+1);
ts = zeros(1, length(tau_d)+1);
ceros = zeros(1,length(tau_d)+1); % 1 cero
poles = zeros(2,length(tau_d)+1); % 2 polos
gain = zeros(1,length(tau_d)+1);
figure('Name','\tau_d variation PD','NumberTitle','off');
for i = 1:length(tau_d)
    num = [Kp*K*tau_d(i) Kp*K];
    denom = [1 p+Kp*K*tau_d(i) Kp*K];
    
    H = tf(num, denom);
    [y,t] = step(H,0.5);
    subplot(3,2,[1,3,5]);
    hold on
    plot(t,y, 'DisplayName', sprintf('kp=%.3f, tau_d=%.2f', Kp, tau_d(i)));
    [ceros(:,i),poles(:,i),gain(:,i)] = tf2zp(num, denom);

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
    zplane(ceros(i),poles(:,i));
    title(sprintf("Pole-Zero Plot for tau_d=%.2f", tau_d(i)))
end

errPermRampa1 = p/K/Kp;

%% VARIANDO Kp, MANTENIENDO TAU_D
% VALORES ESTABLES
Kp = [0.01 0.05 0.5 15];
tau_d = 0.1; % diferentes valores de tau_d para análisis
Mp = zeros(1, length(Kp));
tp = zeros(1, length(Kp)+1);
tr = zeros(1, length(Kp)+1);
ts = zeros(1, length(Kp)+1);
ceros = zeros(1,length(Kp)+1); % 1 cero
poles = zeros(2,length(Kp)+1); % 2 polos
gain = zeros(1,length(Kp)+1);
figure('Name','Kp variation PD','NumberTitle','off');
for i = 1:length(Kp)
    num = [Kp(i)*K*tau_d Kp(i)*K];
    denom = [1 p+Kp(i)*K*tau_d Kp(i)*K];
    
    H = tf(num, denom);
    [y,t] = step(H,8);
    subplot(length(Kp),2,[1,3,5]);
    hold on
    plot(t,y, 'DisplayName', sprintf('kp=%.3f, tau_d=%.2f', Kp(i), tau_d));
    [ceros(:,i),poles(:,i),gain(:,i)] = tf2zp(num, denom);

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
    zplane(ceros(i),poles(:,i));
    title(sprintf("Pole-Zero Plot for Kp=%.2f", Kp(i)))
end

errPermRampa1 = p/K./Kp;


%% Controlador P-D

%% VARIANDO TAU_D, MANTENIENDO KP
% VALORES ESTABLES  
Kp = 2;
tau_d = [0.001 0.0127911276793297 0.0462 0.3]; % diferentes valores de tau_d para análisis
Mp = zeros(1, length(tau_d)+1);
tp = zeros(1, length(tau_d)+1);
tr = zeros(1, length(tau_d)+1);
ts = zeros(1, length(tau_d)+1);
ceros = 0; % 1 cero
poles = zeros(2,length(tau_d)+1); % 2 polos
gain = zeros(1,length(tau_d)+1);
figure('Name','\tau_d variation P-D','NumberTitle','off');
for i = 1:length(tau_d)
    num = Kp*K;
    denom = [1 p+Kp*K*tau_d(i) Kp*K];
    
    H = tf(num, denom);
    [y,t] = step(H,0.5);
    subplot(3,2,[1,3,5]);
    hold on
    plot(t,y, 'DisplayName', sprintf('kp=%.3f, tau_d=%.2f', Kp, tau_d(i)));
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
for i = 1:length(tau_d)
    subplot(length(tau_d),2,i*2);
    zplane(ceros,poles(:,i));
    title(sprintf("Pole-Zero Plot for tau_d=%.2f", tau_d(i)))
end

errPermRampa2 = (p+Kp*K.*tau_d)/K/Kp;


%% VARIANDO KP, MANTENIENDO TAU_D
% VALORES ESTABLES  
Kp = [0.8 5 10];
tau_d = 0.1; % diferentes valores de tau_d para análisis
Mp = zeros(1, length(tau_d)+1);
tp = zeros(1, length(tau_d)+1);
tr = zeros(1, length(tau_d)+1);
ts = zeros(1, length(tau_d)+1);
ceros = 0; % 1 cero
poles = zeros(2,length(tau_d)+1); % 2 polos
gain = zeros(1,length(tau_d)+1);
figure('Name','Kp variation P-D','NumberTitle','off');
for i = 1:length(Kp)
    num = Kp(i)*K;
    denom = [1 p+Kp(i)*K*tau_d Kp(i)*K];
    
    H = tf(num, denom);
    [y,t] = step(H,2);
    subplot(3,2,[1,3,5]);
    hold on
    plot(t,y, 'DisplayName', sprintf('kp=%.3f, tau_d=%.2f', Kp(i), tau_d));
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
for i = 1:length(Kp)
    subplot(length(Kp),2,i*2);
    zplane(ceros,poles(:,i));
    title(sprintf("Pole-Zero Plot for Kp=%.2f", Kp(i)))
end

errPermRampa2 = (p+Kp*K*tau_d)/K./Kp;

