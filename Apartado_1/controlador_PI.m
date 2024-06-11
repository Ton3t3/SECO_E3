clear, clc,close all
%% Controlador PI

% Definición de planta
K = 2652.28;
p = 64.986;
planta = tf(K, [1 p 0]);

%% VARIANDO TAU, MANTENIENDO KP
%VALORES ESTABLES  
kp = 0.8;
tau_i =[1/p + 0.02, 1/p + 0.1, 10];
Mp = zeros(1, length(tau_i)+1);
tp = zeros(1, length(tau_i)+1);
tr = zeros(1, length(tau_i)+1);
ts = zeros(1, length(tau_i)+1);
ceros = zeros(1,length(tau_i)+1); % 1 cero
poles = zeros(3,length(tau_i)+1); % 3 polos
gain = zeros(1,length(tau_i)+1);
figure('Name','tau variation','NumberTitle','off');
for i = 1:length(tau_i)
    %contr_pi = tf([kp kp/tau_i(i)], [1 0]);  %No se usa (v)
    
    %H = planta*contr_pi/(1+planta*contr_pi); %No funciona
    num = [K*kp K*kp/tau_i(i)];
    den = [1 p K*kp K*kp/tau_i(i)];
    H = tf(num,den);
    [y,t] = step(H);
    subplot(3,2,[1,3,5]);
    hold on
    plot(t,y)
    [ceros(i),poles(:,i),gain(i)] = tf2zp(num,den);

    % [Mp(i), tp_index] = max(y);
    % Mp(i)=Mp(i)-1;
    % tp(i)=t(tp_index);
    tr(i) = stepinfo(H).RiseTime;
    ts(i) = stepinfo(H).SettlingTime;
    Mp(i) = stepinfo(H).Peak;
    tp(i) = stepinfo(H).PeakTime;
end
title("Efecto variación \tau_i en controlador PI")
xlabel("t(s)")
ylabel("y(t)")
legend("kp=0.8 \tau_i="+tau_i(1),"kp=0.8 \tau_i="+tau_i(2),"kp=0.8 \tau_i="+tau_i(3))
yline(1,'--')
yline(1.02,':','Color','m')
yline(0.98,':','Color','m')
hold off     
for i = 1:length(tau_i)
    subplot(3,2,i*2);
    zplane(ceros(i),poles(:,i));
    title("Pole-Zero Plot for \tau_i="+tau_i(i))
end

%VALOR CRÍTICAMENTE ESTABLE (TAU = 1/P)
tau_i = 1/p;
figure('Name','tau=1/p','NumberTitle','off');
    contr_pi = tf([kp kp/tau_i], [1 0]);
    
%H = planta*contr_pi/(1+planta*contr_pi); %No funciona
num = [K*kp K*kp/tau_i];
den = [1 p K*kp K*kp/tau_i];
H = tf(num,den);
[y,t] = step(H);
subplot(1,4,[1,2]);
plot(t,y)
[ceros(4),poles(:,4),gain(4)] = tf2zp(num,den);

% [Mp_4, tp_index] = max(y);
% Mp(4)=Mp_4-1;
% tp(4)=t(tp_index);
tr(4) = stepinfo(H).RiseTime;
ts(4) = stepinfo(H).SettlingTime;
Mp(4) = stepinfo(H).Peak;
tp(4) = stepinfo(H).PeakTime;

title("Efecto \tau_i=1/p en controlador PI")
xlabel("t(s)")
ylabel("y(t)")
legend("kp=0.8 \tau_i=1/p")
yline(1, '--')

subplot(1,4,[3,4]);
zplane(ceros(4),poles(:,4));
title("Pole-Zero Plot for \tau_i=1/p")



%% FIJANDO TAU, VARIANDO Kp
tau_i = 1/p + 0.1;
kp = [0.2 0.8625 2];
Mp2 = zeros(1, length(kp));
tp2 = zeros(1, length(kp));
tr2 = zeros(1, length(kp));
ts2 = zeros(1, length(kp));
ceros = zeros(1,length(tau_i)); % 1 cero
poles = zeros(3,length(tau_i)); % 3 polos
gain = zeros(1,length(tau_i));
figure('Name','kp variation','NumberTitle','off');
for i = 1:length(kp)
    num = [K*kp(i) K*kp(i)/tau_i];
    den = [1 p K*kp(i) K*kp(i)/tau_i];
    H = tf(num,den);
    [y,t] = step(H, 1.5);

    subplot(3,2,[1,3,5]);
    hold on
    plot(t,y)
    [ceros(i),poles(:,i),gain(i)] = tf2zp(num,den);

    % [Mp2(i), tp_index] = max(y);
    % Mp2(i)=Mp2(i)-1;
    % tp2(i)=t(tp_index);
    tr2(i) = stepinfo(H).RiseTime;
    ts2(i) = stepinfo(H).SettlingTime;
    Mp2(i) = stepinfo(H).Peak;
    tp2(i) = stepinfo(H).PeakTime;

end
title("Efecto variación Kp en controlador PI")
xlabel("t(s)")
ylabel("y(t)")
legend("\tau_i="+tau_i+" kp="+kp(1),"\tau_i="+tau_i+" kp="+kp(2),"\tau_i="+tau_i+" kp="+kp(3))
yline(1,'--')
yline(1.02,':','Color','m')
yline(0.98,':','Color','m')
hold off
for i = 1:length(kp)
    subplot(3,2,i*2);
    zplane(ceros(i),poles(:,i));
    title("Pole-Zero Plot for kp="+kp(i))
end
%% PLOTS PARA MEMORIA
%tau_i

kp = 0.8;
tau_i = 1/p:0.001:1;
Mp = zeros(1,length(tau_i));
tr = zeros(1,length(tau_i));
tp = zeros(1,length(tau_i));
ts = zeros(1,length(tau_i));

for i =1:length(tau_i)
    num = [K*kp K*kp/tau_i(i)];
    den = [1 p K*kp K*kp/tau_i(i)];
    H = tf(num,den);
    Mp(i)=stepinfo(H).Peak-1;
    tr(i)=stepinfo(H).RiseTime;
    tp(i)=stepinfo(H).PeakTime;
    ts(i)=stepinfo(H).SettlingTime;
end

figure(7)
subplot(2,2,1)
plot(tau_i, Mp)
title("Mp respecto a \tau_i")
% axis([1/p 1 1 2])
xlabel("\tau_i")
ylabel("M_p")
grid on
subplot(2,2,2)
plot(tau_i,tr)
title("t_r respecto a \tau_i")
% axis([1/p 1 0.023 0.045])
xlabel("\tau_i")
ylabel("t_r")
grid on
subplot(2,2,3)
plot(tau_i,tp)
title("t_p respecto a \tau_i")
% axis([1/p 1 0.07 0.1])
xlabel("\tau_i")
ylabel("t_p")
grid on
subplot(2,2,4)
plot(tau_i,ts)
title("t_s respecto a \tau_i")
% axis([1/p 1 0 6])
xlabel("\tau_i")
ylabel("t_s")
grid on

%kp
tau_i = 0.11539;
kp = 0.01:0.01:5;
Mp = zeros(1,length(kp));
tr = zeros(1,length(kp));
tp = zeros(1,length(kp));
ts = zeros(1,length(kp));

for i =1:length(kp)
    num = [K*kp(i) K*kp(i)/tau_i];
    den = [1 p K*kp(i) K*kp(i)/tau_i];
    H = tf(num,den);
    Mp(i)=stepinfo(H).Peak-1;
    tr(i)=stepinfo(H).RiseTime;
    tp(i)=stepinfo(H).PeakTime;
    ts(i)=stepinfo(H).SettlingTime;
end

figure(8)
subplot(2,2,1)
plot(kp, Mp)
title("Mp respecto a kp")
% axis([1/p 1 1 2])
xlabel("kp")
ylabel("M_p")
grid on
subplot(2,2,2)
plot(kp,tr)
title("t_r respecto a kp")
% axis([1/p 1 0.023 0.045])
xlabel("kp")
ylabel("t_r")
grid on
subplot(2,2,3)
plot(kp,tp)
title("t_p respecto a kp")
% axis([1/p 1 0.07 0.1])
xlabel("kp")
ylabel("t_p")
grid on
subplot(2,2,4)
plot(kp,ts)
title("t_s respecto a kp")
% axis([1/p 1 0 6])
xlabel("kp")
ylabel("t_s")
grid on

