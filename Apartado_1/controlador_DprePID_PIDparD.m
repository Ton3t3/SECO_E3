%% 
clear, clc, close all
%% 
% Definición de planta
K = 2652.28;
p = 64.986;
planta = tf(K, [1 p 0]);

%=====================%
%% Controlador D|PID %%
%=====================%

%% VARIANDO TAU_D1, MANTENIENDO TAU_I, KP Y TAU_D2
%VALORES ESTABLES  
kp = 0.8;
tau_d2 = p/(K*kp);
tau_i = 0.02;
tau_d1_crit1 = -p/(K*kp);
tau_d1_crit2 = 1/(K*kp)*(1/tau_i-p);
tau_d1_crit = max(tau_d1_crit1,tau_d1_crit2);
tau_d1 = [tau_d1_crit+0.02,tau_d1_crit+0.1,tau_d1_crit+0.5];

Mp_d1 = zeros(1, length(tau_d1)+1);
tp_d1 = zeros(1, length(tau_d1)+1);
tr_d1 = zeros(1, length(tau_d1)+1);
ts_d1 = zeros(1, length(tau_d1)+1);
ceros_d1 = zeros(2,length(tau_d1)+1); % 2 ceros
poles_d1 = zeros(3,length(tau_d1)+1); % 3 polos
gain_d1 = zeros(1,length(tau_d1)+1);

figure('Name','D|PID:tau_d1 variation','NumberTitle','off');
for i = 1:length(tau_d1)

    tau_d = tau_d1(i)+tau_d2;

    num = [K*kp*tau_d K*kp K*kp/tau_i];
    den = [1 p+K*kp*tau_d1(i) K*kp K*kp/tau_i];
    H = tf(num,den);
    %step(H)
    [y,t] = step(H, 2);
    subplot(3,2,[1,3,5]);
    hold on
    plot(t,y)
    [ceros_d1(:,i),poles_d1(:,i),gain_d1(i)] = tf2zp(num,den);

    tr_d1(i) = stepinfo(H).RiseTime;
    ts_d1(i) = stepinfo(H).SettlingTime;
    Mp_d1(i) = stepinfo(H).Peak;
    tp_d1(i) = stepinfo(H).PeakTime;
end
legend("kp="+kp+" \tau_i="+tau_i+" \tau_{d2}="+tau_d2+" \tau_{d1}="+tau_d1(1),"kp="+kp+" \tau_i="+tau_i+" \tau_{d2}="+tau_d2+" \tau_{d1}="+tau_d1(2), "kp="+kp+" \tau_i="+tau_i+" \tau_{d2}="+tau_d2+" \tau_{d1}="+tau_d1(3))
axis([0 2 -0.5 1.5])
yline(1,'--')
yline(1.02,':','Color','m')
yline(0.98,':','Color','m')
hold off     
for i = 1:length(tau_d1)
    subplot(3,2,i*2);
    zplane(ceros_d1(:,i),poles_d1(:,i));
    title("Pole-Zero Plot for \tau_{d1}="+tau_d1(i))
end

%VALOR CRÍTICAMENTE ESTABLE (TAU_D1 = 1/K·kp · (1/tau_i - p)
tau_d1 = tau_d1_crit;
figure('Name','D|PID:tau_d1=1/K·kp · (1/tau_i -p)','NumberTitle','off');

tau_d = tau_d1+tau_d2;
num = [K*kp*tau_d K*kp K*kp/tau_i];
den = [1 p+K*kp*tau_d1 K*kp K*kp/tau_i];
H = tf(num,den);
[y,t] = step(H);
subplot(1,4,[1,2]);
plot(t,y)
[ceros_d1(:,4),poles_d1(:,4),gain_d1(4)] = tf2zp(num,den);

tr_d1(4) = stepinfo(H).RiseTime;
ts_d1(4) = stepinfo(H).SettlingTime;
Mp_d1(4) = stepinfo(H).Peak;
tp_d1(4) = stepinfo(H).PeakTime;

legend("kp="+kp+" \tau_i="+tau_i+" \tau_{d2}="+tau_d2+" \tau_{d1}="+tau_d1)
yline(1, '--')

subplot(1,4,[3,4]);
zplane(ceros_d1(:,4),poles_d1(:,4));
title("Pole-Zero Plot for \tau_{d1}="+tau_d1)


%% VARIANDO TAU_I, MANTENIENDO TAU_D1, KP Y TAU_D2

%VALORES ESTABLES  
kp = 0.8;
% tau_d1 = -p/(K*kp);
tau_d1 = 0.015;
tau_d2 = p/(K*kp);
tau_i_crit = 1/(p+K*kp*tau_d1);
tau_i = [tau_i_crit+0.002, tau_i_crit+0.006, tau_i_crit+0.009];
tau_d = tau_d1+tau_d2;


Mp_i = zeros(1, length(tau_i)+1);
tp_i = zeros(1, length(tau_i)+1);
tr_i = zeros(1, length(tau_i)+1);
ts_i = zeros(1, length(tau_i)+1);
ceros_i = zeros(2,length(tau_i)+1); % 2 ceros
poles_i = zeros(3,length(tau_i)+1); % 3 polos
gain_i = zeros(1,length(tau_i)+1);

figure('Name','D|PID:tau_i variation','NumberTitle','off');
for i = 1:length(tau_i)

    num = [K*kp*tau_d K*kp K*kp/tau_i(i)];
    den = [1 p+K*kp*tau_d1 K*kp K*kp/tau_i(i)];
    H = tf(num,den);
    %step(H)
    [y,t] = step(H);
    subplot(3,2,[1,3,5]);
    hold on
    plot(t,y)
    [ceros_i(:,i),poles_i(:,i),gain_i(i)] = tf2zp(num,den);

    tr_i(i) = stepinfo(H).RiseTime;
    ts_i(i) = stepinfo(H).SettlingTime;
    Mp_i(i) = stepinfo(H).Peak;
    tp_i(i) = stepinfo(H).PeakTime;
end
legend("kp="+kp+" \tau_{d2}="+tau_d2+" \tau_{d1}="+tau_d1+" \tau_i="+tau_i(1),"kp="+kp+" \tau_{d2}="+tau_d2+" \tau_{d1}="+tau_d1+" \tau_i="+tau_i(2), "kp="+kp+" \tau_{d2}="+tau_d2+" \tau_{d1}="+tau_d1+" \tau_i="+tau_i(3))
yline(1,'--')
yline(1.02,':','Color','m')
yline(0.98,':','Color','m')
hold off     
for i = 1:length(tau_i)
    subplot(3,2,i*2);
    zplane(ceros_i(:,i),poles_i(:,i));
    title("Pole-Zero Plot for \tau_i="+tau_i(i))
end

%VALOR CRÍTICAMENTE ESTABLE (TAU_D1 = 1/K·kp · (1/tau_i - p)
tau_i = tau_i_crit;
figure('Name','D|PID:tau_i=1/(K·kp·tau_d1 + p)','NumberTitle','off');

tau_d = tau_d1+tau_d2;
num = [K*kp*tau_d K*kp K*kp/tau_i];
den = [1 p+K*kp*tau_d1 K*kp K*kp/tau_i];
H = tf(num,den);
[y,t] = step(H);
subplot(1,4,[1,2]);
plot(t,y)
[ceros_i(:,4),poles_i(:,4),gain_i(4)] = tf2zp(num,den);

tr_i(4) = stepinfo(H).RiseTime;
ts_i(4) = stepinfo(H).SettlingTime;
Mp_i(4) = stepinfo(H).Peak;
tp_i(4) = stepinfo(H).PeakTime;

legend("kp="+kp+" \tau_{d2}="+tau_d2+" \tau_{d1}="+tau_d1+" \tau_i="+tau_i)
yline(1, '--')

subplot(1,4,[3,4]);
zplane(ceros_i(:,4),poles_i(:,4));
title("Pole-Zero Plot for \tau_i="+tau_i)

%% VARIANDO Kp, MANTENIENDO TAU_D1, TAU_I Y TAU_D2

%VALORES ESTABLES  
kp = [0.3, 0.8, 2];
tau_d1 = 0.09;
tau_i = 0.02;

tau_d2 = zeros(1,length(kp));
Mp_kp = zeros(1, length(kp));
tp_kp = zeros(1, length(kp));
tr_kp = zeros(1, length(kp));
ts_kp = zeros(1, length(kp));
ceros_kp = zeros(2,length(kp)); % 2 ceros
poles_kp = zeros(3,length(kp)); % 3 polos
gain_kp = zeros(1,length(kp));

figure('Name','D|PID:kp variation','NumberTitle','off');
for i = 1:length(kp)

    tau_d2(i) = p/(K*kp(i));
    tau_d = tau_d1+tau_d2(i);

    num = [K*kp(i)*tau_d K*kp(i) K*kp(i)/tau_i];
    den = [1 p+K*kp(i)*tau_d1 K*kp(i) K*kp(i)/tau_i];
    H = tf(num,den);
    %step(H)
    [y,t] = step(H,2);
    subplot(3,2,[1,3,5]);
    hold on
    plot(t,y)
    [ceros_kp(:,i),poles_kp(:,i),gain_kp(i)] = tf2zp(num,den);

    tr_kp(i) = stepinfo(H).RiseTime;
    ts_kp(i) = stepinfo(H).SettlingTime;
    Mp_kp(i) = stepinfo(H).Peak;
    tp_kp(i) = stepinfo(H).PeakTime;
end
legend("\tau_{d1}="+tau_d1+" \tau_i="+tau_i+" \tau_{d2}="+tau_d2(1)+" kp="+kp(1),"\tau_{d1}="+tau_d1+" \tau_i="+tau_i+" \tau_{d2}="+tau_d2(2)+" kp="+kp(2)," \tau_{d1}="+tau_d1+" \tau_i="+tau_i+" \tau_{d2}="+tau_d2(3)+ " kp="+kp(3))
yline(1,'--')
yline(1.02,':','Color','m')
yline(0.98,':','Color','m')
hold off     
for i = 1:length(kp)
    subplot(3,2,i*2);
    zplane(ceros_kp(:,i),poles_kp(:,i));
    title("Pole-Zero Plot for kp="+kp(i))
end


%=====================%
%% Controlador PID-D %%
%=====================%

%% VARIANDO TAU_D1, MANTENIENDO TAU_I, KP Y TAU_D2
%VALORES ESTABLES  
kp = 0.8;
tau_d2 = -p/(K*kp);
tau_i = 0.02;
tau_d1_crit = 1/(K*kp*tau_i);
tau_d1 = [tau_d1_crit+0.01,tau_d1_crit+0.06,tau_d1_crit+0.12];

Mp_d1P = zeros(1, length(tau_d1)+1);
tp_d1P = zeros(1, length(tau_d1)+1);
tr_d1P = zeros(1, length(tau_d1)+1);
ts_d1P = zeros(1, length(tau_d1)+1);
ceros_d1P = zeros(2,length(tau_d1)+1); % 2 ceros
poles_d1P = zeros(3,length(tau_d1)+1); % 3 polos
gain_d1P = zeros(1,length(tau_d1)+1);

figure('Name','PID-D:tau_d1 variation','NumberTitle','off');
for i = 1:length(tau_d1)

    tau_d = tau_d1(i)+tau_d2;

    num = [K*kp*tau_d1(i) K*kp K*kp/tau_i];
    den = [1 p+K*kp*tau_d K*kp K*kp/tau_i];
    H = tf(num,den);
    %step(H)
    [y,t] = step(H);
    subplot(3,2,[1,3,5]);
    hold on
    plot(t,y)
    [ceros_d1P(:,i),poles_d1P(:,i),gain_d1P(i)] = tf2zp(num,den);

    tr_d1P(i) = stepinfo(H).RiseTime;
    ts_d1P(i) = stepinfo(H).SettlingTime;
    Mp_d1P(i) = stepinfo(H).Peak;
    tp_d1P(i) = stepinfo(H).PeakTime;
end
legend("kp="+kp+" \tau_i="+tau_i+" \tau_{d2}="+tau_d2+" \tau_{d1}="+tau_d1(1),"kp="+kp+" \tau_i="+tau_i+" \tau_{d2}="+tau_d2+" \tau_{d1}="+tau_d1(2), "kp="+kp+" \tau_i="+tau_i+" \tau_{d2}="+tau_d2+" \tau_{d1}="+tau_d1(3))
axis([0 3 -0.5 2])
yline(1,'--')
yline(1.02,':','Color','m')
yline(0.98,':','Color','m')
hold off     
for i = 1:length(tau_d1)
    subplot(3,2,i*2);
    zplane(ceros_d1P(:,i),poles_d1P(:,i));
    title("Pole-Zero Plot for \tau_{d1}="+tau_d1(i))
end

%VALOR CRÍTICAMENTE ESTABLE (TAU_D1 = 1/K·kp · (1/tau_i - p)
tau_d1 = tau_d1_crit;
figure('Name','PID-D:tau_d1=1/K·kp·tau_i','NumberTitle','off');

tau_d = tau_d1+tau_d2;
num = [K*kp*tau_d1 K*kp K*kp/tau_i];
den = [1 p+K*kp*tau_d K*kp K*kp/tau_i];
H = tf(num,den);
[y,t] = step(H);
subplot(1,4,[1,2]);
plot(t,y)
[ceros_d1P(:,4),poles_d1P(:,4),gain_d1P(4)] = tf2zp(num,den);

tr_d1P(4) = stepinfo(H).RiseTime;
ts_d1P(4) = stepinfo(H).SettlingTime;
Mp_d1P(4) = stepinfo(H).Peak;
tp_d1P(4) = stepinfo(H).PeakTime;

legend("kp="+kp+" \tau_i="+tau_i+" \tau_{d2}="+tau_d2+" \tau_{d1}="+tau_d1)
yline(1, '--')

subplot(1,4,[3,4]);
zplane(ceros_d1P(:,4),poles_d1P(:,4));
title("Pole-Zero Plot for \tau_{d1}="+tau_d1)

%% VARIANDO TAU_I, MANTENIENDO TAU_D1, KP Y TAU_D2

%VALORES ESTABLES  
kp = 0.8;
tau_d1 = 0.015;
tau_d2 = -p/(K*kp);
tau_i_crit = 1/(K*kp*tau_d1);
tau_i = [tau_i_crit+0.002, tau_i_crit+0.006, tau_i_crit+0.009];
tau_d = tau_d1+tau_d2;


Mp_iP = zeros(1, length(tau_i)+1);
tp_iP = zeros(1, length(tau_i)+1);
tr_iP = zeros(1, length(tau_i)+1);
ts_iP = zeros(1, length(tau_i)+1);
ceros_iP = zeros(2,length(tau_i)+1); % 2 ceros
poles_iP = zeros(3,length(tau_i)+1); % 3 polos
gain_iP = zeros(1,length(tau_i)+1);

figure('Name','PID-D:tau_i variation','NumberTitle','off');
for i = 1:length(tau_i)

    num = [K*kp*tau_d1 K*kp K*kp/tau_i(i)];
    den = [1 p+K*kp*tau_d K*kp K*kp/tau_i(i)];
    H = tf(num,den);
    %step(H)
    [y,t] = step(H);
    subplot(3,2,[1,3,5]);
    hold on
    plot(t,y)
    [ceros_iP(:,i),poles_iP(:,i),gain_iP(i)] = tf2zp(num,den);

    tr_iP(i) = stepinfo(H).RiseTime;
    ts_iP(i) = stepinfo(H).SettlingTime;
    Mp_iP(i) = stepinfo(H).Peak;
    tp_iP(i) = stepinfo(H).PeakTime;
end
legend("kp="+kp+" \tau_{d2}="+tau_d2+" \tau_{d1}="+tau_d1+" \tau_i="+tau_i(1),"kp="+kp+" \tau_{d2}="+tau_d2+" \tau_{d1}="+tau_d1+" \tau_i="+tau_i(2), "kp="+kp+" \tau_{d2}="+tau_d2+" \tau_{d1}="+tau_d1+" \tau_i="+tau_i(3))
yline(1,'--')
yline(1.02,':','Color','m')
yline(0.98,':','Color','m')
hold off     
for i = 1:length(tau_i)
    subplot(3,2,i*2);
    zplane(ceros_iP(:,i),poles_iP(:,i));
    title("Pole-Zero Plot for \tau_i="+tau_i(i))
end

%VALOR CRÍTICAMENTE ESTABLE (TAU_D1 = 1/K·kp · (1/tau_i - p)
tau_i = tau_i_crit;
figure('Name','PID-D:tau_i=1/(K·kp·tau_d1)','NumberTitle','off');

tau_d = tau_d1+tau_d2;
num = [K*kp*tau_d1 K*kp K*kp/tau_i];
    den = [1 p+K*kp*tau_d K*kp K*kp/tau_i];
H = tf(num,den);
[y,t] = step(H);
subplot(1,4,[1,2]);
plot(t,y)
[ceros_iP(:,4),poles_iP(:,4),gain_iP(4)] = tf2zp(num,den);

tr_iP(4) = stepinfo(H).RiseTime;
ts_iP(4) = stepinfo(H).SettlingTime;
Mp_iP(4) = stepinfo(H).Peak;
tp_iP(4) = stepinfo(H).PeakTime;

legend("kp="+kp+" \tau_{d2}="+tau_d2+" \tau_{d1}="+tau_d1+" \tau_i="+tau_i)
yline(1, '--')

subplot(1,4,[3,4]);
zplane(ceros_iP(:,4),poles_iP(:,4));
title("Pole-Zero Plot for \tau_i="+tau_i)

%% VARIANDO Kp, MANTENIENDO TAU_D1, TAU_I Y TAU_D2

%VALORES ESTABLES  
kp = [0.3, 0.5, 1];
tau_d1 = 0.09;
tau_i = 0.02;

tau_d2 = zeros(1,length(kp));
Mp_kpP = zeros(1, length(kp));
tp_kpP = zeros(1, length(kp));
tr_kpP = zeros(1, length(kp));
ts_kpP = zeros(1, length(kp));
ceros_kpP = zeros(2,length(kp)); % 2 ceros
poles_kpP = zeros(3,length(kp)); % 3 polos
gain_kpP = zeros(1,length(kp));

figure('Name','PID-D:kp variation','NumberTitle','off');
for i = 1:length(kp)

    tau_d2(i) = -p/(K*kp(i));
    tau_d = tau_d1+tau_d2(i);

    num = [K*kp(i)*tau_d1 K*kp(i) K*kp(i)/tau_i];
    den = [1 p+K*kp(i)*tau_d K*kp(i) K*kp(i)/tau_i];
    H = tf(num,den);
    %step(H)
    [y,t] = step(H);
    subplot(3,2,[1,3,5]);
    hold on
    plot(t,y)
    [ceros_kpP(:,i),poles_kpP(:,i),gain_kpP(i)] = tf2zp(num,den);

    tr_kpP(i) = stepinfo(H).RiseTime;
    ts_kpP(i) = stepinfo(H).SettlingTime;
    Mp_kpP(i) = stepinfo(H).Peak;
    tp_kpP(i) = stepinfo(H).PeakTime;
end
legend("\tau_{d1}="+tau_d1+" \tau_i="+tau_i+" \tau_{d2}="+tau_d2(1)+" kp="+kp(1),"\tau_{d1}="+tau_d1+" \tau_i="+tau_i+" \tau_{d2}="+tau_d2(2)+" kp="+kp(2)," \tau_{d1}="+tau_d1+" \tau_i="+tau_i+" \tau_{d2}="+tau_d2(3)+ " kp="+kp(3))
yline(1,'--')
yline(1.02,':','Color','m')
yline(0.98,':','Color','m')
hold off     
for i = 1:length(kp)
    subplot(3,2,i*2);
    zplane(ceros_kpP(:,i),poles_kpP(:,i));
    title("Pole-Zero Plot for kp="+kp(i))
end


%====================%
 %% D|PID vs PID-D %%
%====================%

kp = 0.3;
tau_d1 = 0.09;
tau_i = 0.02;
tau_d2_pidd = -p/(K*kp);
tau_d2_dpid = p/(K*kp);
tau_d_pidd = tau_d1 + tau_d2_pidd;
tau_d_dpid = tau_d1 + tau_d2_dpid;
Mp_controller = zeros(1,2); % 1º PID-D   2º D|PID
tp_controller = zeros(1,2); % 1º PID-D   2º D|PID
tr_controller = zeros(1,2); % 1º PID-D   2º D|PID
ts_controller = zeros(1,2); % 1º PID-D   2º D|PID
figure('Name','PID-D vs D|PID','NumberTitle','off');

num = [K*kp*tau_d1 K*kp K*kp/tau_i];
den = [1 p+K*kp*tau_d_pidd K*kp K*kp/tau_i];
H_pidd = tf(num,den);
[y,t] = step(H_pidd);
subplot(2,2,[1,3])
plot(t,y)
hold on
[ceros_pidd,poles_pidd,gain_pidd] = tf2zp(num,den);

num = [K*kp*tau_d_dpid K*kp K*kp/tau_i];
den = [1 p+K*kp*tau_d1 K*kp K*kp/tau_i];
H_dpid = tf(num, den);
[y,t] = step(H_dpid);
plot(t,y)
[ceros_dpid,poles_dpid,gain_dpid] = tf2zp(num,den);
hold off
title("K_p="+kp+" \tau_{d1}="+tau_d1+" \tau_{d2}= ±"+tau_d2_dpid+" \tau_i="+tau_i)
legend("PID-D Controller","D|PID Controller")
yline(1,'--')
yline(1.02,':','Color','m')
yline(0.98,':','Color','m')

subplot(2,2,2)
zplane(ceros_pidd, poles_pidd);
title("Pole-Zero Plot for PID-D")
subplot(2,2,4)
zplane(ceros_dpid, poles_dpid);
title("Pole-Zero Plot for D|PID")

tr_controller(1) = stepinfo(H_pidd).RiseTime;
ts_controller(1) = stepinfo(H_pidd).SettlingTime;
Mp_controller(1) = stepinfo(H_pidd).Peak;
tp_controller(1) = stepinfo(H_pidd).PeakTime;
tr_controller(2) = stepinfo(H_dpid).RiseTime;
ts_controller(2) = stepinfo(H_dpid).SettlingTime;
Mp_controller(2) = stepinfo(H_dpid).Peak;
tp_controller(2) = stepinfo(H_dpid).PeakTime;
