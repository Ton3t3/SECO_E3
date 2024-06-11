clear, clc, close all
%% 
% DATOS PLANTA %
K = 2652.28;
p = 64.986;
% -·-·-·-·-·-· %

kp = 0.3;
tau_d1 = 0.09;
tau_i = 0.02;


%% Seguimiento D|PID (tau_d2 != p/K·kp)

tau_d2 = 0.5;
tau_d = tau_d1 + tau_d2;

num = [K*kp*tau_d K*kp K*kp/tau_i];
den = [1 p+K*kp*tau_d1 K*kp K*kp/tau_i];
H = tf(num, den);
t = 0:0.001:10;
subplot(1,3,1)
%Escalón
u=ones(size(t));
[y_e,~] = lsim(H,u,t);
dist_e=abs(y_e(9000)-u(9000));
plot(t,y_e,t,u,t(9000),u(9000),'mo')
titulo = sprintf("Error previsto = %d\nError = %d",0,dist_e);
title(titulo)
subplot(1,3,2)
%Rampa
u=t;
[y_r,~] =lsim(H,u,t);
dist_r=abs(y_r(9000)-u(9000));
plot(t,y_r,t,u,t(9000), u(9000),'mo')
titulo = sprintf("Error previsto = %d\nError = %d",0,dist_r);
title(titulo)
subplot(1,3,3)
%Parábola
u=t.^2;
[y_p,~] = lsim(H,u,t);
dist_p=abs(y_p(9000)-u(9000));
% dist_calculada = K*kp*tau_d/(p+K*kp*tau_d1);
dist_calculada = abs((p-K*kp*tau_d2)*tau_i/(K*kp));
plot(t,y_p,t,u,t(9000), u(9000),'mo')
titulo = sprintf("Error previsto = %d\nError = %d",dist_calculada,dist_p);
title(titulo)



%% Seguimiento D|PID (tau_d2 = p/K·kp)
figure(2)
tau_d2 = p/(K*kp);
tau_d = tau_d1 + tau_d2;

num = [K*kp*tau_d1 K*kp K*kp/tau_i];
den = [1 p+K*kp*tau_d K*kp K*kp/tau_i];
H = tf(num, den);
t = 0:0.001:10;
subplot(1,3,1)
%Escalón
u=ones(size(t));
[y_e,~] = lsim(H,u,t);
dist_e=abs(y_e(9000)-u(9000));
plot(t,y_e,t,u,t(9000),u(9000),'mo')
titulo = sprintf("Error previsto = %d\nError = %d",0,dist_e);
title(titulo)
subplot(1,3,2)
%Rampa
u=t;
[y_r,~] =lsim(H,u,t);
dist_r=abs(y_r(9000)-u(9000));
plot(t,y_r,t,u,t(9000), u(9000),'mo')
titulo = sprintf("Error previsto = %d\nError = %d",0,dist_r);
title(titulo)
subplot(1,3,3)
%Parábola
u=t.^2;
[y_p,~] = lsim(H,u,t);
dist_p=abs(y_p(9000)-u(9000));
plot(t,y_p,t,u,t(9000), u(9000),'mo')
titulo = sprintf("Error previsto = %d\nError = %d",0,dist_p);
title(titulo)

%% Seguimiento PID-D (tau_d2 != -p/K·kp)
figure(3)

tau_d2 = 0.2;
tau_d = tau_d1 + tau_d2;

num1 = [K*kp*tau_d1 K*kp K*kp/tau_i];
den1 = [1 p+K*kp*tau_d K*kp K*kp/tau_i];
H = tf(num1, den1);
t = 0:0.001:10;
subplot(1,3,1)
%Escalón
u=ones(size(t));
[y_e,~] = lsim(H,u,t);
dist_e=abs(y_e(9000)-u(9000));
plot(t,y_e,t,u,t(9000),u(9000),'mo')
titulo = sprintf("Error previsto = %d\nError = %d",0,dist_e);
title(titulo)
subplot(1,3,2)
%Rampa
u=t;
[y_r,~] =lsim(H,u,t);
dist_r=abs(y_r(9000)-u(9000));
plot(t,y_r,t,u,t(9000), u(9000),'mo')
titulo = sprintf("Error previsto = %d\nError = %d",0,dist_r);
title(titulo)
subplot(1,3,3)
%Parábola
u=t.^2;
[y_p,~] = lsim(H,u,t);
dist_p=abs(y_p(9000)-u(9000));
% dist_calculada = K*kp*tau_d/(p+K*kp*tau_d1);
dist_calculada = abs((p+K*kp*tau_d2)*tau_i/(K*kp));
plot(t,y_p,t,u,t(9000), u(9000),'mo')
titulo = sprintf("Error previsto = %d\nError = %d",dist_calculada,dist_p);
title(titulo)



%% Seguimiento D|PID (tau_d2 = -p/K·kp)
figure(4)
tau_d2 = -p/(K*kp);
tau_d = tau_d1 + tau_d2;

num1 = [K*kp*tau_d1 K*kp K*kp/tau_i];
den1 = [1 p+K*kp*tau_d K*kp K*kp/tau_i];
H = tf(num1, den1);
t = 0:0.001:10;
subplot(1,3,1)
%Escalón
u=ones(size(t));
[y_e,~] = lsim(H,u,t);
dist_e=abs(y_e(9000)-u(9000));
plot(t,y_e,t,u,t(9000),u(9000),'mo')
titulo = sprintf("Error previsto = %d\nError = %d",0,dist_e);
title(titulo)
subplot(1,3,2)
%Rampa
u=t;
[y_r,~] =lsim(H,u,t);
dist_r=abs(y_r(9000)-u(9000));
plot(t,y_r,t,u,t(9000), u(9000),'mo')
titulo = sprintf("Error previsto = %d\nError = %d",0,dist_r);
title(titulo)
subplot(1,3,3)
%Parábola
u=t.^2;
[y_p,~] = lsim(H,u,t);
dist_p=abs(y_p(9000)-u(9000));
plot(t,y_p,t,u,t(9000), u(9000),'mo')
titulo = sprintf("Error previsto = %d\nError = %d",0,dist_p);
title(titulo)





