clear, clc, close all
%%
% DATOS PLANTA %
K = 2652.28;
p = 64.986;
% -·-·-·-·-·-· %
K = K/23; % REDUCTORA
% VARIABLES %
amort = [0.5,0.6,0.7];
% -·-·-·-·- %

%% ELECCIÓN DE PARÁMETROS: BETA Y BETA2
%ELECCIÓN DE PARÁMETRO BETA
[beta,Mp1,corteA1, corteB1]=curvas_mp(amort(1), 0.15, 0.07);
[~,Mp2,corteA2, corteB2]=curvas_mp(amort(2), 0.15, 0.07);
[~,Mp3,corteA3, corteB3]=curvas_mp(amort(3), 0.15, 0.07);

figure(1)
hold on
plot(beta, Mp1)
plot(beta, Mp2)
plot(beta, Mp3)
%plot(i1(1,:), i1(2,:), "ro")
%plot(i2(1,:), i2(2,:), "bo")
yline(0.07, "LineStyle","-.", "Color",'g')
yline(0.15, "LineStyle","-.", "Color",'g')
hold off
legend("\zeta = "+amort(1)+" \beta ∈ ("+corteA1(1,:)+", "+corteB1(1,:)+")","\zeta = "+amort(2)+" \beta ∈ ("+corteA2(1,:)+", "+corteB2(1,:)+")","\zeta = "+amort(3)+" \beta ∈ ("+corteA3(1,:)+", "+corteB3(1,:)+")")
xlabel("\beta")
ylabel("M_p")
title("Relación de \beta y M_p")

%ELECCIÓN DE PARÁMETRO BETA_2
[beta2,ts1a, corte1a] = curvas_beta2_ts(amort(1), corteA1(1,:), 0.55*p);
[~,ts1b, corte1b] = curvas_beta2_ts(amort(1), corteB1(1,:), 0.55*p);
[~,ts2a, corte2a] = curvas_beta2_ts(amort(2), corteA2(1,:), 0.55*p);
[~,ts2b, corte2b] = curvas_beta2_ts(amort(2), corteB2(1,:), 0.55*p);
[~,ts3a, corte3a] = curvas_beta2_ts(amort(3), corteA3(1,:), 0.55*p);
[~,ts3b, corte3b] = curvas_beta2_ts(amort(3), corteB3(1,:), 0.55*p);

figure(2)
hold on
plot(beta2, ts1a, beta2, ts1b)
plot(beta2, ts2a, beta2, ts2b)
plot(beta2, ts3a, beta2, ts3b)
% plot(corte1a(1,:), corte1a(2,:), "ro")
% plot(corte1b(1,:), corte1b(2,:), "ro")
yline(0.55*p, "LineStyle","-.", "Color",'g')
hold off
legend("\zeta = "+amort(1)+" => \beta = "+corteA1(1,:)+" => \beta_{2,max} = "+corte1a(1,:),"\zeta = "+amort(1)+" => \beta = "+corteB1(1,:)+" => \beta_{2,max} = "+corte1b(1,:),"\zeta = "+amort(2)+" => \beta = "+corteA2(1,:)+" => \beta_{2,max} = "+corte2a(1,:),...
       "\zeta = "+amort(2)+" => \beta = "+corteB2(1,:)+" => \beta_{2,max} = "+corte2b(1,:),"\zeta = "+amort(3)+" => \beta = "+corteA3(1,:)+" => \beta_{2,max} = "+corte3a(1,:),"\zeta = "+amort(3)+" => \beta = "+corteB3(1,:)+" => \beta_{2,max} = "+corte3b(1,:))
xlabel("\beta_2") 
ylabel("p·t_s")
title("Relación de \beta_2 y t_s")

%COMPROBACIÓN DE TIEMPO DE SUBIDA
[tr1a, bool1a] = comprobador_tr(amort(1), corteA1(1,:), corte1a(1,:), 0.35*p, p);
[tr1b, bool1b] = comprobador_tr(amort(1), corteB1(1,:), corte1b(1,:), 0.35*p, p);
[tr2a, bool2a] = comprobador_tr(amort(2), corteA2(1,:), corte2a(1,:), 0.35*p, p);
[tr2b, bool2b] = comprobador_tr(amort(2), corteB2(1,:), corte2b(1,:), 0.35*p, p);
[tr3a, bool3a] = comprobador_tr(amort(3), corteA3(1,:), corte3a(1,:), 0.35*p, p);
[tr3b, bool3b] = comprobador_tr(amort(3), corteB3(1,:), corte3b(1,:), 0.35*p, p);

figure(3)
Hs = transfer_function(amort(1), corteA1(1,:), corte1a(1,:));
[y1a,t1a] = step(transfer_function(amort(1), corteA1(1,:), corte1a(1,:)),50);
[y1b,t1b] = step(transfer_function(amort(1), corteB1(1,:), corte1b(1,:)),50);
[y2a,t2a] = step(transfer_function(amort(2), corteA2(1,:), corte2a(1,:)),50);
[y2b,t2b] = step(transfer_function(amort(2), corteB2(1,:), corte2b(1,:)),50);
[y3a,t3a] = step(transfer_function(amort(3), corteA3(1,:), corte3a(1,:)),50);
[y3b,t3b] = step(transfer_function(amort(3), corteB3(1,:), corte3b(1,:)),50);
hold on
plot(t1a,y1a,t1b,y1b,t2a,y2a,t2b,y2b,t3a,y3a,t3b,y3b)
yline(0.98, "g-.")
yline(1.02, "g-.")
yline(1,"Linestyle","--","Color","black")
hold off
legend("\zeta = "+amort(1)+" \beta = "+corteA1(1,:)+" \beta_{2,max} = "+corte1a(1,:)+" => t_r = "+tr1a,"\zeta = "+amort(1)+" \beta = "+corteB1(1,:)+" \beta_{2,max} = "+corte1b(1,:)+" => t_r = "+tr1b,"\zeta = "+amort(2)+" \beta = "+corteA2(1,:)+" \beta_{2,max} = "+corte2a(1,:)+" => t_r = "+tr2a,...
       "\zeta = "+amort(2)+" \beta = "+corteB2(1,:)+" \beta_{2,max} = "+corte2b(1,:)+" => t_r = "+tr2b,"\zeta = "+amort(3)+" \beta = "+corteA3(1,:)+" \beta_{2,max} = "+corte3a(1,:)+" => t_r = "+tr3a,"\zeta = "+amort(3)+" \beta = "+corteB3(1,:)+" \beta_{2,max} = "+corte3b(1,:)+" => t_r = "+tr3b)
xlabel("p·t")
ylabel("y(pt)")

%% CÁLCULO DE PARÁMETROS CARACTERÍSTICOS: KP, TAU_D1 Y TAU_I (ASUMIENDO TAU_D2 = P/K*KP)
%SE CALCULARÁ LOS PARÁMETROS DE UNA FUNCIÓN DE TRANSFERENCIA DE LAS
%   CALCULADAS ANTERIORMENTE
beta = corteA1(1,:);
beta_2 = corte1a(1,:);
amortig = amort(1);

kp = (p^2*(2*beta+1/amortig^2))/(beta_2^2*K);
tau_d1 = (beta_2*(beta-beta_2+2))/(p*(2*beta+1/amortig^2));
tau_d2 = p/(K*kp);
tau_i = (beta_2*amortig^2*(2*beta+1/amortig^2))/(beta*p);

num = [p+K*kp*tau_d1 K*kp K*kp/tau_i];
den = [1 p+K*kp*tau_d1 K*kp K*kp/tau_i];
H_dpid = tf(num, den);
[y,t] = step(H_dpid);
figure(4)
subplot(1,2,1)
hold on
plot(t,y)
yline(stepinfo(H_dpid).Peak,"b-.")
xline(stepinfo(H_dpid).SettlingTime,"m-.")
xline(stepinfo(H_dpid).RiseTime, "r-.")
yline(0.98,"g:")
yline(1.02,"g:")
yline(1,"-")
hold off
legend("Kp = "+kp+" \tau_{d1} = "+tau_d1+" \tau_{d2} = "+tau_d2+" \tau_i = "+tau_i,"M_p = "+(stepinfo(H_dpid).Peak-1),"t_s = "+stepinfo(H_dpid).SettlingTime,"t_r = "+stepinfo(H_dpid).RiseTime)
xlabel("t(s)")
ylabel("y(t)")
[ceros_dpid,poles_dpid,~] = tf2zp(num,den);
subplot(1,2,2)
zplane(ceros_dpid, poles_dpid);
title("Pole-Zero Plot for D|PID")

%AL VER QUE FUNCIONA, PROBAMOS A EJECUTAR UNA VERSIÓN CON VALORES DISTINTOS
%A LOS UMBRALES
umbral_Mp = 0.11;
umbral_ts = 0.3;
%COEF_AMORTIGUAMIENTO = 0.5
[~,~,beta, ~]=curvas_mp(amort(1), umbral_Mp, umbral_Mp);
beta = beta(1,:);
[~,~, beta_2] = curvas_beta2_ts(amort(1), beta, umbral_ts*p);
beta_2 = beta_2(1,:);
amortig = amort(1);

kp = (p^2*(2*beta+1/amortig^2))/(beta_2^2*K);
tau_d1 = (beta_2*(beta-beta_2+2))/(p*(2*beta+1/amortig^2));
tau_d2 = p/(K*kp);
tau_i = (beta_2*amortig^2*(2*beta+1/amortig^2))/(beta*p);

num = [p+K*kp*tau_d1 K*kp K*kp/tau_i];
den = [1 p+K*kp*tau_d1 K*kp K*kp/tau_i];
H = tf(num, den);
[y,t] = step(H);
figure(5)
subplot(3,2,1)
hold on
plot(t,y)
yline(stepinfo(H).Peak,"b-.")
xline(stepinfo(H).SettlingTime,"m-.")
xline(stepinfo(H).RiseTime, "r-.")
yline(0.98,"g:")
yline(1.02,"g:")
yline(1,"-")
hold off
legend("Kp = "+kp+" \tau_{d1} = "+tau_d1+" \tau_{d2} = "+tau_d2+" \tau_i = "+tau_i,"M_p = "+(stepinfo(H).Peak-1),"t_s = "+stepinfo(H).SettlingTime,"t_r = "+stepinfo(H).RiseTime)
xlabel("t(s)")
ylabel("y(t)")
title("Step for \zeta = "+amortig)
[ceros_dpid,poles_dpid,~] = tf2zp(num,den);
subplot(3,2,2)
zplane(ceros_dpid, poles_dpid);
title("Pole-Zero Plot for \zeta = "+amortig)

%COEF_AMORTIGUAMIENTO = 0.6
[~,~,beta, ~]=curvas_mp(amort(2), umbral_Mp, umbral_Mp);
beta = beta(1,:);
[~,~, beta_2] = curvas_beta2_ts(amort(2), beta, umbral_ts*p);
beta_2 = beta_2(1,:);
amortig = amort(2);

kp = (p^2*(2*beta+1/amortig^2))/(beta_2^2*K);
tau_d1 = (beta_2*(beta-beta_2+2))/(p*(2*beta+1/amortig^2));
tau_d2 = p/(K*kp);
tau_i = (beta_2*amortig^2*(2*beta+1/amortig^2))/(beta*p);

num = [p+K*kp*tau_d1 K*kp K*kp/tau_i];
den = [1 p+K*kp*tau_d1 K*kp K*kp/tau_i];
H = tf(num, den);
[y,t] = step(H);
subplot(3,2,3)
hold on
plot(t,y)
yline(stepinfo(H).Peak,"b-.")
xline(stepinfo(H).SettlingTime,"m-.")
xline(stepinfo(H).RiseTime, "r-.")
yline(0.98,"g:")
yline(1.02,"g:")
yline(1,"-")
hold off
legend("Kp = "+kp+" \tau_{d1} = "+tau_d1+" \tau_{d2} = "+tau_d2+" \tau_i = "+tau_i,"M_p = "+(stepinfo(H).Peak-1),"t_s = "+stepinfo(H).SettlingTime,"t_r = "+stepinfo(H).RiseTime)
xlabel("t(s)")
ylabel("y(t)")
title("Step for \zeta = "+amortig)
[ceros_dpid,poles_dpid,~] = tf2zp(num,den);
subplot(3,2,4)
zplane(ceros_dpid, poles_dpid);
title("Pole-Zero Plot for \zeta = "+amortig)

%COEF_AMORTIGUAMIENTO = 0.7
[~,~,beta, ~]=curvas_mp(amort(3), umbral_Mp, umbral_Mp);
beta = beta(1,:);
[~,~, beta_2] = curvas_beta2_ts(amort(3), beta, umbral_ts*p);
beta_2 = beta_2(1,:);
amortig = amort(3);

kp = (p^2*(2*beta+1/amortig^2))/(beta_2^2*K);
tau_d1 = (beta_2*(beta-beta_2+2))/(p*(2*beta+1/amortig^2));
tau_d2 = p/(K*kp);
tau_i = (beta_2*amortig^2*(2*beta+1/amortig^2))/(beta*p);

num = [p+K*kp*tau_d1 K*kp K*kp/tau_i];
den = [1 p+K*kp*tau_d1 K*kp K*kp/tau_i];
H = tf(num, den);
[y,t] = step(H);
figure(5)
subplot(3,2,5)
hold on
plot(t,y)
yline(stepinfo(H).Peak,"b-.")
xline(stepinfo(H).SettlingTime,"m-.")
xline(stepinfo(H).RiseTime, "r-.")
yline(0.98,"g:")
yline(1.02,"g:")
yline(1,"-")
hold off
legend("Kp = "+kp+" \tau_{d1} = "+tau_d1+" \tau_{d2} = "+tau_d2+" \tau_i = "+tau_i,"M_p = "+(stepinfo(H).Peak-1),"t_s = "+stepinfo(H).SettlingTime,"t_r = "+stepinfo(H).RiseTime)
xlabel("t(s)")
ylabel("y(t)")
title("Step for \zeta = "+amortig)
[ceros_dpid,poles_dpid,~] = tf2zp(num,den);
subplot(3,2,6)
zplane(ceros_dpid, poles_dpid);
title("Pole-Zero Plot for \zeta = "+amortig)



%% FUNCIONES EMPLEADAS
function[tr, bool] = comprobador_tr(amort, beta, beta_2, margen, p)
    H = transfer_function(amort, beta, beta_2);
    tr = stepinfo(H).RiseTime;
    bool = (tr <= margen);
    tr = tr/p;
end

function[beta_2,ts, corte] = curvas_beta2_ts(amort, beta, margen)
    beta_2 = 0:0.01:60;
    ts = zeros(1,length(beta_2));
    for i = 1:length(beta_2)
        H = transfer_function(amort, beta, beta_2(i));
        ts(i) = stepinfo(H).SettlingTime;
    end
    x = ones(1,length(beta_2));
    x1 = margen.*x;
    corte = InterX([beta_2;ts],[beta_2;x1]);
end


function[beta,Mp, corte_A, corte_B] = curvas_mp(amort, margen_A, margen_B)
    beta = 0:0.01:30;
    Mp = zeros(1,length(beta));       
    beta_2 = 0.25; %beta_2 arbitrario, al no influir en el valor de Mp


    for i = 1:length(beta)
        H = transfer_function(amort, beta(i), beta_2);
        Mp(i) = stepinfo(H).Peak -1;
    end
    x = ones(1,length(beta));
    x1 = margen_A.*x;
    x2 = margen_B.*x;
    corte_A = InterX([beta;Mp],[beta;x1]);
    corte_B = InterX([beta;Mp],[beta;x2]);

end

function[H] = transfer_function(amort, beta, beta_2)
    Q = beta^2-2*beta+1/amort^2;
    r1 = (1/beta_2 * (beta* (1/amort^2-4)+2/amort^2))/Q;
    r2 = (1/(beta_2*amort)^2*(1/amort^2-2*beta))/Q;
    r3 = (beta^3/beta_2)/Q;


    % H = (r1*s + r2)/den1 + r3/den2
    % H1 = (r1*s + r2)/den1
    % H2 = r3/den2

    num1 = [r1 r2];
    den1 = [1 2/beta_2 1/(beta_2*amort)^2];
    H1 = tf(num1,den1);

    num2 = r3;
    den2 = [1 beta/beta_2];
    H2 = tf(num2, den2);

    H = H1 + H2;
end