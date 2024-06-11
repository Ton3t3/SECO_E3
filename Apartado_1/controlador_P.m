
%% Controlador P
clear, clc, clf

% ESTANDARIZADO

%% VARIANDO Kps
% VALORES ESTABLES
% Definición sistema de 2o orden
K = 2652.28;
p = 64.986;
Kp=[p^2/4/K 0.8534, 5, 10];

Mp = zeros(1, length(Kp)+1);
tp = zeros(1, length(Kp)+1);
tr = zeros(1, length(Kp)+1);
ts = zeros(1, length(Kp)+1);
% No ceros
poles = zeros(2,length(Kp)+1); % 2 polos
gain = zeros(1,length(Kp)+1);
figure('Name','Kp variation P','NumberTitle','off');
for i = 1:length(Kp)
    num = Kp(i)*K;
    denom = [1 p Kp(i)*K];
    
    H = tf(num, denom);
    [y,t] = step(H,0.20);
    subplot(length(Kp),2,[1,3,5]);
    hold on
    plot(t,y, 'DisplayName', sprintf('kp=%.3f', Kp(i)));
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
    zplane(0,poles(:,i));
    title(sprintf("Pole-Zero Plot for Kp=%.2f", Kp(i)))
end

errPermRampa1 = p/K./Kp;


%% Controlador P 
 
% PRIMERA PRUEBA. CÁLCULOS DE PARÁMETROS ANALÍTICOS

% % Definición sistema de 2o orden
% K = 2652.28;
% p = 64.986;
% planta = tf(K, [1 p 0]);
% Kp=[0.8534, 5, 10];
% 
% % Régimen transitorio
% wn = sqrt(Kp*K);
% amort = p./(2*wn);
% tol = 0.02;
% ts2 = log(1./(tol.*sqrt(1-amort.^2)))./(wn.*amort);
% 
% Mp = zeros(1, length(Kp));
% tp = zeros(1, length(Kp));
% for i=1:length(Kp)
%     prop = tf(Kp(i),1);
%     H = planta*prop/(1+planta*prop);
%     [y,t] = step(H, 0.2);
% 
%     figure(1)
%     plot(t, y)
%     hold on
% 
%     [Mp(i), tp_index] = max(y);
%     Mp(i)=Mp(i)-1;
%     tp(i)=t(tp_index);
% 
% end
% 
% phi0 = acos(amort);
% wd = pi./tp;
% tr = (pi-phi0)./wd;
% 
% % Régimen permanente
% regPerm = 4*amort.^2/p;
% % 
% % t = 1:1000+1000;
% % u_trapecio = 1 .* ones(2000,1);
% % u_trapecio(1:1000) = 1/1000 * t(1:1000);
% % 
% % [y_trian, ~, ~] = lsim(H, u_trapecio, t/1000);
% % plot(u_trapecio)
% % hold on
% % plot(y_trian);
% % hold off
% 
% 
% figure(2)
% prop = tf(-1,1);
% Ha = prop*planta;
% rlocus(Ha)







