close all;
clear all;

K=2;
T1=9;
T2=6;
alfa1=0.15;
alfa2=-0.16;
alfa3=0.25;
alfa4=0.43;
Tsim=70;
u_p_vector = [-0.70, -0.20, 0.50, 0.80];
%% Zadanie 3 - porównanie modeli.
u=1;
u_p = 0;
x1_p = 0;
x2_p = 0;
% Continious model.
sim('zad1', Tsim);

% Discrete model, T = 0.1.
T = 0.1;
sim('zad2', Tsim);
figure(1);
plot(x1_ciagle.time, x1_ciagle.signals.values, x1_dyskr.time,x1_dyskr.signals.values);
xlabel('t[s]');
ylabel('u');
legend('Model ciagly', 'Model dyskretny');
title('Odpowiedz skokowa modelu ciaglego oraz dyskretnego dla T=0.1');

% Discrete model, T = 1.
T = 1;
sim('zad2', Tsim);
figure(2);
plot(x1_ciagle.time, x1_ciagle.signals.values, x1_dyskr.time,x1_dyskr.signals.values);
xlabel('t[s]');
ylabel('u');
legend('Model ciagly', 'Model dyskretny');
title('Odpowiedz skokowa modelu ciaglego oraz dyskretnego dla T=1');

% Discrete model, T = 10.
T = 5;
sim('zad2', Tsim);
figure(3);
plot(x1_ciagle.time, x1_ciagle.signals.values, x1_dyskr.time,x1_dyskr.signals.values);
xlabel('t[s]');
ylabel('u');
legend('Model ciagly', 'Model dyskretny');
title('Odpowiedz skokowa modelu ciaglego oraz dyskretnego dla T=5');
%% Zadanie 4 - charakterystyka statyczna.
y_stat = [];
y_stat_wz = [];
staticCharacteristic = @(u) K*(alfa1*u + alfa2*u^2+alfa3*u^3+alfa4*u^4);

for u=-1:0.1:1
    T=0.1;
    sim('zad2', Tsim);
    y_stat = [y_stat x1_dyskr.signals.values(end)];
    y_stat_wz = [y_stat_wz staticCharacteristic(u)];
end

figure(4);
plot(-1:0.1:1, y_stat);
hold on;
plot(-1:0.1:1, y_stat_wz);
xlabel('u');
ylabel('y_{stat}');
legend('Ch-ka symulacyjna', 'Ch-ka analityczna');
title('Charakterystka statyczna modelu dyskretnego nieliniowego');

%% Zadanie 5 - charakterystyka statyczna zlinearyzowana.
u_p = 0.5;
y_stat_wz = [];
for u=-1:0.1:1
   linearStaticCharacteristic = @(u) K*(T1+T2)/(T1*T2)*(alfa1 + 2*alfa2*u_p+3*alfa3*u_p^2+4*alfa4*u_p^3)*(u-u_p); 
   y_stat_wz = [y_stat_wz linearStaticCharacteristic(u)];
end
figure(5);
plot(-1:0.1:1, y_stat_wz);
xlabel('u');
ylabel('y_stat');
title('Charakterystka statyczna modelu dyskretnego zlinearyzowanego w punkcie u_p=0.5');

%% Zadanie 6 - Zlinearyzowana charakterystyka statyczna na tle charakterystyki dla czterech punktów linearyzacji.
i = 0;
%syms 'K'; 'T1'; 'T2'; 'T'; 'alfa1'; 'alfa2'; 'alfa3'; 'alfa4'; 'u_p';
staticCharacteristic = @(u) K*(alfa1*u + alfa2*u^2+alfa3*u^3+alfa4*u^4);
y_stat_lin = [];
for u=-1:0.1:1
    T=0.1;
    sim('zad2', Tsim);
    y_stat_lin = [y_stat_lin staticCharacteristic(u)];
end

for u_p=u_p_vector
    y_stat_wz = [];
    char_title = "Ch-ki statyczne modeli nieliniowego oraz zlinearyzowanego w punkcie %0.2f";
    %Linearization in different work points.
    linearStaticCharacteristic = @(u) K*(alfa1*(u_p+(u-u_p)) + alfa2*(u_p^2 + 2*u_p*(u-u_p)) +alfa3*(u_p^3 + 3*(u_p^2)*(u-u_p))+alfa4 * (u_p^4 + 4*(u_p^3)*(u-u_p)));
    for u=-1:0.1:1
       y_stat_wz = [y_stat_wz linearStaticCharacteristic(u)];
    end
    figure(6+i);
    plot(-1:0.1:1, y_stat_wz);
    hold on;
    plot(-1:0.1:1, y_stat_lin);
    hold on;
    scatter(u_p, linearStaticCharacteristic(u_p),'g');
    title(sprintf(char_title, u_p));
    legend('Ch-ka liniowa','Ch-ka nieliniowa','Punkt linearyzacji');
    xlabel('u');
    ylabel('y_{stat}');
    i = i + 1;
end

%% Zadanie 9 - porównanie modeli dynamicznych dyskretnych - nieliniowych i zlinearyzowanych.
% Discrete non-linear model, T = 1.
T = 1;
Tsim = 70;
u_steps = [-0.15, -0.10, -0.05, 0.05, 0.1, 0.15];
countX1p = @(u) K*(alfa1*u + alfa2*u^2+alfa3*u^3+alfa4*u^4);
countX2p = @(u) (T1+T2)/(T1*T2)*countX1p(u);
char_title = "Modele dynamiczne nieliniowy oraz zlinearyzowany w punkcie u=%0.2f";

% Discrete linear model.
% Choose linearization point.
i = 0;
for u_p=u_p_vector
    x1_p = countX1p(u_p);
    x2_p = countX2p(u_p);
    
    figure(10+i);
    
    % Add results for all plots.
    for u_i=u_steps
        u = u_p + u_i;
        sim('zad2', Tsim);
        L(1) = plot(x1_dyskr.time,x1_dyskr.signals.values, 'Color',[0.7 0 0]);
        hold on;
    end
    
    %Result for different step value.
    for u_i=u_steps
        u = u_p + u_i;
        sim('zad8',Tsim);
        L(2) = plot(x1_dyskr.time,x1_dyskr.signals.values, 'Color',[0 0.5 0]);
        hold on;
    end
    title(sprintf(char_title, u_p));
    legend(L, {'Model nieliniowy', 'Model liniowy'});
    xlabel('t[s]');
    ylabel('y');
    i = i + 1;
end

%% Zadanie 10 - Transmitancja.
T = 1;
%syms 'z' 'K' 'T1' 'T2' 'T' 'alfa1' 'alfa2' 'alfa3' 'alfa4' 'u_p'
syms 'z';

for u_p=u_p_vector

    A = [(1-T*(T1+T2)/(T1*T2)) T;
        -T/(T1*T2), 1];
    B = [0; (K*T)/(T1*T2)*(alfa1 + alfa2*2*u_p +alfa3*3*u_p^2 + alfa4*4*u_p^3)];
    C = [1 0];
    D = 0;
    Gz = collect(C*(eye(2)*z-A)^-1*B + D,z)
end

%% Zadanie 11 - Wzmocnienie statyczne transmitancji
T = 1;
%syms 'z' 'K' 'T1' 'T2' 'T' 'alfa1' 'alfa2' 'alfa3' 'alfa4' 'u_p'
syms 'z';

K_stat_wz = [];
K_stat = [];
for u_p=-1:0.1:1
    A = [(1-T*(T1+T2)/(T1*T2)) T;
        -T/(T1*T2), 1];
    B = [0; (K*T)/(T1*T2)*(alfa1 + alfa2*2*u_p +alfa3*3*u_p^2 + alfa4*4*u_p^3)];
    C = [1 0];
    D = 0;
    Gz = collect(C*((eye(2)*z-A)^-1)*B + D,z);
        
    K_stat_wz = [K_stat_wz double(subs(Gz, 1))];
    %K_stat_wz = subs(Gz,'z',1)
end
figure(14);
plot(-1:0.1:1, K_stat_wz);
title('Wzmocnienie statyczne transmitancji w roznych punktach pracy');
xlabel('u_p');
ylabel('K');

%% Zadanie dodatkowe - Rozwa¿yæ, ¿e wzmocnienie statyczne transmitancji odpowiada wzmocnieniu statycznego obiektu zlinearyzowanego
T = 1;
syms 'z';
Tsim=1000;
step = 0.1;
x1_p = 0;
x2_p = 0;
K_stat_wz = [];
K_stat = [];
for u_p=u_p_vector(1:3)
    
    %Gain from transfer function.
    A = [(1-T*(T1+T2)/(T1*T2)) T;
        -T/(T1*T2), 1];
    B = [0; (K*T)/(T1*T2)*(alfa1 + alfa2*2*u_p +alfa3*3*u_p^2 + alfa4*4*u_p^3)];
    C = [1 0];
    D = 0;
    Gz = collect(C*((eye(2)*z-A)^-1)*B + D);
    K_stat_wz = [K_stat_wz double(subs(Gz, 1))];
    
    %Gain from simulation.
    u = u_p;
    sim('zad8');
    u = u_p + step;
    temp1 = x1_dyskr.signals.values(end);
    sim('zad8');
    temp2 = x1_dyskr.signals.values(end);
    K_stat = [K_stat (temp2-temp1)/step];
end
figure(15);
scatter(K_stat, K_stat_wz);
hold on;
plot([0:0.1:1],[0:0.1:1]);
legend('wartoœci wzmocnieñ','y=x');
title('Porównanie wzmocnienia statycznego w roznych punktach pracy');
xlabel('K - modelu zlinearyzowanego');
ylabel('K - transmitancji');