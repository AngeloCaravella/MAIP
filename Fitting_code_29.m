%% Angelo Caravella%
clear all
close all
clc
%% Dati di input

 x = load('C:\Users\Angelo CARAVELLA\Desktop\Dati\5_rivisto\X.csv');
 y = load('C:\Users\Angelo CARAVELLA\Desktop\Dati\5_rivisto\Y.csv');
 z = load('C:\Users\Angelo CARAVELLA\Desktop\Dati\5_rivisto\Z.csv');
N=length(x);
h_mn = [x';y';z']; %matrice 3*N


%%%%%%%%%%%%%%%%%%%ALGORITMO DI COMPENSAZIONE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  M_best = [0.49;0;0;
%      0.47 ;0;0.49
%     ]; 
% b_best = [1; -1.8; 0.6]; % Vettore zero 3x1
% valori_=[ M_best ;b_best];

phi=8;
while phi<=8
% Funzione obiettivo
fun = @(X_sol) sum(arrayfun(@(n) norm( [X_sol(1), X_sol(2), X_sol(3); X_sol(2), X_sol(4), X_sol(5);
    X_sol(3), X_sol(5), X_sol(6)]*(h_mn(:, n) - X_sol(7:9)) )^2 - 1, 1:N).^2);

%valori_iniziali=[M_best; b_best];
p_rand=rand(6,1);
valori_in= [p_rand(1); 0; 0;
    p_rand(2) ;0;p_rand(3); 
    p_rand(4); -p_rand(5); p_rand(6) ];

X_sol =fmincon(fun, valori_in);

% Estrai Mest e best dai risultati
Mest = [X_sol(1), X_sol(2), X_sol(3); X_sol(2), X_sol(4), X_sol(5);X_sol(3), X_sol(5), X_sol(6) ];
best = X_sol(7:9);

% Calcola h per ogni colonna di h_mn
% Inizializza h come una matrice vuota
h = zeros(3, N);

% Calcola h per ogni colonna di h_mn utilizzando un ciclo for
for n = 1:N
    h(:, n) = Mest *(h_mn(:, n) - best);

end

% Separa la matrice h in 3 vettori colonna
hx = h( 1,:);
hy = h( 2,:);
hz = h( 3,:);

% Calcolo dell'elevation
elevation= (180/pi) * acos(hz);

ele=std(elevation);
mean_e=mean(elevation);
norm_std_e=ele/mean_e;

% Calcolo dell'azimuth (φ)
azimuth= (180/pi) * atan2(hy, hx);
 azi=std(azimuth);
 mean_a=mean(azimuth);
norm_std_a=azi/mean_a;

phi=sqrt(norm_std_a^2+norm_std_e^2);


 end

N=length(hz)

norme_h = arrayfun(@(n) norm(h(:, n)), 1:N);

distanza = abs(norme_h - 1);
sts=std(distanza)


%Crea un grafico tridimensionale
figure;
plot3(hx, hy, hz, 'r'); % 'r' per rosso
hold on; % mantiene il grafico esistente

plot3(x', y', z', 'b'); % 'b' per blu
xlabel('x');
ylabel('y');
zlabel('z');
title('Confronto rosso compensati, blu non compensati');
grid on;
hold off; % rilascia il grafico
figure;
plot3(hx, hy, hz, 'r'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


matrice_cov_dati_originali = cov(h_mn');

matrice_cov = cov(h');

%%%%%%%%%%%%%%PRIMO PLOT %%%%%%%%%%%%%%%%
% Calcolo dell' elevation
elevation= (180/pi) * acos(x);

% Calcolo dell'azimuth (φ)
azimuth= (180/pi) * atan2(y, x);

x= cosd(azimuth) .* cosd(elevation);
y = sind(azimuth) .* cosd(elevation);
z = sind(elevation);

% Plot 3D
figure;
quiver3(zeros(size(x)), zeros(size(y)), zeros(size(z)), real(x), real(y), real(z), 'LineWidth', 1, 'AutoScale','on');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Spostamento del sensore in 3D con dati non compensati');
axis equal; % Imposta gli assi con la stessa scala
grid on;
%%%%%%%%%%%%%%%%%SECONDO PLOT %%%%%%%%%%%%%%%%%%%
elevation= (180/pi) * acos(hz);

% Calcolo dell'azimuth (φ)
azimuth= (180/pi) * atan2(hy, hx);
ele=std(elevation)

azi=std(azimuth)
phi=sqrt(ele^2+azi^2)


x_comp = cosd(azimuth) .* cosd(elevation);
y_comp = sind(azimuth) .* cosd(elevation);
z_comp = sind(elevation);

%Plot 3D
figure;
%Plot della parte reale
quiver3(zeros(size(x_comp)), zeros(size(y_comp)), zeros(size(z_comp)), real(x_comp), real(y_comp), real(z_comp), 'LineWidth', 1, 'Color', 'b', 'AutoScale', 'on');

xlabel('X');
ylabel('Y');
zlabel('Z');
title('Spostamento del sensore in 3D');
axis equal; % Imposta gli assi con la stessa scala
grid on;
legend('Parte reale', 'Parte immaginaria');


