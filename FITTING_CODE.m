%% Angelo Caravella%
clear all
close all
clc
% Dati di input

 x = load('C:\Users\Engel\Desktop\Dati\6\X.csv');
 y = load('C:\Users\Engel\Desktop\Dati\6\Y.csv');
 z = load('C:\Users\Engel\Desktop\Dati\6\Z.csv');

h_mn = [x';y';z']; % Genera una matrice 3xN di dati casuali

N=length(x);

M_best = [0.49;0;0;0.47 ;0;0.49
    ]; % Matrice identità 3x3
b_best = [1; -1.8; 0.6]; % Vettore zero 3x1

% Funzione obiettivo
fun = @(X_sol) sum(arrayfun(@(n) norm( [X_sol(1), X_sol(2), X_sol(3); X_sol(2), X_sol(4), X_sol(5);
    X_sol(3), X_sol(5), X_sol(6)]*(h_mn(:, n) - X_sol(7:9)) )^2 - 1, 1:N).^2);


X_sol = fminsearch(fun, [M_best; b_best]);


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
h1 = h( 1,:);
h2 = h( 2,:);
h3 = h( 3,:);

%Crea un grafico tridimensionale
figure;
plot3(h1, h2, h3, 'r'); % 'r' per rosso
hold on; % mantiene il grafico esistente

plot3(x', y', z', 'b'); % 'b' per blu
xlabel('x');
ylabel('y');
zlabel('z');
title('Confronto rosso compensati, blu non compensati');
grid on;
hold off; % rilascia il grafico


figure;
plot3(h1, h2, h3, 'g'); 
