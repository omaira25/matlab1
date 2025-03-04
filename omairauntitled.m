clc; clear; close all;

% Parámetros
T = 2; % Período
L = T/2; 
n_max = 6; % Número de términos en la serie

% Definir la función original en un período
t = linspace(0, 2*T, 1000);
f_original = zeros(size(t));

for i = 1:length(t)
    if mod(t(i), T) < 1
        f_original(i) = mod(t(i), T);
    else
        f_original(i) = 1;
    end
end

% Cálculo de la serie de Fourier truncada
f_approx = ones(size(t)) * (3/4); % a_0/2

for n = 1:n_max
    a_n = (pi*n*sin(2*pi*n) + cos(pi*n) - 1) / (pi^2 * n^2);
    b_n = (-pi*n*cos(2*pi*n) + sin(pi*n)) / (pi^2 * n^2);
    
    f_approx = f_approx + a_n * cos(n*pi*t/L) + b_n * sin(n*pi*t/L);
end

% Graficar función original y serie truncada
figure;
plot(t, f_original, 'r', 'LineWidth', 2); hold on;
plot(t, f_approx, 'b--', 'LineWidth', 2);
xlabel('t'); ylabel('f(t)');
title('Serie de Fourier truncada en n=6');
legend('Función original', 'Serie de Fourier');
grid on;