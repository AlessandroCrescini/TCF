%% PROBLEM 3, question 4

clear all; clc; close all;

syms x

n = 99;
A = linspace(0.01,0.1,n);
alpha = 0.35;
z = 0;
B = 2;
q_up = 0.7;
q_down = 0.1;

investment = zeros(n,1);

for i = 1:n
    numeric_sol = vpasolve( exp(z)*(x^alpha) - (B/(q_up-q_down) + 1/q_up)*x + A(i)/q_up ==0, x);
    investment(i) = max(0, numeric_sol);
end


figure(1)
plot(A, investment, '-b', 'LineWidth', 2);
grid on
xlabel('$A$','Interpreter', 'latex'); ylabel('Optimal $I$','Interpreter', 'latex');
ylim([0.09,0.15]);
