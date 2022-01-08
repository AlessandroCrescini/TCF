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
    
    % Solve numerically for optimal I
    numeric_sol = vpasolve( exp(z)*(x^alpha) - (B/(q_up-q_down) + 1/q_up)*x + A(i)/q_up ==0, x);
    investment(i) = max(0, numeric_sol);
end


figure(1)
plot(A, investment, '-b', 'LineWidth', 2);
grid on
xlabel('$A$','Interpreter', 'latex'); ylabel('Optimal $I$','Interpreter', 'latex');
ylim([0,0.15]); xlim([0,0.15]);


%% PROBLEM 3, question 5

clear; clc; close all;

syms x

alpha = 0.35;
z = 0;
B = 2;
q_up = 0.7;
q_down = 0.1;

N = 1000;
u = rand(N,1);

% Choose the CEO type
A = rand(N,1)*0.03;            % assets smaller than 0.03
% A = 0.06 + rand(N,1)*0.03;   % assets greater than 0.06

% NB. Since A is uniform, to achieve more statistics and save time one can
% generate uniformly in the interesting intervals (0,0.03) or (0.06,0.09).
% With other distirbutions this doesn't hold in general

success_project = find(u<=q_up);   % useful to reduce the number of iterations in the next loop

profits = zeros(N,1);

for i = 1:length(success_project)
    index = success_project(i);
    numeric_sol = vpasolve( exp(z)*(x^alpha) - (B/(q_up-q_down) + 1/q_up)*x + A(index)/q_up ==0, x);
    investment = max(0, numeric_sol);
    profits(index) = B*investment/(q_up-q_down);
end

disp(mean(profits));


%% PROBLEM 3, question 6

% This section doesn't run particularly fast. 
% Tried an approximation to 3rd degree pol. (instead of using VpaSolve) to accelerate. Here is the 'correct' solution only

clear; clc; close all;

syms x

n = 10;
alpha = 0.35;
z = 0;
B = 2;
q_up = linspace(0.66,0.99,n);
q_down = 0.1;

N = 1000;
u = rand(N,1);

% Choose the CEO type
A_low = rand(N,1)*0.03;            % assets smaller than 0.03
A_high = 0.06 + rand(N,1)*0.03;    % assets greater than 0.06

% Same considerations as in Question 3.5

mean_profit_L = zeros(n,1);
mean_profit_H = zeros(n,1);

% First loop for the values of q_up, nested loop replicates profits calculations

for k = 1:n
    
    success_project = find(u<=q_up(k));   % Boolean vector indentifying successful projects

    profits_L = zeros(N,1);
    profits_H = zeros(N,1);
    
    for i = 1:length(success_project)
        
        index = success_project(i);
        
        numeric_sol = vpasolve( exp(z)*(x^alpha) - (B/(q_up(k)-q_down) + 1/q_up(k))*x + A_low(index)/q_up(k) ==0, x);
        investment = max(0, numeric_sol);
        profits_L(index) = B*investment/(q_up(k)-q_down);
        
        numeric_sol = vpasolve( exp(z)*(x^alpha) - (B/(q_up(k)-q_down) + 1/q_up(k))*x + A_high(index)/q_up(k) ==0, x);
        investment = max(0, numeric_sol);
        profits_H(index) = B*investment/(q_up(k)-q_down);
        
    end
    
    mean_profit_L(k) = mean(profits_L);
    mean_profit_H(k) = mean(profits_H);

end


figure(1)
plot(q_up, mean_profit_L, '--b', 'LineWidth', 2);
grid on
xlabel('$\bar{q}$','Interpreter', 'latex'); 
ylabel('Mean realized profits'); hold on
plot(q_up, mean_profit_H, '-b', 'LineWidth', 2);
leg = legend('A<0.03','A>0.06'); set(leg,'Location','best');

