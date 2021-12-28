%% PROBLEM 2, question 3

clear all; clc; close all;

alpha = 0.35;
z = 0;
q_up = 0.7;

N = 10000;

% in case of success:
investment = (alpha*exp(z)*q_up)^(1/(1-alpha));
profit = exp(z) * investment^(alpha);

u = rand(N,1);
success_events = (u <= q_up);   % boolean vector: 1 if success, 0 in fail

average_profit = mean( profit.*success_events );

%% PROBLEM 2, question 4

clear all; clc; close all;

alpha = 0.35;
z = 0;

n = 10;
q_up = linspace(0.11,0.99,n);

% in case of success:
investment = (alpha*exp(z).*q_up).^(1/(1-alpha));
profit = exp(z).*(investment.^alpha);

N = 10000;
u = rand(N,1);

average_profit = zeros(n,1);
variance_profit = zeros(n,1);

for i=1:n
    success_events = (u <= q_up(i));
    average_profit(i) = mean(profit(i).*success_events);
    variance_profit(i) = var(profit(i).*success_events);
end

figure(1)
plot(q_up, average_profit, '-b', 'LineWidth', 2);
grid on
xlabel('$\bar{q}$','Interpreter', 'latex'); hold on
plot(q_up, variance_profit, '-r', 'LineWidth', 2);
leg = legend('Mean','Variance'); set(leg,'Location','best');
    
    
    
    