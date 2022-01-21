%% Problem Set 4, Exercise 2

clear; clc; close all;

% Initialization
alpha = 0.35;
delta = 0.15;
mu_z = 0;
rho_z = 0.7;
sigma_e = 0.15;
psi0 = 0.01;
psi1 = 0.005;
r = 0.04;

tolerance = 0.001;
n_z = 5;
m = 4;
n_k = 15;

% Discretization of the productivity shock z
[ z, Pi ] = tauchen(mu_z, rho_z, sigma_e, n_z, m);
z = exp(z)';

% Grid for capital k
k_max = (  (z(end)*alpha/(1+r)) / ( 1- (1-delta)/(1+r)  ) )^(1/(1-alpha));
k_min = 0.1*(  (z(1)*alpha/(1+r)) / ( 1- (1-delta)/(1+r)  ) )^(1/(1-alpha));
k = linspace(k_min, k_max, n_k)';


% The value function is represented by a n_k x n_z matrix
% V(i,j) is the element V(k,z) with k(i) and z(j)

V_guess = zeros(n_k,n_z);
V_opt = ones(n_k,n_z);

% Iterate until the algorithm convergence limit goes below the threshold
% For each iteration step, V_opt(k,z) is computed element-by-element

while max(max(abs(V_guess-V_opt))) > tolerance
    
    V_guess = V_opt;
    
    for i = 1:n_k
        for j = 1:n_z
            argument = z(j)*(k(i)^alpha) - (k - (1-delta)*k(i)) - (psi0/(2*k(i)).*(k - (1-delta)*k(i)).^2) ...
                - psi1*k(i).*(k~=(1-delta)*k(i))   + (V_guess*Pi(:,j))./(1+r);
            V_opt(i,j) = max(argument);
        end
    end
    
end


figure(1)
plot(z, V_opt(1,:), '-r', 'LineWidth', 2);hold on
plot(z, V_opt(end,:), '-b', 'LineWidth', 2);hold on
grid on
xlabel('Productivity Shock $z$','Interpreter','latex');
ylabel('Value Function $V(k,z)$','Interpreter','latex');
leg = legend('$k=0.07$','$k=9.32$','Interpreter','latex'); set(leg,'Location','best');

figure(2)
plot(k, V_opt(:,1), '-r', 'LineWidth', 2);hold on
plot(k, V_opt(:,end), '-b', 'LineWidth', 2);hold on
grid on
xlabel('Current Period Capital $k$','Interpreter','latex');
ylabel('Value Function $V(k,z)$','Interpreter','latex');
leg = legend('$z=0.43$','$z=2.32$','Interpreter','latex'); set(leg,'Location','best');
    









