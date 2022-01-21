%% Problem Set 4, Exercise 3

clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To get the plot about alpha, delta, ... please run only              %
% this section, then the correspondent loop section.                   %
%                                                                      %
% Each section replies exactly the code for Ex.2, slightly rearranged  %
% to consider different values for alpha, delta, ...                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization
N = 10;

alpha = linspace(0.15,0.35,N);
delta = linspace(0.05,0.15,N);
mu_z = 0;
rho_z = linspace(0.5,0.9,N);
sigma_e = linspace(0.05,0.5,N);
psi0 = linspace(0,0.2,N);
psi1 = linspace(0,0.1,N);
r = linspace(0.01,0.2,N);

tolerance = 0.001;
n_z = 5;
m = 4;
n_k = 15;

%% LOOP on: alpha
V = NaN(N,1);
investment = NaN(N,1);
for index = 1:N
    
    delta = 0.15;
    rho_z = 0.7;
    sigma_e = 0.15;
    psi0 = 0.01;
    r = 0.04;
    psi1 = 0.005;
    
    [ z, Pi ] = tauchen(mu_z, rho_z, sigma_e, n_z, m);
    z = exp(z)';
    k_max = (  (z(end)*alpha(index)/(1+r)) / ( 1- (1-delta)/(1+r)  ) )^(1/(1-alpha(index)));
    k_min = 0.1*(  (z(1)*alpha(index)/(1+r)) / ( 1- (1-delta)/(1+r)  ) )^(1/(1-alpha(index)));
    k = linspace(k_min, k_max, n_k)';

    V_guess = zeros(n_k,n_z);
    V_opt = ones(n_k,n_z);
    
    K_values = NaN(n_k,n_z); 

    while max(max(abs(V_guess-V_opt))) > tolerance
    
        V_guess = V_opt;
    
        for i = 1:n_k
            for j = 1:n_z
                argument = z(j)*(k(i)^alpha(index)) - (k - (1-delta)*k(i)) - (psi0/(2*k(i)).*(k - (1-delta)*k(i)).^2) ...
                    - psi1*k(i).*(k~=(1-delta)*k(i))   + (V_guess*Pi(:,j))./(1+r);
                [V_opt(i,j), maxposition] = max(argument);
                K_values(i,j) = k(maxposition);
            end
        end
    
    end

    V(index) = mean(mean(V_opt));
    investment(index) = mean(mean((K_values - (1-delta).*k)));
end

figure(1)
plot(alpha, V, '-b', 'LineWidth', 2);hold on
grid on
xlabel('Profit Curvature $\alpha$','Interpreter','latex');
ylabel('Average Value Function $V(k,z)$','Interpreter','latex');

figure(2)
plot(alpha, investment, '--b', 'LineWidth', 2);hold on
grid on
xlabel('Profit Curvature $\alpha$','Interpreter','latex');
ylabel('Average Optimal Investment $I$','Interpreter','latex');

%% LOOP on: delta
V = NaN(N,1);
investment = NaN(N,1);
for index = 1:N
    
    alpha = 0.35;
    rho_z = 0.7;
    sigma_e = 0.15;
    psi0 = 0.01;
    psi1 = 0.005;
    r = 0.04;
    
    [ z, Pi ] = tauchen(mu_z, rho_z, sigma_e, n_z, m);
    z = exp(z)';
    k_max = (  (z(end)*alpha/(1+r)) / ( 1- (1-delta(index))/(1+r)  ) )^(1/(1-alpha));
    k_min = 0.1*(  (z(1)*alpha/(1+r)) / ( 1- (1-delta(index))/(1+r)  ) )^(1/(1-alpha));
    k = linspace(k_min, k_max, n_k)';

    V_guess = zeros(n_k,n_z);
    V_opt = ones(n_k,n_z);

    K_values = NaN(n_k,n_z); 

    while max(max(abs(V_guess-V_opt))) > tolerance
    
        V_guess = V_opt;
    
        for i = 1:n_k
            for j = 1:n_z
                argument = z(j)*(k(i)^alpha) - (k - (1-delta(index))*k(i)) - (psi0/(2*k(i)).*(k - (1-delta(index))*k(i)).^2) ...
                    - psi1*k(i).*(k~=(1-delta(index))*k(i))   + (V_guess*Pi(:,j))./(1+r);
                [V_opt(i,j), maxposition] = max(argument);
                K_values(i,j) = k(maxposition);
            end
        end
    
    end

    V(index) = mean(mean(V_opt));
    investment(index) = mean(mean((K_values - (1-delta(index)).*k)));
end

figure(1)
plot(delta, V, '-b', 'LineWidth', 2);hold on
grid on
xlabel('Capital Depreciation Rate $\delta$','Interpreter','latex');
ylabel('Average Value Function $V(k,z)$','Interpreter','latex');

figure(2)
plot(delta, investment, '--b', 'LineWidth', 2);hold on
grid on
xlabel('Capital Depreciation Rate $\delta$','Interpreter','latex');
ylabel('Average Optimal Investment $I$','Interpreter','latex');

%% LOOP on: rho
V = NaN(N,1);
investment = NaN(N,1);
for index = 1:N
    
    alpha = 0.35;
    delta = 0.15;
    sigma_e = 0.15;
    psi0 = 0.01;
    psi1 = 0.005;
    r = 0.04;
    
    [ z, Pi ] = tauchen(mu_z, rho_z(index), sigma_e, n_z, m);
    z = exp(z)';
    k_max = (  (z(end)*alpha/(1+r)) / ( 1- (1-delta)/(1+r)  ) )^(1/(1-alpha));
    k_min = 0.1*(  (z(1)*alpha/(1+r)) / ( 1- (1-delta)/(1+r)  ) )^(1/(1-alpha));
    k = linspace(k_min, k_max, n_k)';

    V_guess = zeros(n_k,n_z);
    V_opt = ones(n_k,n_z);

    K_values = NaN(n_k,n_z); 

    while max(max(abs(V_guess-V_opt))) > tolerance
    
        V_guess = V_opt;
    
        for i = 1:n_k
            for j = 1:n_z
                argument = z(j)*(k(i)^alpha) - (k - (1-delta)*k(i)) - (psi0/(2*k(i)).*(k - (1-delta)*k(i)).^2) ...
                    - psi1*k(i).*(k~=(1-delta)*k(i))   + (V_guess*Pi(:,j))./(1+r);
                [V_opt(i,j), maxposition] = max(argument);
                K_values(i,j) = k(maxposition);
            end
        end
    
    end

    V(index) = mean(mean(V_opt));
    investment(index) = mean(mean((K_values - (1-delta).*k)));
end

figure(1)
plot(rho_z, V, '-b', 'LineWidth', 2);hold on
grid on
xlabel('Auto-Regressive coeff. $\rho_z$','Interpreter','latex');
ylabel('Average Value Function $V(k,z)$','Interpreter','latex');

figure(2)
plot(rho_z, investment, '--b', 'LineWidth', 2);hold on
grid on
xlabel('Auto-Regressive coeff. $\rho_z$','Interpreter','latex');
ylabel('Average Optimal Investment $I$','Interpreter','latex');

%% LOOP on: sigma_e
V = NaN(N,1);
investment = NaN(N,1);
for index = 1:N
    
    alpha = 0.35;
    delta = 0.15;
    rho_z = 0.7;
    psi0 = 0.01;
    r = 0.04;
    psi1 = 0.005;
    
    [ z, Pi ] = tauchen(mu_z, rho_z, sigma_e(index), n_z, m);
    z = exp(z)';
    k_max = (  (z(end)*alpha/(1+r)) / ( 1- (1-delta)/(1+r)  ) )^(1/(1-alpha));
    k_min = 0.1*(  (z(1)*alpha/(1+r)) / ( 1- (1-delta)/(1+r)  ) )^(1/(1-alpha));
    k = linspace(k_min, k_max, n_k)';

    V_guess = zeros(n_k,n_z);
    V_opt = ones(n_k,n_z);

    K_values = NaN(n_k,n_z); 

    while max(max(abs(V_guess-V_opt))) > tolerance
    
        V_guess = V_opt;
    
        for i = 1:n_k
            for j = 1:n_z
                argument = z(j)*(k(i)^alpha) - (k - (1-delta)*k(i)) - (psi0/(2*k(i)).*(k - (1-delta)*k(i)).^2) ...
                    - psi1*k(i).*(k~=(1-delta)*k(i))   + (V_guess*Pi(:,j))./(1+r);
                [V_opt(i,j), maxposition] = max(argument);
                K_values(i,j) = k(maxposition);
            end
        end
    
    end

    V(index) = mean(mean(V_opt));
    investment(index) = mean(mean((K_values - (1-delta).*k)));
end

figure(1)
plot(sigma_e, V, '-b', 'LineWidth', 2);hold on
grid on
xlabel('Stochastic noise st.dev. $\sigma_\epsilon$','Interpreter','latex');
ylabel('Average Value Function $V(k,z)$','Interpreter','latex');

figure(2)
plot(sigma_e, investment, '--b', 'LineWidth', 2);hold on
grid on
xlabel('Stochastic noise st.dev. $\sigma_\epsilon$','Interpreter','latex');
ylabel('Average Optimal Investment $I$','Interpreter','latex');

%% LOOP on: psi_0
V = NaN(N,1);
investment = NaN(N,1);
for index = 1:N
    
    alpha = 0.35;
    delta = 0.15;
    rho_z = 0.7;
    sigma_e = 0.15;
    r = 0.04;
    psi1 = 0.005;
    
    [ z, Pi ] = tauchen(mu_z, rho_z, sigma_e, n_z, m);
    z = exp(z)';
    k_max = (  (z(end)*alpha/(1+r)) / ( 1- (1-delta)/(1+r)  ) )^(1/(1-alpha));
    k_min = 0.1*(  (z(1)*alpha/(1+r)) / ( 1- (1-delta)/(1+r)  ) )^(1/(1-alpha));
    k = linspace(k_min, k_max, n_k)';

    V_guess = zeros(n_k,n_z);
    V_opt = ones(n_k,n_z);

    K_values = NaN(n_k,n_z); 

    while max(max(abs(V_guess-V_opt))) > tolerance
    
        V_guess = V_opt;
    
        for i = 1:n_k
            for j = 1:n_z
                argument = z(j)*(k(i)^alpha) - (k - (1-delta)*k(i)) - (psi0(index)/(2*k(i)).*(k - (1-delta)*k(i)).^2) ...
                    - psi1*k(i).*(k~=(1-delta)*k(i))   + (V_guess*Pi(:,j))./(1+r);
                [V_opt(i,j), maxposition] = max(argument);
                K_values(i,j) = k(maxposition);
            end
        end
    
    end

    V(index) = mean(mean(V_opt));
    investment(index) = mean(mean((K_values - (1-delta).*k)));
end

figure(1)
plot(psi0, V, '-b', 'LineWidth', 2);hold on
grid on
xlabel('Variable Investment Cost $\psi_0$','Interpreter','latex');
ylabel('Average Value Function $V(k,z)$','Interpreter','latex');

figure(2)
plot(psi0, investment, '--b', 'LineWidth', 2);hold on
grid on
xlabel('Variable Investment Cost $\psi_0$','Interpreter','latex');
ylabel('Average Optimal Investment $I$','Interpreter','latex');

%% LOOP on: psi_1
V = NaN(N,1);
investment = NaN(N,1);
for index = 1:N
    
    alpha = 0.35;
    delta = 0.15;
    rho_z = 0.7;
    sigma_e = 0.15;
    r = 0.04;
    psi0 = 0.01;
    
    [ z, Pi ] = tauchen(mu_z, rho_z, sigma_e, n_z, m);
    z = exp(z)';
    k_max = (  (z(end)*alpha/(1+r)) / ( 1- (1-delta)/(1+r)  ) )^(1/(1-alpha));
    k_min = 0.1*(  (z(1)*alpha/(1+r)) / ( 1- (1-delta)/(1+r)  ) )^(1/(1-alpha));
    k = linspace(k_min, k_max, n_k)';

    V_guess = zeros(n_k,n_z);
    V_opt = ones(n_k,n_z);

    K_values = NaN(n_k,n_z); 

    while max(max(abs(V_guess-V_opt))) > tolerance
    
        V_guess = V_opt;
    
        for i = 1:n_k
            for j = 1:n_z
                argument = z(j)*(k(i)^alpha) - (k - (1-delta)*k(i)) - (psi0/(2*k(i)).*(k - (1-delta)*k(i)).^2) ...
                    - psi1(index)*k(i).*(k~=(1-delta)*k(i))   + (V_guess*Pi(:,j))./(1+r);
                [V_opt(i,j), maxposition] = max(argument);
                K_values(i,j) = k(maxposition);
            end
        end
    
    end

    V(index) = mean(mean(V_opt));
    investment(index) = mean(mean((K_values - (1-delta).*k)));
end

figure(1)
plot(psi1, V, '-b', 'LineWidth', 2);hold on
grid on
xlabel('Fixed Investment Cost $\psi_1$','Interpreter','latex');
ylabel('Average Value Function $V(k,z)$','Interpreter','latex');

figure(2)
plot(psi1, investment, '--b', 'LineWidth', 2);hold on
grid on
xlabel('Fixed Investment Cost $\psi_1$','Interpreter','latex');
ylabel('Average Optimal Investment $I$','Interpreter','latex');

%% LOOP on: r
V = NaN(N,1);
investment = NaN(N,1);
for index = 1:N
    
    alpha = 0.35;
    delta = 0.15;
    rho_z = 0.7;
    sigma_e = 0.15;
    psi0 = 0.01;
    psi1 = 0.005;
    
    [ z, Pi ] = tauchen(mu_z, rho_z, sigma_e, n_z, m);
    z = exp(z)';
    k_max = (  (z(end)*alpha/(1+r(index))) / ( 1- (1-delta)/(1+r(index))  ) )^(1/(1-alpha));
    k_min = 0.1*(  (z(1)*alpha/(1+r(index))) / ( 1- (1-delta)/(1+r(index))  ) )^(1/(1-alpha));
    k = linspace(k_min, k_max, n_k)';

    V_guess = zeros(n_k,n_z);
    V_opt = ones(n_k,n_z);

    K_values = NaN(n_k,n_z); 

    while max(max(abs(V_guess-V_opt))) > tolerance
    
        V_guess = V_opt;
    
        for i = 1:n_k
            for j = 1:n_z
                argument = z(j)*(k(i)^alpha) - (k - (1-delta)*k(i)) - (psi0/(2*k(i)).*(k - (1-delta)*k(i)).^2) ...
                    - psi1*k(i).*(k~=(1-delta)*k(i))   + (V_guess*Pi(:,j))./(1+r(index));
                [V_opt(i,j), maxposition] = max(argument);
                K_values(i,j) = k(maxposition);
            end
        end
    
    end

    V(index) = mean(mean(V_opt));
    investment(index) = mean(mean((K_values - (1-delta).*k)));
end

figure(1)
plot(r, V, '-b', 'LineWidth', 2);hold on
grid on
xlabel('Discount Rate $r$','Interpreter','latex');
ylabel('Average Value Function $V(k,z)$','Interpreter','latex');

figure(2)
plot(r, investment, '--b', 'LineWidth', 2);hold on
grid on
xlabel('Discount Rate $r$','Interpreter','latex');
ylabel('Average Optimal Investment $I$','Interpreter','latex');
