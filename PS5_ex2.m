%% Problem Set 5

%% SOLVING THE HENNESSY - WHITED (2005) MODEL 

clear; clc; close all;

% Initialization
alpha = 0.65;
delta = 0.15;
mu_z = 0;
rho_z = 0.7;
sigma_e = 0.15;
psi0 = 0.15;
psi1 = 0.005;
r = 0.02;
tau = 0.25;
chi = 0.5;
lambda0 = 0.1;
lambda1 = 0.15;

tolerance = 0.001;
n_z = 9;
m = 4;
n_k = 31;
n_d = 15;

% Discretization of the productivity shock z
[ z, Pi ] = tauchen(mu_z, rho_z, sigma_e, n_z, m);
z = exp(z)';

% Grid for capital k
k_max = (  (z(end)*alpha/(1+r)) / ( 1- (1-delta)/(1+r)  ) )^(1/(1-alpha));
k = [ NaN(1,n_k-1), k_max]';
for i = (n_k-1):-1:1
    k(i) = (1-delta)*k(i+1);
end

% Grid for debt d
d_max = 0.3*(1-chi)*k_max;
d_min = -0.1*(1-chi)*k_max;
d = linspace(d_min, d_max, n_d);


% The value function is represented by a  n_k x n_z x n_d  3D-matroid
% V(i,j,l) is the element V(k,z,d) with k(i), z(j) and d(l)

V_guess = zeros(n_k,n_z,n_d);
V_opt = ones(n_k,n_z,n_d);

K_values = NaN(n_k,n_z,n_d);  % To store the optimal future capital k'
D_values = NaN(n_k,n_z,n_d);  % To store the optimal future debt d'


% Iterate until the algorithm convergence limit goes below the threshold
% For each iteration step, V_opt(k,z,d) is computed element-by-element
% For the maximization grid search, we loop on a matrix having 
% k' on the rows, d' on the columns

while max(abs(V_guess-V_opt),[],'all') > tolerance
    
    V_guess = V_opt;
    
    for l = 1:n_d
        for i = 1:n_k
            for j = 1:n_z
                
                I = repmat(k,1,n_d) - (1-delta)*k(i);
                
                e = (1-tau)*z(j)*(k(i)^alpha) - I - (psi0/(2*k(i))).*(I.^2) ...
                    - psi1*k(i).*(repmat(k,1,n_d)~=(1-delta)*k(i)) -  d(l) ...
                    + repmat(d,n_k,1)./(1+r*(1-tau));
                
                Lambda = (lambda0 + lambda1.*abs(e)).*(e<0);
                
                J = NaN(n_k,n_d);
                for p = 1:n_d
                    J(:,p) = V_guess(:,:,p)*Pi(:,j);
                end
                
                argument = e - Lambda + J./(1+r);
                
                % Verify the borrowing constraint d' <= (1-chi)k'
                max_value = max(argument, [], 'all');
                [row, column] = find(argument==max_value);
                while d(column) > (1-chi)*k(row)
                    argument(row,column) = NaN;
                    max_value = max(argument, [], 'all');
                    [row, column] = find(argument==max_value);
                end
                
                V_opt(i,j,l) = max_value;
                K_values(i,j,l) = k(row);
                D_values(i,j,l) = d(column);
                
            end
        end
    end
    
end


%% SIMULATION

N = 10000;              % number of simulations
T = 120;                % simulation length

rng(1);
u = rand(N,T);          % random to simulate shocks
F = cumsum(Pi,1);       % Empirical cdf of transition probabilities

z_start = z(fix(n_z/2)); % Middle point of the grid


% z_sim contains the starting point for z and the following (N,T) simulated
% values. d_sim and k_sim contains their starting points and the optimal
% policies for all the (N,T) scenarios. V_sim collects the firm value V in each
% scenario

z_sim = [ z_start.*ones(N,1), NaN(N,T) ];

% Simulation of the shocks
for i = 2:(T+1)
    for j = 1:N
        initial_state = find(z==z_sim(j,i-1));
        aux = F(:,initial_state).*(F(:,initial_state) > u(j,i-1));        % Boolean array
        
        aux(aux==0) = NaN;
        [temp, final_state] = min(aux);
        z_sim(j,i) = z(final_state);
    end
end


% Initialization of k_sim, d_sim, V_sim

d_start = (d(1)+d(end))/2;
k_start = (k(1)+k(end))/2;

d_sim = [ d_start.*ones(N,1), NaN(N,T) ];
k_sim = [ k_start.*ones(N,1), NaN(N,T) ];

% To semplify the code, I wrote a dedicate function for the linear
% interpolation in another M-file

V_start = interpolate2d(V_opt,k,d,z,k_start,d_start,z_start);
V_sim = [ V_start.*ones(N,1), NaN(N,T) ];


% Computing k_sim, d_sim and V_sim

for t = 2:(T+1)
    for n = 1:N
    
        k_sim(n,t) = interpolate2d(K_values,k,d,z, ...
            k_sim(n,t-1), d_sim(n,t-1), z_sim(n,t-1));
     
        d_sim(n,t) = interpolate2d(D_values,k,d,z, ...
            k_sim(n,t-1), d_sim(n,t-1), z_sim(n,t-1));
     
        V_sim(n,t) = interpolate2d(V_opt,k,d,z, ...
            k_sim(n,t), d_sim(n,t), z_sim(n,t));
        
    end
    
end


%% COMPUTING MOMENTS

% Definition of the interesting variables

net_profits = (1-tau).*z_sim.*(k_sim.^alpha);
net_profitability = net_profits./k_sim;
investment = k_sim(:,2:end) - (1-delta).*k_sim(:,1:end-1);
investment_rate = investment./k_sim(:,1:end-1);
debt = d_sim;
book_leverage = d_sim./k_sim;
tobin = (V_sim + d_sim)./k_sim;

% Consider simulated points from the last 20 years only

net_profits = net_profits(:,end-19:end);
net_profitability = net_profitability(:,end-19:end);
investment = investment(:,end-19:end);
investment_rate = investment_rate(:,end-19:end);
debt = debt(:,end-19:end);
book_leverage = book_leverage(:,end-19:end);
tobin = tobin(:,end-19:end);

% Computing the moments
Mean = [mean(net_profits,'all'); mean(net_profitability,'all'); mean(investment,'all'); ...
    mean(investment_rate,'all'); mean(debt,'all'); mean(book_leverage,'all'); mean(tobin,'all') ];
 
Variance = [var(net_profits,0,'all'); var(net_profitability,0,'all'); var(investment,0,'all'); ...
    var(investment_rate,0,'all'); var(debt,0,'all'); var(book_leverage,0,'all'); var(tobin,0,'all') ];

Minimum = [min(net_profits,[],'all'); min(net_profitability,[],'all'); min(investment,[],'all'); ...
    min(investment_rate,[],'all'); min(debt,[],'all'); min(book_leverage,[],'all'); min(tobin,[],'all') ];

Quantile_25 = [quantile(net_profits,0.25,'all'); quantile(net_profitability,0.25,'all'); quantile(investment,0.25,'all'); ...
    quantile(investment_rate,0.25,'all'); quantile(debt,0.25,'all'); quantile(book_leverage,0.25,'all'); quantile(tobin,0.25,'all') ];

Median = [median(net_profits,'all'); median(net_profitability,'all'); median(investment,'all'); ...
    median(investment_rate,'all'); median(debt,'all'); median(book_leverage,'all'); median(tobin,'all') ];

Quantile_75 = [quantile(net_profits,0.75,'all'); quantile(net_profitability,0.75,'all'); quantile(investment,0.75,'all'); ...
    quantile(investment_rate,0.75,'all'); quantile(debt,0.75,'all'); quantile(book_leverage,0.75,'all'); quantile(tobin,0.75,'all') ];

Maximum = [max(net_profits,[],'all'); max(net_profitability,[],'all'); max(investment,[],'all'); ...
    max(investment_rate,[],'all'); max(debt,[],'all'); max(book_leverage,[],'all'); max(tobin,[],'all') ];


% Create a table, then format it for TeX environment
% The command 'table2latex' comes from official 2022 release by V.M.Cagigal

moments = table(Mean, Variance, Minimum, Quantile_25, Median, Quantile_75, Maximum);

table2latex(moments);
