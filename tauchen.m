
function [s, Pi] = tauchen(mu,rho,sigma,N,m)

% This function discretizes a continuous AR(1) process with the Tauchen
% (1986) algorithm. This code is rearranged from J.H.Lang (2010)

% INPUT
% The AR(1) process is: 
% x(t) = mu + rho*x(t-1) + epsilon, with epsilon ~ N(0,sigma^2)
%
% N  : number of points of the grid
% m  : grid width with respect of mu, in unit of std.dev.

% OUTPUT
% s  : column vector of the discretized x, in ascending order
% Pi : NxN matrix of transition probabilities.
%      Final state in rows, initial state in columns


s       = NaN(N,1);
Pi      = NaN(N,N);

mu_x = mu/(1-rho);
sigma_x = sigma/sqrt(1-rho^2);

s(1)    = mu_x - m*sigma_x;
s(N)    = mu_x + m*sigma_x;
step    = (s(N)-s(1))/(N-1);


% Compute the intermediate values of the grid
for i = 2:(N-1)
   s(i) = s(i-1) + step; 
end


% Compute the transition probabilities matrix
for i = 1:N
    for j = 1:N
        Pi(j,i) = normcdf(s(j)+0.5*step, rho*s(i), sigma_x) - ...   
            normcdf(s(j)-0.5*step, rho*s(i), sigma_x);
    end
end



