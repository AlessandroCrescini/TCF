%% PROBLEM 2, question 6

rng(1);

alpha = 0.35;
r = 0.04;

N = 10000;
z = randn(N,1);

investment = ((1+r)./(exp(z).*alpha)).^(1/(alpha-1)) ;
av_investment = mean(investment);

figure(1)
histogram(investment,'Normalization','probability');
grid on
xlim([0,5]);
xlabel('$I$','Interpreter', 'latex'); 
ylabel('pdf');
hold on;
xline(av_investment, '-r', {'Mean'});


%% PROBLEM 2, question 7

rng(1);

n = 10;
alpha = linspace(0.15,0.55,n);
r = 0.04;

N = 10000;
z = randn(N,1);

av_investment = zeros(n,1);

for i=1:n
    investment = ((1+r)./(exp(z).*alpha(i))).^(1/(alpha(i)-1)) ;
    av_investment(i) = mean(investment);
    
end


figure(1)
plot(alpha, av_investment, '-b');
grid on
xlabel('$\alpha$','Interpreter', 'latex'); 
ylabel('Mean Investment');



%% PROBLEM 2, question 8

rng(1);

alpha = 0.35;
n = 10;
r = linspace(0,0.1,n);

N = 10000;
z = randn(N,1);

av_investment = zeros(n,1);

for i=1:n
    investment = ((1+r(i))./(exp(z).*alpha)).^(1/(alpha-1)) ;
    av_investment(i) = mean(investment);
    
end


figure(1)
plot(r, av_investment, '-b');
grid on
xlabel('$r$','Interpreter', 'latex'); 
ylabel('Mean Investment');






%% PROBLEM 2, question 9

rng(1);

n = 10;
alpha = 0.35;
r = 0.04;
mu = linspace(-1,1,n);

N = 10000;

av_investment = zeros(n,1);

for i=1:n
    z = randn(N,1) + mu(i);
    investment = ((1+r)./(exp(z).*alpha)).^(1/(alpha-1)) ;
    av_investment(i) = mean(investment);
    
end


figure(1)
plot(mu, av_investment, '-b');
grid on
xlabel('$\mu$','Interpreter', 'latex'); 
ylabel('Mean Investment');

%% PROBLEM 2, question 10

rng(1);

n = 10;
alpha = 0.35;
r = 0.04;
sigma = linspace(0.5,1.5,n);

N = 10000;

av_investment = zeros(n,1);

for i=1:n
    z = randn(N,1).*sigma(i);
    investment = ((1+r)./(exp(z).*alpha)).^(1/(alpha-1)) ;
    av_investment(i) = mean(investment);
    
end


figure(1)
plot(sigma, av_investment, '-b');
grid on
xlabel('$\sigma_z$','Interpreter', 'latex'); 
ylabel('Mean Investment');









