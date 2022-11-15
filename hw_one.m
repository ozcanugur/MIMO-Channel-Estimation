% ADSP - HW1
% Ugur Ozcan

clear all
clc

%% 1.1 Noise Generation

N = 2;
L = round(linspace(10,200,10));
rho = linspace(0,0.99,10);

for i = 1 : length(rho)
    for j = 1 : length(L)
        for mc = 1 : 500
            Cw = [1, rho(i); rho(i) 1];
            w = sqrtm(Cw) * randn(N, L(j));
            Cw_est = (w * w') / L(j);
            errC = (Cw - Cw_est).^2;
            mse(i,j,mc) = mean(errC(:));
        end
    end
end

figure; hold on
plot(L, mean(mse,3).', 'LineWidth', 2)
xlabel('Number of samples [L]');
ylabel('MSE');


figure; hold on
plot(rho, mean(mse(:,[1, 10],:),3), 'LineWidth', 2)
xlabel('Correlation coefficient [\rho]');
ylabel('MSE');

%% 
clear


N = 5;
SNRdb = -10:2:30;
%Q = linspace(1, 50, N);
Q = 10;
alpha = linspace(0, .99, N);
beta = 0.9;
M = 4;

bits = randi(2,1,1e3)-1;
training_bits = bits(1:Q);

SNR = db2pow(SNRdb);


for i=1:M
    for j=1:M
        for n=1:N
            H{i,j,n} = alpha(n)^(abs(i-j))*beta;
        end
    end
end


H_est = ones(M,M)

