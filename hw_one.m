% ADSP - HW1
% Ugur Ozcan

clear all
clc

%% 1.1 Noise Generation

N = 4;
L = round(linspace(10,200,10));
rho = linspace(0,0.99,10);

% rhoyu düzelt
for i = 1 : length(rho)
    for j = 1 : length(L)
        for mc = 1 : 500
            Cw = toeplitz([1, rho(i)*ones(1,N-1)]);
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

Cw = toeplitz([1, 0*ones(1,M-1)]);
W= sqrtm(Cw) * randn(M, Q);

for i=1:M
    for j=1:M
        for n=1:N
            H(i,j) = alpha(n)^(abs(i-j))*beta;
        end
    end
end
%%

Q = linspace(1, 50, N);

for q = Q
    for it = 100
        X = randi(2,M,Q)-1;
        h = H(:);
        w = reshape(W',[numel(W) 1]);
        X_bp = kron(eye(4),X');
        
        y = X_bp*h + w;
        
        h_hat_mle = (X_bp'*(Cw\X_bp))\(X_bp'*(Cw\y));
        err_mle(Q,it) = mean((h_hat_mle - h).^2);
    end
end

err = mean(err_mle);
