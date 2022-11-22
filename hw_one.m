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

% %% 
% clear
% 
% 
% N = 5;
% SNRdb = -10:2:30;
% %Q = linspace(1, 50, N);
% Q = 10;
% alpha = linspace(0, .99, N);
% beta = 0.9;
% M = 4;
% 
% bits = randi(2,1,1e3)-1;
% training_bits = bits(1:Q);
% 
% SNR = db2pow(SNRdb);
% 
% Cw = toeplitz([1, 0*ones(1,M-1)]);
% W= sqrtm(Cw) * randn(M, Q);
% 
% for i=1:M
%     for j=1:M
%         for n=1:N
%             H(i,j) = alpha(n)^(abs(i-j))*beta;
%         end
%     end
% end
%%
% 
% Q = linspace(10, 50, N);
% 
% for q = Q
%     for it = 100
%         X = randi(2,M,q)-1;
%         h = H(:);
%         w = reshape(W',[numel(W) 1]);
%         X_bp = kron(eye(4),X');
%         
%         y = X_bp*h + w;
%         
%         h_hat_mle = (X_bp'*(Cw\X_bp))\(X_bp'*(Cw\y));
%         err_mle(Q,it) = mean((h_hat_mle - h).^2);
%     end
% end
% 
% err = mean(err_mle);

%% 
clear
clc


N = 5;
SNRdb = -10:2:30;
alpha = linspace(0, .99, N);
Q = linspace(10, 50, 5);
beta = 0.9;
M = 4;
rho = 0.15;
K = 1;

SNR = 10.^(SNRdb/10);


for i=1:M
    for j=1:M
        for a=1:length(alpha)
            for k=1:K
                H{a}(i,j,k) = alpha(a)^(abs(i-j))*beta^(k);
            end
        end
    end
end



for input=1:M       
    for output=1:M  
        C(output,input)=rho^(abs(output-input));
    end
end


for r = 1:length(rho)
    for q = 1:length(Q)
        CW{r,q} = toeplitz([1, zeros(1,Q(q)-1), rho(r), zeros(1,Q(q)-1), rho(r)^2, zeros(1,Q(q)-1), rho(r)^3, zeros(1,Q(q)-1)]);
    end
end

for s = 1:length(SNR)
    for a = 1:length(alpha)
        for q = 1:length(Q)
            for it = 1:10
    %           X = 2*randi([0,1],M,Q(q))-1;
                X = sqrt(SNR(s)).*randn(M,Q(q));
    %           X = sqrt(SNR(s)).*X;
                h = H{a}(:);
                W = sqrtm(C) * randn(M, Q(q));
                w = reshape(W',[numel(W) 1]);
                X_bp = kron(eye(4),X');
                
                y = X_bp*h + w;
                
                h_hat_mle = (X_bp'*(CW{1,q}\X_bp))\(X_bp'*(CW{1,q}\y));
                err_mle(q,s,a,it) = nanmean((h_hat_mle - h).^2);
    
    %           crb_mat{s,q} =  inv(X_bp'*inv(CW{1,q})*X_bp);
                crb_mat{q,s,a} =  inv(X_bp'*X_bp);
    
                CRB(q,s,a,it) = trace(crb_mat{q,s,a})/length(crb_mat{q,s,a});
            end
        end
    end
end
% err = sum(err_mle,3)/numel(err_mle(1,1,:));
err = nanmean(err_mle,4);
CRB2 = nanmean(CRB,4);

for q = 1:length(Q)
    CRBa(q,:)=10.^(-SNRdb/10)/Q(q);
end

%%

loglog(SNR,err(:,:,1))
hold on
%loglog(SNR,CRB2,"DisplayName","CRB")
loglog(SNR,CRBa,"DisplayName","CRBa")

%%
loglog(SNR,reshape(err(5,:,:),21,5))
hold on
%loglog(SNR,CRB2,"DisplayName","CRB")
% loglog(SNR,CRBa,"DisplayName","CRBa")


%%

N = size(X',1);
M = size(X',2);
k=3

X_conv = zeros(N+k, k*M );
for i = 1 : k
    X_conv( i : i+N-1, (i-1)*M+1 : i*M  ) = X';
end