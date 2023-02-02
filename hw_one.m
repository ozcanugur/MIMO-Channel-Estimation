% ADSP - HW1
% Ugur Ozcan

%% 1.1 Noise Generation

clear all
clc

N = 4;
M = 4;
L = round(linspace(10,200,10));
rho = linspace(0,0.99,10);

for i = 1 : length(rho)
    for j = 1 : length(L)
        for input=1:M
            for output=1:N
                Cw(output,input)=rho(i)^(abs(output-input));
            end
        end
        for mc = 1 : 500
            w = sqrtm(Cw) * randn(N, L(j));
            Cw_est = (w * w') / L(j);
            errC = (Cw - Cw_est).^2;
            mse(i,j,mc) = mean(errC(:));
        end
    end
end

figure; 
loglog(L, mean(mse(1,:,:),3).', 'LineWidth', 2, "DisplayName",strcat("\rho=",string(rho(1))))
hold on
loglog(L, mean(mse(5,:,:),3).', 'LineWidth', 2, "DisplayName",strcat("\rho=",string(rho(5))))
loglog(L, mean(mse(10,:,:),3).', 'LineWidth', 2, "DisplayName",strcat("\rho=",string(rho(10))))
xlabel('Number of samples [L]');
ylabel('MSE');
grid
legend

figure;
plot(rho, mean(mse(:,1,:),3), 'LineWidth', 2, "DisplayName",strcat("L=",string(L(1))))
hold on
plot(rho, mean(mse(:,5,:),3), 'LineWidth', 2, "DisplayName",strcat("L=",string(L(5))))
plot(rho, mean(mse(:,10,:),3), 'LineWidth', 2, "DisplayName",strcat("L=",string(L(10))))
xlabel('Correlation coefficient [\rho]');
ylabel('MSE');
grid
legend


%% 1.2 MIMO estimation 1

clear
clc

SNRdb = -10:2:30;
alpha = linspace(0, .99, 5);
Q = linspace(10, 50, 5);
beta = 0.9;
M = 4;
rho = 0.9;
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
        CW{r,q} = toeplitz([1, zeros(1,Q(q)-1),...
            rho(r), zeros(1,Q(q)-1), rho(r)^2,...
            zeros(1,Q(q)-1), rho(r)^3, zeros(1,Q(q)-1)]);
    end
end

for s = 1:length(SNR)
    for a = 1:length(alpha)
        for q = 1:length(Q)
            for it = 1:200
                X = sqrt(SNR(s)).*randn(M,Q(q));
                h = H{a}(:);
                W = sqrtm(C) * randn(M, Q(q));
                w = reshape(W',[numel(W) 1]);
                X_bp = kron(eye(4),X');
                y = X_bp*h + w;
                h_hat_mle = (X_bp'*(CW{1,q}\X_bp))\(X_bp'*(CW{1,q}\y));
                err_mle(q,s,a,it) = nanmean((h_hat_mle - h).^2);
                crb_mat{q,s,a} =  inv(X_bp'*(CW{1,q}\X_bp));
                CRB(q,s,a,it) = trace(crb_mat{q,s,a})/length(crb_mat{q,s,a});
            end
        end
    end
    s
end

err = nanmean(err_mle,4);
CRB2 = nanmean(CRB,4);

%% Figure for fixed Alpha=0

semilogy(SNRdb,err(:,:,1),'LineWidth', 2)
hold on
semilogy(SNRdb,CRB2(:,:,1),'--','LineWidth', 1)
legend([strcat("Q= ",string(Q)) strcat("CRB, Q=",string(Q))])
%loglog(SNR,CRBa,"DisplayName","CRBa")
title("SNR vs Error values for Changing Pilot Numbers when \alpha=0")
xlabel('SNR [dB]');
ylabel('MSE');
grid

%% Figure for fixed Q=10

semilogy(SNRdb,reshape(err(1,:,:),21,5),'LineWidth', 2)
hold on
semilogy(SNRdb,reshape(CRB2(1,:,:),21,5),'--','LineWidth', 1)
legend([strcat("\alpha= ",string(alpha)) strcat("CRB, \alpha=",string(alpha))])
%loglog(SNR,CRBa,"DisplayName","CRBa")
title("SNR vs Error values for Changing Pilot Numbers when Q=10")
xlabel('SNR [dB]');
ylabel('MSE');
grid

%% 1.2 MIMO estimation 2

clear
clc

SNRdb = -10:2:30;
alpha = linspace(0, .99, 5);
Q = linspace(10, 50, 5);
beta = 0.1;
M = 4;
rho = 0.1;
K = 4;

SNR = 10.^(SNRdb/10);

for input=1:M
    for output=1:M
        C(output,input)=rho^(abs(output-input));
    end
end

for r = 1:length(rho)
    for q = 1:length(Q)
        CW{r,q} = toeplitz([1, zeros(1,Q(q)-2+K),...
            rho(r), zeros(1,Q(q)-2+K), rho(r)^2,...
            zeros(1,Q(q)-2+K), rho(r)^3, zeros(1,Q(q)-2+K)]);
    end
end

CRB = zeros(numel(SNR),numel(alpha),numel(Q));
err = zeros(numel(SNR),numel(alpha),numel(Q));

IT = 100;
status = 0;
for s = 1:length(SNR)
    for a = 1:length(alpha)
        for q = 1:length(Q)
            for it = 1:IT
                X = sqrt(SNR(s)).*randn(M,Q(q));
                X_conv = zeros(Q(q)+K-1, K*M );
                for i = 1 : K
                    X_conv( i : i+Q(q)-1, (i-1)*M+1 : i*M  ) = X';
                end
                H = toeplitz(alpha(a).^(0:M-1));
                h2 = H(:);
                env = beta.^(0:K-1);
                H_mem = h2*env;
                h_mem = zeros(4*M^2,1);
                for i = 0:M-1
                    H_current = H_mem( (M*i+1 : M*(i+1)) , : );
                    h_mem( M*4*i+1 : M*4*(i+1)  ) = H_current(:);
                end
                W = sqrtm(C) * randn(M, Q(q)+K-1);
                w = reshape(W',[numel(W) 1]);
                h = h_mem;
                X_bp = kron(eye(4),X_conv);
                y = X_bp*h + w;
                h_hat_ls = pinv(X_bp)*y;
                err_ls(s,a,q,it) = mean((h_hat_ls - h).^2);
                crb_mat{s,a,q} =  inv(X_bp'*X_bp);
                mean_crb = trace(crb_mat{s,a,q})/size(crb_mat{s,a,q},1);  
                CRB(s,a,q) = CRB(s,a,q) + (mean_crb / IT);
                err(s,a,q) = err(s,a,q) + (err_ls(s,a,q,it)/IT);
            end
        end
        if it*s*a*q/(length(Q)*length(SNR)*length(alpha)*IT)*100>status
            status = it*s*a*q/(length(Q)*length(SNR)*length(alpha)*IT)*100;
        end
        disp(status)
    end
end

loglog(SNR,reshape(err(:,3,:),21,length(Q)))

%% Figure for fixed Alpha=0
semilogy(SNRdb,reshape(err(:,1,:),21,length(Q)),'LineWidth', 2)
hold on
semilogy(SNRdb,reshape(CRB(:,1,:),21,length(Q)),'--','LineWidth', 1)
legend([strcat("Q= ",string(Q)) strcat("CRB, Q=",string(Q))])
%loglog(SNR,CRBa,"DisplayName","CRBa")
title("SNR vs Error values for Changing Pilot Numbers for \beta=0.1 \alpha=0")
xlabel('SNR [dB]');
ylabel('MSE');
grid

%% Figure for fixed Q=50
semilogy(SNRdb,reshape(err(:,:,5),21,length(alpha)),'LineWidth', 2)
hold on
semilogy(SNRdb,reshape(CRB(:,:,5),21,length(alpha)),'--','LineWidth', 1)
legend([strcat("\alpha= ",string(alpha)) strcat("CRB, \alpha=",string(alpha))])
%loglog(SNR,CRBa,"DisplayName","CRBa")
title("SNR vs Error values for Changing Pilot Numbers for \beta=0.1 Q=50")
xlabel('SNR [dB]');
ylabel('MSE');
grid

%% 1.3 MIMO Deconvolution (Ideal Lower Bound)

clear
clc

P = 201;
Q = linspace(10, P-1, 20);
rho=.1;
alpha=.5;
beta = 1;

SNRdb = -10:2:30;
M = 4;
K = 1;

SNR = 10.^(SNRdb/10);

for i=1:M
    for j=1:M
        for a=1:length(alpha)
            for k=1:K
                H1{a}(i,j,k) = alpha(a)^(abs(i-j))*beta^(k);
            end
        end
    end
end

for input=1:M
    for output=1:M
        C(output,input)=rho^(abs(output-input));
    end
end

for s = 1:length(SNR)
    for a = 1:length(alpha)
        for q = 1:length(Q)
            for it = 1:100
                Xall = randn(M,P);
                X = Xall(:,1:Q(q));
                Cw = C./SNR(s);
                W = sqrtm(Cw) * randn(M, Q(q));
                H = H1{a};
                Y = H*X + W;
                % MLE Estimator of X
                X_ml = (H'*(Cw\H))\(H'*(Cw\Y));
                % Winer Filter (LMMSE)
                X_mmse = (H'*(Cw\H)+eye(M))\(H'*(Cw\Y));
                X_ls = H\Y;
                err_ls(s,q,it) = mean((X_ls(:) - X(:)).^2);
                err_mle(s,q,it) = mean((X_ml(:) - X(:)).^2);
                err_mmse(s,q,it) = mean((X_mmse(:) - X(:)).^2);
            end
        end
    end
    s
end

err = mean(err_mle,3);

%% Figure for fixed Alpha=0
semilogy(SNRdb,err(:,:),'LineWidth', 2)
hold on
semilogy(SNRdb,CRB2(:,:),'--','LineWidth', 1)
legend([strcat("Q= ",string(Q)) strcat("CRB, Q=",string(Q))])
%loglog(SNR,CRBa,"DisplayName","CRBa")
title("SNR vs Error values for Changing Pilot Numbers for \beta=0.1 \alpha=0")
xlabel('SNR [dB]');
ylabel('MSE');
grid

%% 1.3 MIMO Deconvolution

clear
clc

P = 201;
Q = linspace(10, P-1, 20);
rho=.1;
alpha=.5;
beta = 1;

SNRdb = -10:2:30;
M = 4;
K = 1;

SNR = 10.^(SNRdb/10);

for i=1:M
    for j=1:M
        for a=1:length(alpha)
            for k=1:K
                H1{a}(i,j,k) = alpha(a)^(abs(i-j))*beta^(k);
            end
        end
    end
end

for input=1:M
    for output=1:M
        C(output,input)=rho^(abs(output-input));
    end
end

for s = 1:length(SNR)
    
    for a = 1:length(alpha)
        for q = 1:length(Q)
            %Xall = sqrt(SNR(s)).*randn(M,P);
            
            for it = 1:200
                Xall = randn(M,P);
                X = Xall(:,1:Q(q));
                Cw = C./SNR(s);
                W = sqrtm(Cw) * randn(M, P);
                
          
                H = H1{a};
                Y = H*Xall + W;
                
                [H_est,filter_err(s,q,it)] = MemorylessFilterEstimator(Q(q),SNR(s),X);

                % MLE Estimator of X
                X_ml = (H_est'*(Cw\H_est))\(H_est'*(Cw\Y));
                
                % Winer Filter (LMMSE)
                X_mmse = (H_est'*(Cw\H_est)+eye(M))\(H_est'*(Cw\Y));
                
                % Least Squares Estimator (LS)
                X_ls = H\Y;

                
                err_ls(s,q,it) = mean((X_ls(:) - Xall(:)).^2);
                err_mle(s,q,it) = mean((X_ml(:) - Xall(:)).^2);
                err_mmse(s,q,it) = mean((X_mmse(:) - Xall(:)).^2);
             
            end
        end
    end
    s
end

MSE_h = mean(filter_err,3);
MSE_ls = mean(err_ls,3);
MSE_mle = mean(err_mle,3);
MSE_mmse = mean(err_mmse,3);

decision_metric = MSE_h+MSE_mle;


%% Figure MLE

semilogy(SNRdb,MSE_mle,'LineWidth', 2)
hold on
legend(strcat("Q= ",string(Q)))
%loglog(SNR,CRBa,"DisplayName","CRBa")
title("SNR vs Deconvolution MSE for Changing Pilot Numbers with MLE")
xlabel('SNR [dB]');
ylabel('MSE');
grid

%% Figure MMSE

semilogy(SNRdb,MSE_mmse,'LineWidth', 2)
hold on
legend(strcat("Q= ",string(Q)))
title("SNR vs Deconvolution MSE for Changing Pilot Numbers with MMSE")
xlabel('SNR [dB]');
ylabel('MSE');
grid

%% Figure Decision Metric

semilogy(SNRdb,decision_metric,'LineWidth', 2)
hold on
legend(strcat("Q= ",string(Q)))
%loglog(SNR,CRBa,"DisplayName","CRBa")
title("SNR vs MSE_H+MSE_X for Changing Pilot Numbers")
xlabel('SNR [dB]');
ylabel('Decision Metric');
grid
%%

function [H_est,err_mle] = MemorylessFilterEstimator(Q,SNR,X)

    P = 201;
    rho=.1;
    alpha=.5;
    M = 4;
    K = 1;
     
    for i=1:M
        for j=1:M
            H(i,j) = alpha^(abs(i-j));
        end
    end
    
    for input=1:M
        for output=1:M
            C(output,input)=rho^(abs(output-input));
        end
    end
    
    C = C./SNR;

    for r = 1:length(rho)
        for q = 1:length(Q)
            CW{r,q} = toeplitz([1, zeros(1,Q(q)-1), rho(r),...
                zeros(1,Q(q)-1), rho(r)^2, zeros(1,Q(q)-1),...
                rho(r)^3, zeros(1,Q(q)-1)]);
        end
    end

    h = H(:);
    W = sqrtm(C) * randn(M, Q);
    w = reshape(W',[numel(W) 1]);
    X_bp = kron(eye(4),X');
    y = X_bp*h + w;
    h_hat_mle = (X_bp'*(CW{1,q}\X_bp))\(X_bp'*(CW{1,q}\y));
    H_est = reshape(h_hat_mle,4,4);
    err_mle = nanmean((h_hat_mle - h).^2);

end