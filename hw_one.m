% ADSP - HW1
% Ugur Ozcan

clear all
clc

%% 1.1 Noise Generation

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

figure; hold on
plot(L, mean(mse(1,:,:),3).', 'LineWidth', 2, "DisplayName",strcat("\rho=",string(rho(1))))
plot(L, mean(mse(5,:,:),3).', 'LineWidth', 2, "DisplayName",strcat("\rho=",string(rho(5))))
plot(L, mean(mse(10,:,:),3).', 'LineWidth', 2, "DisplayName",strcat("\rho=",string(rho(10))))
xlabel('Number of samples [L]');
ylabel('MSE');
grid
legend

figure; hold on
plot(rho, mean(mse(:,1,:),3), 'LineWidth', 2, "DisplayName",strcat("L=",string(L(1))))
plot(rho, mean(mse(:,5,:),3), 'LineWidth', 2, "DisplayName",strcat("L=",string(L(5))))
plot(rho, mean(mse(:,10,:),3), 'LineWidth', 2, "DisplayName",strcat("L=",string(L(10))))
xlabel('Correlation coefficient [\rho]');
ylabel('MSE');
grid
legend

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


%% deneme3

N = size(X',1);
M = size(X',2);
k=3

X_conv = zeros(N+k, k*M );
for i = 1 : k
    X_conv( i : i+N-1, (i-1)*M+1 : i*M  ) = X';
end


%%
clear
clc



N = 5;
SNRdb = -10:2:30;
alpha = linspace(0, .99, N);
Q = linspace(10, 50, 5);
beta = 0.1;
M = 4;
rho = 0.15;
K = 2;

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


% for r = 1:length(rho)
%     for q = 1:length(Q)
%         CW{r,q} = toeplitz([1.*ones(1,K), zeros(1,Q(q)-1), rho(r).*ones(1,K), zeros(1,Q(q)-1), rho(r)^2.*ones(1,K), zeros(1,Q(q)-1), rho(r)^3.*ones(1,K), zeros(1,Q(q)-1)]);
%     end
% end
for r = 1:length(rho)
    for q = 1:length(Q)
        CW{r,q} = toeplitz([1, zeros(1,Q(q)-2+K), rho(r), zeros(1,Q(q)-2+K), rho(r)^2, zeros(1,Q(q)-2+K), rho(r)^3, zeros(1,Q(q)-2+K)]);
    end
end


IT = 20
status = 0
for s = 1:length(SNR)
    for a = 1:length(alpha)
        for q = 1:length(Q)
            for it = 1:IT
                X = sqrt(SNR(s)).*randn(M,Q(q));
                X_conv = zeros(Q(q)+K-1, K*M );
                for i = 1 : K
                    X_conv( i : i+Q(q)-1, (i-1)*M+1 : i*M  ) = X';
                    h1(i,:) = reshape(H{a}(:,:,i),[1,numel(H{a}(:,:,i))]);
                end

                h = h1(:);
                W = sqrtm(C) * randn(M, Q(q)+K-1);
                w = reshape(W',[numel(W) 1]);

                X_bp = kron(eye(4),X_conv);
                y = X_bp*h + w;

                h_hat_mle = (X_bp'*(CW{1,q}\X_bp))\(X_bp'*(CW{1,q}\y));
                err_mle(s,a,q,it) = mean((h_hat_mle - h).^2);

            end
        end
        if it*s*a*q/(length(Q)*length(SNR)*length(alpha)*IT)*100>status
            status = it*s*a*q/(length(Q)*length(SNR)*length(alpha)*IT)*100;
        end
        disp(status)
    end
end

err = mean(err_mle,4);


%loglog(SNR,err(:,1,:))
%hold on
loglog(SNR,reshape(err(:,4,:),21,5))
hold on
% loglog(SNR,CRB2,"DisplayName","CRB")


%% MIMO DEC.
clear
clc

P = 201;
Q = linspace(10, P-1, 20);
rho=.1;
alpha=.5;



SNRdb = -10:2:30;
beta = 0.99;
M = 4;
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


CW2 = toeplitz([1, zeros(1,P-1), rho, zeros(1,P-1), rho^2, zeros(1,P-1), rho^3, zeros(1,P-1)]);


for s = 1:length(SNR)
    for q = 20 %1:length(Q)
        CW = toeplitz([1, zeros(1,Q(q)-1), rho, zeros(1,Q(q)-1), rho^2, zeros(1,Q(q)-1), rho^3, zeros(1,Q(q)-1)]);
        for it = 1:1
            Xall = sqrt(SNR(s)).*randn(M,P);
            X = Xall(:,1:Q(q));
            x = X(:);
            xall = Xall(:);
            h = H{a}(:);
            Wall = sqrtm(C) * randn(M, P);
            W = Wall(:,1:Q(q));
            w = reshape(W',[numel(W) 1]);
            wall = reshape(Wall',[numel(Wall) 1]);
            X_bp = kron(eye(4),X');
            X_bpall = kron(eye(M),Xall');

            y = X_bp*h + w;

            h_mle = (X_bp'*(CW\X_bp))\(X_bp'*(CW\y));
            H_est = reshape(h_mle,[4 4]);
            %                 x_zf = y_all*pinv(h);
            %                 err_zf(q,s,a,it) = mean(mean((x_zf - X_bpall).^2));


            %HH = kron(eye(Q(q)),H_est');
            HHall = kron(eye(P),H_est');
            HHdn = kron(eye(P),H{a}');
            %y2 = HH*x+w;
            y_all2 = HHall*xall+wall;
            x_mle = (HHall'*(CW2\HHall))\(HHall'*(CW2\y_all2));
            err_mle(q,s,it) = mean(mean((x_mle - xall).^2));
            CC=sqrt(SNR(s))*eye(size(HHall));
            %             x_mmse = CC*HHall'*inv(HHall*CC*HHall'+CW2)*y_all2;
            x_mmse = inv(HHall'*inv(CW2)*HHall+inv(CC))*HHall'*(CW2\y_all2);
            err_mmse(q,s,it) = mean(mean((x_mmse - xall).^2));
            x_dn = inv(HHdn'*inv(CW2)*HHdn)*(HHdn'*(CW2\y_all2));
            err_dn(q,s,it) = mean(mean((x_dn - xall).^2));
        end
    end

end
%%
loglog(SNR,mean(err_dn,3))
hold on
loglog(SNR,mean(err_mle,3))