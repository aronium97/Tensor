clc
close all
clear all

L = 54;
F = 210;
T = 100;

nodes = [[0,0.3];[0.18,0.28];[0.19, 0.05];[0.41,0.39];[0.6,0.21];[1,0.21];[0.91,0.4];[0.8,0.4];[0.39,0.51];[0.59,0.59];[0.6,0.6];[0.42,0.79];[0.76, 0.86];[0.92,0.98];[0.95,0.85]];

rho = 5;

% make random variables
sigma = 10^-2;
r = 5;
p = 0.05;
pi_missing_data = 0.8; % missing data pi, 1= no missing data

% make v_l_t
V = randn(L,T).*sigma^2;
w = randn(r,T);
U = randn(F,r)*1/F;
Z = U*w;
helpDistri = rand(F,T);
A = (helpDistri<(p/2))*(-1) + (helpDistri>=(p/2) & helpDistri<p)*(1) + (helpDistri>=p)*0;

for t = 1:T
    yAvailability = rand(L,1) < pi_missing_data;
    omega_t(:,:,t) = eye(L).*yAvailability;

    for l = 1:L
        omega_l(t,t,l) = yAvailability(l);
    end
end


for t = 1:T
    yAvailability = ones(L,1);
    omega_t1(:,:,t) = eye(L).*yAvailability;

    for l = 1:L
        omega_l1(t,t,l) = yAvailability(l);
    end
end



R = getR(L,F,nodes);

Y = (R*Z + R*A + V);


K_bsca = 100; % num iterations bcsa algorithm
K_bcd = 4; % num iterations batch bcd algorithm

lambda1 = 1 % am ehesten mit mu_soft connected % batch bcd: max(max(abs(R'*V)));
lambdastar = 1 % batch bcd: (sqrt(T) + sqrt(F)*sqrt(pi))*sigma%norm(V,1);

mu_soft_bsca = 1;  
mu_soft_bcd = 10;

% init P and Q at random
% X = LxT = PQ'
Q = randn(T,rho); 
P = randn(L,rho);%5*R*A*Q;%randn(L,rho);



[obj_value_bsca, time_value_bsca] = bsca_missing_data(P, Q, A, K_bsca, R, Y, omega_t, omega_l, lambdastar, lambda1, mu_soft_bsca);
plot(time_value_bsca,obj_value_bsca, "--")
set(gca, 'YScale', 'log')

hold on

[obj_value_bcd, time_value_bsd] = batch_bcd(P, Q, A, K_bcd, R, Y, omega_t, omega_l, lambdastar, lambda1, mu_soft_bcd);
plot(time_value_bsd, obj_value_bcd, "--")
set(gca, 'YScale', 'log')

hold on

[obj_value_bsca, time_value_bsca] = bsca_missing_data(P, Q, A, K_bsca, R, Y, omega_t1, omega_l1, lambdastar, lambda1, mu_soft_bsca);
plot(time_value_bsca, obj_value_bsca, "--")
set(gca, 'YScale', 'log')

hold on

[obj_value_bcd, time_value_bsd] = batch_bcd(P, Q, A, K_bcd, R, Y, omega_t1, omega_l1, lambdastar, lambda1, mu_soft_bcd);
plot(time_value_bsd, obj_value_bcd, "--")
set(gca, 'YScale', 'log')


ylabel("obj value")
xlabel("seconds")
legend("BSCA (p=" + string(pi_missing_data) + ")", "BATCH BCD (p=" + string(pi_missing_data) + ")", "BSCA (p=1)", "BATCH BCD (p=1)")



function [obj_value, timeValue] = batch_bcd(P, Q, A, K, R, Y, omega_t, omega_l, lambdastar, lambda1, mu_soft)

    T = size(Y,2);
    L = size(P,1);
    F = size(R,2);
    rho = size(Q,2);

    tic;
    obj_value(1) = getObj(Y,P,Q,R,A, lambdastar, lambda1)
    timeValue(1) = toc;

    for k = 1:K
        % update the anomaly map
        for f = 1:F
            ys = [];
            for t = 1:T
                if f == 1
                    % sum 2:
                    sum2 = 0;
                    for f_s = (f+1):F
                        sum2 = sum2 + R(:,f_s)*A(f_s, t);
                    end
                    % hole expression:
                    ys(:,t) = omega_t(:,:,t)*(Y(:,t) -  P*Q(t,:)' - sum2);
                else
                    % sum 1:
                    sum1 = 0;
                    for f_s = 1:(f-1)
                        sum1 = sum1 + R(:,f_s)*A_new(f_s, t);
                    end
                    % sum 2:
                    sum2 = 0;
                    for f_s = (f+1):F
                        sum2 = sum2 + R(:,f_s)*A(f_s, t);
                    end
                    % hole expression:
                    ys(:,t) = omega_t(:,:,t)*(Y(:,t) -  P*Q(t,:)' - sum1 - sum2);
                end
            end
    
            for t = 1:T
                A_new(f, t)  = sign(R(:,f)'*ys(:,t))*max(0, abs(R(:,f)'*ys(:,t)) - mu_soft) / norm(R(:,f),2);
                if isnan(A_new(f,t))
                    error("nan values in A. maybe 0 columns in R?")
                end
            end
            A_new(f,:);
        end
        A = A_new;
    
        % update the nominal traffic subspace:
        for l = 1:L
            P(l,:) = inv(lambdastar*eye(rho) + Q'*omega_l(:,:,l)*Q) * Q'*omega_l(:,:,l)*(Y(l,:)' - A'*R(l,:)');
        end
    
        % update the projection coefficients
        for t = 1:T
            Q(t,:) = inv(lambdastar*eye(rho) + P'*omega_t(:,:,t)*P)*P'*omega_t(:,:,t)*(Y(:,t) - R*A(:,t));
        end
    
        elapsedTime = toc;
        obj_value(k+1) = getObj(Y,P,Q,R,A, lambdastar, lambda1)%0.5*norm(Y-P*Q'-R*A,'fro').^2 + lambdastar/2*(norm(P,'fro').^2 + norm(Q,'fro').^2) + lambda1*norm(A,1) 
        timeValue(k+1) = elapsedTime + timeValue(k); 
    end

end




function [obj_value, timeValue] = bsca_missing_data(P, Q, A, K, R, Y, omega_t, omega_l, lambdastar, lambda1, mu_soft)

    T = size(Y,2);
    L = size(P,1);
    rho = size(Q,2);

    tic;
    obj_value(1) = getObj(Y,P,Q,R,A, lambdastar, lambda1)
    timeValue(1) = toc;

    for k = 1:K
    
        tic;

        % update the anomaly map
        Q = Q';
        % START:
        %------------------------
        for t = 1:T
            P_snake_t = omega_t(:,:,t)*P*Q(:,t);
            D_snake_t = omega_t(:,:,t)*R;
            Y_snake_t = omega_t(:,:,t)*Y(:,t);
        
            b = -P_snake_t - Y_snake_t;
            % A = D_snake_t
            % x = A(:,t)
        
            % formula from yang
            ATA = diag(D_snake_t'*D_snake_t);
            ATA(ATA==0) = 10^-10; % will be removed bny soft operator anyway
            
            r = ATA.*A(:,t) - D_snake_t'*(D_snake_t*A(:,t) - b);

            Bx = -(ATA.^-1).*soft_operator(r,mu_soft);

            if (sum(Bx) == 0)
                warning("Bx sum is zero")
            end
    
            gamma = min(1,max(0,-(D_snake_t*A(:,t) - b)'*D_snake_t*(Bx - A(:,t)) + mu_soft*(norm(Bx, 1) - norm(A(:,t),1)) / ( (D_snake_t*(Bx - A(:,t)))' * (D_snake_t*(Bx - A(:,t))) )));
    
            A_new(:,t) = A(:,t) + gamma*(Bx -A(:,t));
        end
    
        A = A_new;
        
        % END
        Q = Q';
        %------------------------
    
    
        % update the nominal traffic subspace:
        for l = 1:L
            P(l,:) = inv(lambdastar*eye(rho) + Q'*omega_l(:,:,l)*Q) * Q'*omega_l(:,:,l)*(Y(l,:)' - A'*R(l,:)');
        end
    
        % update the projection coefficients
        for t = 1:T
            Q(t,:) = inv(lambdastar*eye(rho) + P'*omega_t(:,:,t)*P)*P'*omega_t(:,:,t)*(Y(:,t) - R*A(:,t));
        end
    
        elapsedTime = toc;
        obj_value(k+1) = getObj(Y,P,Q,R,A, lambdastar, lambda1)
        timeValue(k+1) = elapsedTime + timeValue(k);
    end
end



% how should we measure?
function [obj_value] = getObj(Y,P,Q,R,A, lambdastar, lambda1)
    %obj_value = 0.5*norm(Y-P*Q'-R*A,'fro').^2 + lambdastar/2*(norm(P,'fro').^2 + norm(Q,'fro').^2) + lambda1*norm(A,1);
    obj_value = 0.5*norm(Y-P*Q'-R*A,'fro').^2;
end



function result = soft_operator(x_in, a)
    result = max(x_in - a*ones(size(x_in)), zeros(size(x_in))) - max(-x_in - a*ones(size(x_in)), zeros(size(x_in)));
end

