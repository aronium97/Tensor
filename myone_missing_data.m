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
p = 0.005;

% make v_l_t
V = randn(L,T).*sigma^2;
w = randn(r,T);
U = randn(F,r)*1/F;
Z = U*w;
helpDistri = rand(F,T);
A = (helpDistri<(p/2))*(-1) + (helpDistri>=(p/2) & helpDistri<p)*(1) + (helpDistri>=p)*0;

omega_t = eye(L);
omega_l = eye(T);

R = getR(L,F,nodes);
%scatter(nodes(:,1),nodes(:,2))

Y = omega_t*(R*Z + R*A + V);

K = 200; % num iterations

lambda1 = 10%max(max(abs(R'*V)));
lambdastar = 10%(sqrt(T) + sqrt(F)*sqrt(pi))*sigma%norm(V,1);

mu_soft = 200;

% init P and Q at random
% X = LxT = PQ'
Q = randn(T,rho); 
P = randn(L,rho);%5*R*A*Q;%randn(L,rho);






obj_value = bsca_missing_data(P, Q, A, K, R, Y, omega_t, omega_l, lambdastar, lambda1, mu_soft);

plot(obj_value, "--")
set(gca, 'YScale', 'log')


function obj_value = bsca_missing_data(P, Q, A, K, R, Y, omega_t, omega_l, lambdastar, lambda1, mu_soft)

    T = size(Y,2);
    L = size(P,1);
    rho = size(Q,2);


    obj_value(1) = getObj(Y,P,Q,R,A, lambdastar, lambda1)

    for k = 1:K
        % update the anomaly map
        Q = Q';
        % START:
        %------------------------
        for t = 1:T
            P_snake_t = omega_t*P*Q(:,t);
            D_snake_t = omega_t*R;
            Y_snake_t = omega_t*Y(:,t);
        
            b = -P_snake_t - Y_snake_t;
            % A = D_snake_t
            % x = A(:,t)
        
            % formula from yang
            ATA = diag(D_snake_t'*D_snake_t);
            
            r = ATA.*A(:,t) - D_snake_t'*(D_snake_t*A(:,t) - b);
    
            mu_soft = 100;

            Bx = (ATA.^-1).*soft_operator(r,mu_soft);
    
            gamma = min(1,max(0,-(D_snake_t*A(:,t) - b)'*D_snake_t*(Bx - A(:,t)) + mu_soft*(norm(Bx, 1) - norm(A(:,t),1)) / ( (D_snake_t*(Bx - A(:,t)))' * (D_snake_t*(Bx - A(:,t))) )));
    
            A_new(:,t) = A(:,t) + gamma*(Bx -A(:,t));
        end
    
        A = A_new;
        
        % END
        Q = Q';
        %------------------------
    
    
        % update the nominal traffic subspace:
        for l = 1:L
            P(l,:) = inv(lambdastar*eye(rho) + Q'*omega_l*Q) * Q'*omega_l*(Y(l,:)' - A'*R(l,:)');
        end
    
        % update the projection coefficients
        for t = 1:T
            Q(t,:) = inv(lambdastar*eye(rho) + P'*omega_t*P)*P'*omega_t*(Y(:,t) - R*A(:,t));
        end
    
        obj_value(k+1) = getObj(Y,P,Q,R,A, lambdastar, lambda1)
    end
end





function [obj_value] = getObj(Y,P,Q,R,A, lambdastar, lambda1)
    obj_value = 0.5*norm(Y-P*Q'-R*A,'fro').^2 + lambdastar/2*(norm(P,'fro').^2 + norm(Q,'fro').^2) + lambda1*norm(A,1);
end

function result = soft_operator(x_in, a)
    result = max(x_in - a*ones(size(x_in)), zeros(size(x_in))) - max(-x_in - a*ones(size(x_in)), zeros(size(x_in)));
end

