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
scatter(nodes(:,1),nodes(:,2))

Y = omega_t*(R*Z + R*A + V);

K = 70; % num iterations

lambda1 = 10%max(max(abs(R'*V)));
lambdastar = 10%(sqrt(T) + sqrt(F)*sqrt(pi))*sigma%norm(V,1);

% init P and Q at random
% X = LxT = PQ'
Q = randn(T,rho); 
P = randn(L,rho);%5*R*A*Q;%randn(L,rho);



obj_value(1) = getObj(Y,P,Q,R,A, lambdastar, lambda1)%0.5*norm(Y-P*Q'-R*A,'fro').^2 + lambdastar/2*(norm(P,'fro').^2 + norm(Q,'fro').^2) + lambda1*norm(A,1)




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
                ys(:,t) = omega_t*(Y(:,t) -  P*Q(t,:)' - sum2);
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
                ys(:,t) = omega_t*(Y(:,t) -  P*Q(t,:)' - sum1 - sum2);
            end
        end

        for t = 1:T
            A_new(f, t)  = sign(R(:,f)'*ys(:,t))*max(0, abs(R(:,f)'*ys(:,t)) - lambda1) / norm(R(:,f),2);
            if isnan(A_new(f,t))
                error("nan values in A. maybe 0 columns in R?")
            end
        end
        A_new(f,:);
    end
    A = A_new;

    % update the nominal traffic subspace:
    for l = 1:L
        P(l,:) = inv(lambdastar*eye(rho) + Q'*omega_l*Q) * Q'*omega_l*(Y(l,:)' - A'*R(l,:)');
    end

    % update the projection coefficients
    for t = 1:T
        Q(t,:) = inv(lambdastar*eye(rho) + P'*omega_t*P)*P'*omega_t*(Y(:,t) - R*A(:,t));
    end

    obj_value(k+1) = getObj(Y,P,Q,R,A, lambdastar, lambda1)%0.5*norm(Y-P*Q'-R*A,'fro').^2 + lambdastar/2*(norm(P,'fro').^2 + norm(Q,'fro').^2) + lambda1*norm(A,1)
end

plot(obj_value, "--")
set(gca, 'YScale', 'log')



function [obj_value] = getObj(Y,P,Q,R,A, lambdastar, lambda1)
    obj_value = 0.5*norm(Y-P*Q'-R*A,'fro').^2 + lambdastar/2*(norm(P,'fro').^2 + norm(Q,'fro').^2) + lambda1*norm(A,1);
end

