clear; clc;

format long;

%test

% parameters
N = 200;
K = 400;
I = 400; % X: N*K, D: N*I, S: I*K
rho_real = 5;
rho = 3; % rank of X; P: N*rho, Q: rho*K

% number of samples in Monte Carlo simulations
Sample = 20;

% generate the data
D = randn(N, I);
for n = 1: 1: N
    D(n, :) = D(n, :) / norm(D(n, :));
end

S0 = sprandn(I, K, 0.05); % density   

% make real data
P0 = sqrt(100/I) * randn(N, rho_real);
Q0 = sqrt(100/K) * randn(rho_real, K);
X0 = P0 * Q0; % perfect X
sigma = 0.01;
V = sigma * randn(N, K); % noise

Y = X0 + D * S0 + V; % observation

% own parameters
c_lambda = 2.5 * 10^-1;
lambda = c_lambda * norm(Y); 
c_mu = 2 * 10^-3; 
mu = c_mu / 10 * norm(D'*(Y), inf);

% initial point (common for all algorithms)
initial_P = randn(N, rho);
initial_Q = randn(rho, K);
initial_S = zeros(I,K);


MaxIter_bSCA = 50;

result = FUN_admm(initial_P, initial_Q, initial_S, MaxIter_bSCA, D, Y, lambda, mu, rho)
plot(result)


hold on
result1 = bsca(initial_P, initial_Q, initial_S, MaxIter_bSCA, D, Y, lambda, mu)
plot(result1)


result2 = sca(initial_P, initial_Q, initial_S, MaxIter_bSCA, D, Y, lambda, mu)
plot(result2)

legend("admm","bsca","sca")

function result = objective_function(Y, P, Q, D, S, lambda, mu)
    result = 0.5 * norm(Y - P * Q - D * S, 'fro') ^ 2 + 0.5 * lambda * (norm(P, 'fro') ^ 2 + norm(Q, 'fro') ^ 2) + mu * norm(S(:), 1);
end

function result = soft_operator(x_in, a)
    result = max(x_in - a*ones(size(x_in)), zeros(size(x_in))) - max(-x_in - a*ones(size(x_in)), zeros(size(x_in)));
end

function obj_value = bsca(initial_P, initial_Q, initial_S, MaxIter_bSCA, D, Y, lambda, mu)

    P = initial_P; 
    Q = initial_Q; 
    S = initial_S;

    diagDD = diag(diag(D'*D));
    invdiagDD = inv(diagDD);
    obj_value = [];
    for i = 1:MaxIter_bSCA
        % P
        bp_z = (Y - D*S)*Q'*inv(Q*Q' + lambda*eye(size(Q*Q')));
        P = bp_z;

        % Q
        bq_z = inv(P'*P + lambda*eye(size(P'*P)))*P'*(Y-D*S);
        Q = bq_z;

        % S
        bs_z = invdiagDD * soft_operator(diagDD*S - D'*(D*S - Y + P*Q), mu);
        gamma = min(max(-(trace((P*Q + D*S - Y)'*D*(bs_z - S)) + mu*(norm(bs_z,1) - norm(S,1))) / norm(D*(bs_z - S),'fro')^2, 0),1);
        S_new = S + gamma*(bs_z - S);
        S = S_new;
        disp("iteration done")

        % save value
        s(i) = objective_function(Y, P, Q, D, S, lambda, mu);
    end   

    disp('hi')
    
    obj_value = s;
end




function s = sca(initial_P, initial_Q, initial_S, MaxIter_bSCA, D, Y, lambda, mu)

    P = initial_P; 
    Q = initial_Q; 
    S = initial_S;

    diagDD = diag(diag(D'*D));
    invdiagDD = inv(diagDD);
    obj_value = [];
    for i = 1:MaxIter_bSCA
        % P
        bp_z = (Y - D*S)*Q'*inv(Q*Q' + lambda*eye(size(Q*Q')));
        deltaP = bp_z - P;

        % Q
        bq_z = inv(P'*P + lambda*eye(size(P'*P)))*P'*(Y-D*S);
        deltaQ = bq_z - Q;

        % S
        bs_z = invdiagDD * soft_operator(diagDD*S - D'*(D*S - Y + P*Q), mu);
        deltaS = bs_z - S;


        % Step size
        A = deltaP * deltaQ;
        B = P * deltaQ + deltaP * Q + D * deltaS;
        C = P * Q + D * S - Y;
        
        a = 2 * sum(sum(A.^2, 1));
        b = 3 * sum(sum(A.*B, 1));
        c = sum(sum(B.^2, 1)) + 2 * sum(sum(A.*C, 1)) + lambda * sum(sum(deltaP.^2, 1)) + lambda * sum(sum(deltaQ.^2, 1));
        d = sum(sum(B.*C, 1)) + lambda * sum(sum(deltaP.*P, 1)) + lambda * sum(sum(deltaQ.*Q, 1)) + mu * (norm(bs_z(:), 1) - norm(S(:), 1));

        % calculating the stepsize by closed-form expression
        Sigma1 = (-(b / 3 / a) ^ 3 + b * c / 6 / a^2 - d / 2 / a);
        Sigma2 = c / 3 / a - (b / 3 / a) ^ 2;
        Sigma3 = Sigma1 ^ 2 + Sigma2 ^ 3;
        Sigma3_sqrt = sqrt(Sigma3);
        if Sigma3 >= 0
            gamma = nthroot(Sigma1 + Sigma3_sqrt, 3)...
                + nthroot(Sigma1 - Sigma3_sqrt, 3)...
                - b / 3 / a;
        else
            C1 = 1; C1(4) = - (Sigma1 + Sigma3_sqrt);
            C2 = 1; C2(4) = - (Sigma1 - Sigma3_sqrt);
            R = real(roots(C1) + roots(C2)) - b/3/a * ones(3,1);
            gamma = min(R(R>0));
            clear C1 C2 R;
        end
        

        % update blocks
        P = P + gamma * deltaP; 
        Q = Q + gamma * deltaQ; 
        S = S + gamma * deltaS;

        disp("iteration done")

        % save value
        s(i) = objective_function(Y, P, Q, D, S, lambda, mu);
    end   

    disp('hi')
    
    obj_value = s;
end







function result = FUN_admm(initial_P, initial_Q, initial_S, MaxIter_admm, D, Y, lambda, mu, rho)
    % ADMM algorithm
    %global val0;
    %global Y;
    %gglobal D;
    global lambda;
    global mu;
    %global rho;
    %global I;
    %global K;

    L  = initial_P; 
    A  = initial_S;
    M = zeros(size(A));
    B = M;
    R = D;

    c  = 10^4;

    result = [];
    for t = 1:1:MaxIter_admm
        
        % update local variables
        M = M + mu*(B + A);

        % update first group of local primal variables
        Q = (Y'*L - B'*R'*L)*inv(L'*L + (1+lambda*eye(rho)));

        A = inv(c*3)*soft_operator(M + c*B, mu/c);
        
        % update second group of local primal variables
        L = (Y - R*B)*Q*inv(Q'*Q + lambda*eye(rho));

        % update auxiliary local primal variables:
        B = inv(R'*R + c*eye(size(A,1))) * (R'*(Y - L*Q') - M + c*A);

        disp('iteration admm done')
        Q_ = Q';
        result(t) = objective_function(Y, L, Q_, R, A, lambda, mu);

        
    end

    X_admm = L * Q';
    S_admm = A;

    disp(result)

    
end




