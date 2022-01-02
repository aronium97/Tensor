clc
close all
clear all

L = 48;
F = 210;
T = 100;

%lambda1 = 10;
%lambdastar = 10;
rho = 5;

% make random variables
sigma = 10^-9;
r = 2;
p = 0.005;

% make v_l_t
V = randn(L,T).*sigma^2;
w = randn(r,T);
U = randn(F,r)*1/F;
Z = U*w;
helpDistri = rand(F,T);
A = (helpDistri<(p/2))*(-1) + (helpDistri>=(p/2) & helpDistri<p/2)*(1) + (helpDistri>=p/2)*0;

omega_t = eye(L);
omega_l = eye(T);

R = zeros(L,F);

% node matrix
nodes = [[0,0.3];[0.18,0.28];[0.19, 0.05];[0.41,0.39];[0.6,0.21];[1,0.21];[0.91,0.4];[0.8,0.4];[0.39,0.51];[0.59,0.59];[0.6,0.6];[0.42,0.79];[0.79, 0.9];[0.92,0.98];[0.95,0.85]];
noOfLinks = 0;
links = [];
% check for no of links
for nodenumber = 1:length(nodes)-1
    xyStart = nodes(nodenumber,:);
    for nodenumberTarget = (nodenumber+1):length(nodes)
        xyTarget = nodes(nodenumberTarget,:);
        if norm(xyStart-xyTarget,2) <= 0.35
            noOfLinks = noOfLinks+1;
            if length(links)==0
                links =[[[nodenumber],[nodenumberTarget]]];
            else
                links =[ links; [[nodenumber],[nodenumberTarget]]];  
            end
        end
    end
end
% mnake it bidirektional
links = [links;links(:,2),links(:,1)];
noOfLinks = length(links);
if not(noOfLinks == L)
    print("error: noOfLinks")
    exit()
end

% check for flows
routes = {};
flownumber = 0; % f
for nodenumber = 1:length(nodes)
    xyStart = nodes(nodenumber,:);
    for nodenumberTarget = 1:length(nodes)
        if not(nodenumberTarget==nodenumber)
            flownumber = flownumber + 1;
            xyOurs = xyStart;
            xyTarget = nodes(nodenumberTarget,:);
            % we have a flow now. now lets get the shortest path
       
            % find the best route
            nodesOfRoute = [nodenumber];
            bestNeighboor = nodenumber;
            unallowedTransitions = ["blub"];
            running = true;
            while not(bestNeighboor == nodenumberTarget) & running
                minDis = 20;
                for potentialNeighboorNode = 1:length(nodes)
                    xyNeighboor = nodes(potentialNeighboorNode,:);
                    potentialTransition = string(nodesOfRoute(end)) + string(potentialNeighboorNode);
                    if not(ismember(potentialNeighboorNode, nodesOfRoute)) & not(ismember(potentialTransition, unallowedTransitions)) & norm(xyOurs-xyNeighboor,2) <= 0.35 & norm(xyTarget-xyNeighboor,2) < minDis
                        minDis = norm(xyTarget-xyNeighboor,2);
                        bestNeighboor = potentialNeighboorNode;
                    end  
                end
                if not(nodesOfRoute(end) == bestNeighboor)                     
                    nodesOfRoute = [nodesOfRoute, bestNeighboor];
                    xyOurs = nodes(nodesOfRoute(end),:);
                    
                    % check for link pair and activate r
                    linkPair = [nodesOfRoute(end-1),nodesOfRoute(end)];
                    linknumber = find(sum(links == linkPair,2)==2); % l
                    R(linknumber,flownumber) = 1;
                else
                    if length(nodesOfRoute) ==1
                        running = false;
                        nodenumber
                        nodenumberTarget
                        "-----"
                    else

                        lastTransition = string(nodesOfRoute(end-1)) + string(nodesOfRoute(end));
                        unallowedTransitions = [unallowedTransitions; lastTransition];
    
                        % check for link pair and deactivate r
                        linkPair = [nodesOfRoute(end-1),nodesOfRoute(end)];
                        linknumber = find(sum(links == linkPair,2)==2); % l
                        R(linknumber,flownumber) = 0;
    
                        % delete the route
                        nodesOfRoute = nodesOfRoute(1:end-1);
                        xyOurs = nodes(nodesOfRoute(end),:);
                        bestNeighboor = nodesOfRoute(end);
                    end
                    
                    
                end
            end
        end

    end
end



Y = omega_t*(R*Z + R*A + V);

K = 12; % num iterations

lambda1 = 10%max(max(abs(R'*V)));
lambdastar = 10%(sqrt(T) + sqrt(F)*sqrt(pi))*sigma%norm(V,1);

% init P and Q at random
% X = LxT = PQ'
Q = randn(T,rho); 
P = randn(L,rho);%5*R*A*Q;%randn(L,rho);

obj_value(1) = 0.5*norm(Y-P*Q'-R*A,'fro').^2 + lambdastar/2*(norm(P,'fro').^2 + norm(Q,'fro').^2) + lambda1*norm(A,1)

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
            A_new(f, t) = sign(R(:,f)'*ys(:,t))*max(0, abs(R(:,f)'*ys(:,t)) - lambda1) / norm(R(:,f),2);
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

    obj_value(k+1) = 0.5*norm(Y-P*Q'-R*A,'fro').^2 + lambdastar/2*(norm(P,'fro').^2 + norm(Q,'fro').^2) + lambda1*norm(A,1)
end

plot(obj_value)
set(gca, 'YScale', 'log')

