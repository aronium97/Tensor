function [R] = getR(L,F, nodes)

R = zeros(L,F);

% node matrix
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
    print("noOfLinks")
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


end