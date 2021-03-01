function [path, totalCost, farthestPreviousHop, farthestNextHop] = ...
	 dijkstra(n, netCostMatrix, s, d, farthestPreviousHop, farthestNextHop)
 %
%  SYNTAX:  [path, totalCost] = dijkstra(n, CostMatrix, s, d)
%
%  PURPOSE:  Find the least cost path between two nodes in a network. 
%
%  INPUT:
%    n = the number of nodes in the network
%    CostMatrix = matrix whose (i,j)th entry defines the cost of moving from
%	node i to node j or conversely.  Use Inf if impossible.
%    s = source node index
%    d = destination node index
%
%  OUTPUT:
%    path = the list of nodes in the path from source to destination
%    totalCost = the total cost of the path
% 
%  RESTRICTIONS: Node numbers must be consecutive, starting at 1.
 
%  METHOD:  Uses recursive call to itself, passing the variables 
%	farthestPreviousHop and farthestNextHop.  I don't know what these do.
%
%  WRITTEN BY:  Mike Walker 	ON 	4/25/07
%	(based on algorithm downloaded from the web, no author in the file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)dijkstra.m	1.2 02/19/09

if(nargin < 5)
   for i = 1:n,
      % initialize the farthest node to be itself;
      farthestPreviousHop(i) = i;     % used to compute the RTS/CTS range;
      farthestNextHop(i) = i;
   end;
end;

% all the nodes are un-visited;
visited(1:n) = 0;

distance(1:n) = inf;    % it stores the shortest distance between each node and the source node;
parent(1:n) = 0;

distance(s) = 0;
for i = 1:(n-1),
    temp = [];
    for h = 1:n,
         if visited(h) == 0   % in the tree;
             temp=[temp distance(h)];
         else
             temp=[temp inf];
         end
     end;
     [t, u] = min(temp);    % it starts from node with the shortest distance to the source;
     visited(u) = 1;       % mark it as visited;
     for v = 1:n,           % for each neighbors of node u;
         if ( ( netCostMatrix(u, v) + distance(u)) < distance(v) )
             distance(v) = distance(u) + netCostMatrix(u, v);   % update the shortest distance when a shorter path is found;
             parent(v) = u;                                     % update its parent;
         end;             
     end;
end;

path = [];
if parent(d) ~= 0   % if there is a path!
    t = d;
    path = [d];
    while t ~= s
        p = parent(t);
        path = [p path];
        
        if netCostMatrix(t, farthestPreviousHop(t)) < netCostMatrix(t, p)
            farthestPreviousHop(t) = p;
        end;
        if netCostMatrix(p, farthestNextHop(p)) < netCostMatrix(p, t)
            farthestNextHop(p) = t;
        end;

        t = p;      
    end;
end;

totalCost = distance(d);

return;
