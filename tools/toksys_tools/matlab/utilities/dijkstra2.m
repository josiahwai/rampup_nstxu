function [dist, path] = dijkstra(nodes, segments, start_id, finish_id)
%DIJKSTRA Calculates the shortest path between two points on
% a map using Dijkstra's Algorithm
% 
% DIJKSTRA requires a list of nodes and a list of segments as inputs,
% as well as a starting node. (These are created randomly if they are not
% specified as inputs into the function
% 
% Optionally, an ending node can also be used as an input to terminate the
% algorithm early (otherwise, the algorithm will find the shortest paths
% and distances to all of the nodes on the map from the starting node)
%
% AUTHOR: Joseph Kirk (c) 11/2006
% EMAIL: jdkirk630 at gmail dot com
%
% For more information about Dijkstra's Algorithm check out:
%       http://en.wikipedia.org/wiki/Dijkstra's_algorithm or
%       http://mathworld.wolfram.com/DijkstrasAlgorithm.html

if (nargin < 3) % SETUP
    % (GENERATE NODES AND SEGMENTS IF NOT GIVEN AS INPUTS)
    % Create a random set of nodes/vertices, and connect some of them with
    % edges/segments. Then graph the resulting map.
    num_nodes = 40; L = 100; max_seg_length = 30; ids = (1:num_nodes)';
    nodes = [ids L*rand(num_nodes, 2)]; % create random nodes
    h = figure(1); plot(nodes(:, 2), nodes(:, 3), 'k.') % plot the nodes
    text(nodes(num_nodes, 2), nodes(num_nodes, 3), ...
        [' ' num2str(ids(num_nodes))], 'Color', 'b', 'FontWeight', 'b')
    hold on
    c = 1; segments = [];
    for i = 1:num_nodes-1 % create edges between some of the nodes
        text(nodes(i, 2), nodes(i, 3), ...
            [' ' num2str(ids(i))], 'Color', 'b', 'FontWeight', 'b')
        for j = i+1:num_nodes
            d = sqrt(sum((nodes(i, 2:end) - nodes(j, 2:end)).^2));
            if and(d < max_seg_length, rand < 0.6)
                plot([nodes(i, 2) nodes(j, 2)], [nodes(i, 3) nodes(j, 3)], 'k.-')
                % add this link to the segments list
                segments = [segments; c nodes(i, 1) nodes(j, 1)];
                c = c + 1;
            end
        end
    end
    axis([0 L 0 L])
    % Calculate Shortest Path Using Dijkstra's Algorithm
    % Get random starting/ending nodes, compute the shortest distance and path.
    start_id = ceil(num_nodes*rand)
    finish_id = ceil(num_nodes*rand)
    [distance, path] = dijkstra(nodes, segments, start_id, finish_id)
    % If a Shortest Path exists, Plot it on the Map.
    figure(h)
    for k = 2:length(path)
        m = find(nodes(:, 1) == path(k-1));
        n = find(nodes(:, 1) == path(k));
        plot([nodes(m, 2) nodes(n, 2)], [nodes(m, 3) nodes(n, 3)], 'ro-')
    end
    title(['Shortest Distance from ' num2str(start_id) ' to ' ...
        num2str(finish_id) ' = ' num2str(distance)])
    hold off

else %--------------------------------------------------------------------------
    % MAIN FUNCTION - DIJKSTRA'S ALGORITHM

    % initializations
    node_ids = nodes(:, 1);
    [num_map_pts, cols] = size(nodes);
    table = sparse(num_map_pts, 2);
    shortest_distance = Inf(num_map_pts, 1);
    settled = zeros(num_map_pts, 1);
    path = num2cell(NaN(num_map_pts, 1));
    col = 2;
    pidx = find(start_id == node_ids);
    shortest_distance(pidx) = 0;
    table(pidx, col) = 0;
    settled(pidx) = 1;
    path(pidx) = {start_id};
    if (nargin < 4) % compute shortest path for all nodes
        while_cmd = 'sum(~settled) > 0';
    else % terminate algorithm early
        while_cmd = 'settled(zz) == 0';
        zz = find(finish_id == node_ids);
    end
    while eval(while_cmd)
        % update the table
        table(:, col-1) = table(:, col);
        table(pidx, col) = 0;
        % find neighboring nodes in the segments list
        jj = find(node_ids(pidx) == segments(:, 2));
        kk = find(node_ids(pidx) == segments(:, 3));
        neighbor_ids = [segments(jj, 3); segments(kk, 2)];
        % calculate the distances to the neighboring nodes and keep track of the paths
        for k = 1:length(neighbor_ids)
            cidx = find(neighbor_ids(k) == node_ids);
            if ~settled(cidx)
                d = sqrt(sum((nodes(pidx, 2:cols) - nodes(cidx, 2:cols)).^2));
                if (table(cidx, col-1) == 0) || ...
                        (table(cidx, col-1) > (table(pidx, col-1) + d))
                    table(cidx, col) = table(pidx, col-1) + d;
                    tmp_path = path(pidx);
                    path(cidx) = {[tmp_path{1} neighbor_ids(k)]};
                else
                    table(cidx, col) = table(cidx, col-1);
                end
            end
        end
        % find the minimum non-zero value in the table and save it
        nidx = find(table(:, col));
        ndx = find(table(nidx, col) == min(table(nidx, col)));
        if isempty(ndx)
            break
        else
            pidx = nidx(ndx(1));
            shortest_distance(pidx) = table(pidx, col);
            settled(pidx) = 1;
        end
    end
    if (nargin < 4) % return the distance and path arrays for all of the nodes
        dist = shortest_distance';
        path = path';
    else % return the distance and path for the ending node
        dist = shortest_distance(zz);
        path = path(zz);
        path = path{1};
    end
end

