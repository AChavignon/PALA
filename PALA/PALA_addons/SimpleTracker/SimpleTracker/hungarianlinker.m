function [ target_indices target_distances unassigned_targets total_cost ] = hungarianlinker(source, target, max_distance)
%HUNGARIANLINKER link two lists of points based on the hungarian algorithm.
%
% target_indices = HUNGARIANLINKER(source, target) finds for each point in
% 'source' the closest point in 'target'. These 2 inputs must be arrays
% with one point per row, and have their cartesian coordinates in each
% column (1D, 2D, 3D, ...). Source to target assignment is based on the
% famous hungarian algorithm using its excellent implementation by the
% excellent Yi Cao. The two arrays might not have the same number of
% points.
%
% The indices of the 'target' points are returned in an array
% 'target_indices', so that each row in 'source' matches the corresponding
% row in 'target(target_indices, :)'.
%
% The linking is exclusive: one source point is linked to at most one
% target point, and conversely. The linking is globally optimal: the sum of
% the square distance is minimized, contrary to the naive nearest neighbor
% approach.
%
% target_indices = HUNGARIANLINKER(source, target, max_distance) adds
% a condition on distance. If the nearest neighbor is found to be at a
% distance larger than the given 'max_distance', they are not linked, and
% the 'target_indices' receive the value -1 for this source point. The same
% happens if all target points are exhausted.
% 
% [ target_indices target_distances ] = HUNGARIANLINKER(source, target)
% additionaly return the distance to the matched target point. Un-matched
% source points have a distance value set to NaN.
%
% [ target_indices target_distances unmatched_targets ] =
%                                         HUNGARIANLINKER(source, target) 
% additionaly return the indices of the points in 'target' that have not
% been linked.
%
% [ target_indices target_distances unmatched_targets total_cost ] =
%                                         HUNGARIANLINKER(source, target) 
% additionaly return the globally optimized value of the square distance
% sum.
%
% The matching algorithm used here is one of the best available and ensures
% that the resulting assignment is a optimum. However the price to pay is
% an increased complexity. The cost for the naive nearest neighbor approach
% roughly scales as O(p^2) where p is the number of source points. The
% munkres implementation of the hungarian algorithm by Yi Cao is in O(p^3),
% and is the best so far.
%
% EXAMPLE:
% 
% n_points = 20;
% source = 10 * rand(n_points, 2);
% target = source + rand(n_points, 2);
% target_indices = hungarianlinker(source, target);
% colors = hsv(n_points);
% figure
% hold on
% for i = 1 :n_points
%    plot(source(i,1), source(i,2), 'o', 'Color', colors(i,:))
%    plot(target(target_indices(i),1), target(target_indices(i),2), 's', ...
%       'Color', colors(i,:))
%    plot( [ source(i,1) target(target_indices(i),1) ] , ...
%       [ source(i,2)  target(target_indices(i),2) ], ...
%        'Color', colors(i,:))
% end
% 
%
% Jean-Yves Tinevez <jeanyves.tinevez@gmail.com>.
% However all credits should go to Yi Cao, which did the hard job of
% implementing the Munkres algorithm; this file is merely a wrapper for it.

    if nargin < 3
        max_distance = Inf;
    end

    n_source_points = size(source, 1);
    n_target_points = size(target, 1);
    
    D = NaN(n_source_points, n_target_points);
    
    % Build distance matrix
    for i = 1 : n_source_points
        
        % Pick one source point
        current_point = source(i, :);
        
        % Compute square distance to all target points
        diff_coords = target - repmat(current_point, n_target_points, 1);
        square_dist = sum(diff_coords.^2, 2);
        
        % Store them
        D(i, :) = square_dist;
        
    end
    
    % Deal with maximal linking distance: we simply mark these links as already
    % treated, so that they can never generate a link.
    D ( D > max_distance * max_distance ) = Inf;
    
    % Find the optimal assignment is simple as calling Yi Cao excellent FEX
    % submission.
    [ target_indices total_cost ] = munkres(D);
    % Set unmatched sources to -1
    target_indices ( target_indices  == 0 ) = -1;
    
    % Collect distances
    target_distances = NaN(numel(target_indices), 1);
    for i = 1 : numel(target_indices)
        if target_indices(i) < 0
            continue
        end
        
        target_distances(i) = sqrt ( D ( i , target_indices(i)) );
        
    end
    
    unassigned_targets = setdiff ( 1 : n_target_points , target_indices );
    
    
end