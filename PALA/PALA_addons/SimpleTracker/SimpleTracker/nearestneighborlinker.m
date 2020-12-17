function [ target_indices target_distances unassigned_targets ] = nearestneighborlinker(source, target, max_distance)
%NEARESTNEIGHBORLINKER link two lists of points based on nearest neighbor.
%
% target_indices = NEARESTNEIGHBORLINKER(source, target) finds for each
% point in 'source' the closest point in 'target'. These 2 inputs must be
% arrays with one point per row, and have their cartesian coordinates in
% each column (1D, 2D, 3D, ...). Nearest neighbor matching is based on
% euclidean distance. The two arrays might not have the same number of
% points.
%
% The indices of the 'target' points are returned in an array
% 'target_indices', so that each row in 'source' matches the corresponding
% row in 'target(target_indices, :)'.
%
% The linking is exclusive: one source point is linked to at most one
% target point, and conversely. The linking is only locally optimal: the
% two closest points amongst the two sets are sought for first, then the
% second closest pair, excluding the first, etc... This ensures that the
% resulting linking will not depend on the order of the points in each set.
%
% target_indices = NEARESTNEIGHBORLINKER(source, target, max_distance) adds
% a condition on distance. If the nearest neighbor is found to be at a
% distance larger than the given 'max_distance', they are not linked, and
% the 'target_indices' receive the value -1 for this source point. The same
% happens if all target points are exhausted.
% 
% [ target_indices target_distances ] = 
%                                   NEARESTNEIGHBORLINKER(source, target)
% additionaly return the distance to the matched target point. Un-matched
% source points have a distance value set to NaN.
%
% [ target_indices target_distances unmatched_targets ]= 
%                                   NEARESTNEIGHBORLINKER(source, target)
% additionaly return the indices of the points in 'target' that have not
% been linked.
%
% This is the cheapest (in term of accuracy) algorithm for linking that can
% be made. In particular, it is not guaranteed (and it is generally not the
% case) that the returned linking is an optimum for the sum of distances.
% Each source point is matched regardless of the others, there is no global
% optimization here (the Hungarian algorithm does that). Also, there exists
% refinement to nearest neighbor searches, such as the use of KD-trees;
% this contribution is exempt of such developments.
%
% EXAMPLE:
% 
% n_points = 20;
% source = 10 * rand(n_points, 2);
% target = source + rand(n_points, 2);
% target_indices = nearestneighborlinker(source, target);
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
% VERSION HISTORY
%
% * v1.0 - November 2011 - Initial release.
% * v1.1 - May 2012 - Fix a severe bug thanks to Dave Cade
%
% Jean-Yves Tinevez < jeanyves.tinevez@gmail.com> November 2011 - 2012

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
    
    target_indices = -1 * ones(n_source_points, 1);
    target_distances = NaN(n_source_points, 1);
    
    % Parse distance matrix
    while ~all(isinf(D(:)))
        
        [ min_D closest_targets ] = min(D, [], 2); % index of the closest target for each source points
        [ ~, sorted_index ] = sort(min_D);
        
        for i = 1 : numel(sorted_index)
            
            source_index =  sorted_index(i);
            target_index =  closest_targets ( sorted_index(i) );
            
            % Did we already assigned this target to a source?
            if any ( target_index == target_indices )
                
                % Yes, then exit the loop and change the distance matrix to
                % prevent this assignment
                break
                
            else
                
                % No, then store this assignment
                target_indices( source_index ) = target_index;
                target_distances ( source_index ) = sqrt ( min_D (  sorted_index(i) ) );
                
                % And make it impossible to find it again by putting the target
                % point to infinity in the distance matrix
                D(:, target_index) = Inf;
                % And the same for the source line
                D(source_index, :) = Inf;
                
                if all(isinf(D(:)))
                    break
                end
                
            end
            
        end
        
    end
    
    unassigned_targets = setdiff ( 1 : n_target_points , target_indices );
        
end