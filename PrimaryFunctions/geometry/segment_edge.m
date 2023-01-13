function [new_points] = segment_edge(Points,unit_segment_length)
triangle_current = [1 2 3];
triangle_points = Points;
edges = [triangle_current(1),triangle_current(2);...
    triangle_current(2),triangle_current(3);...
    triangle_current(3),triangle_current(1)];
% Points_for_triangulation = triangle_points;
new_points = [];
for edg = 1:3
    edge_current = edges(edg,:);
    edge_current = sort(edge_current,2);
    edge_current_vect = triangle_points(triangle_current==edge_current(2),:) - triangle_points(triangle_current==edge_current(1),:);
    edge_points = [triangle_points(triangle_current==edge_current(2),:);triangle_points(triangle_current==edge_current(1),:)];
    unit_edge_current_vect = edge_current_vect./vecnorm(edge_current_vect,2,2);
%     segmented_noter_old = 0;
    
    %     if segmented_noter_old == 0
    %                     [common_TR_ID] = find_opposite_triangles_to_edge(obj,edge_current);
    %                     neighbor_current = obj_list(common_TR_ID(1));
    x_coor_current_edge = edge_points(:,1);
    y_coor_current_edge = edge_points(:,2);
    z_coor_current_edge = edge_points(:,3);
    num_points = floor(vecnorm(edge_current_vect,2,2)/unit_segment_length);
    if abs(num_points*unit_segment_length-vecnorm(edge_current_vect,2,2))<unit_segment_length/2
        num_points = num_points-1;
    end
    ind_multipliers = (1:num_points)*unit_segment_length;
    
    xmin = x_coor_current_edge(1); xmax = x_coor_current_edge(2);
    if xmin(1)~=xmax(1)
        x_new_points = xmin(1) - unit_edge_current_vect(1)*ind_multipliers;
    else
        x_new_points = xmin(1)*ones(1,(num_points));
    end
    
    ymin = y_coor_current_edge(1); ymax = y_coor_current_edge(2);
    if ymin(1)~=ymax(1)
        y_new_points = ymin(1) - unit_edge_current_vect(2)*ind_multipliers;
    else
        y_new_points = ymin(1)*ones(1,(num_points));
    end
    
    zmin = z_coor_current_edge(1); zmax = z_coor_current_edge(2);
    if zmin(1)~=zmax(1)
        z_new_points = zmin(1) - unit_edge_current_vect(3)*ind_multipliers;
    else
        z_new_points = zmin(1)*ones(1,(num_points));
    end
    current_generated_points = [x_new_points' y_new_points' z_new_points'];
    %         if edg == 1
    %             obj.EdgeSegmentationInfo.Edge1.Points = [x_new_points' y_new_points' z_new_points'];
    %         elseif edg == 2
    %             obj.EdgeSegmentationInfo.Edge2.Points = [x_new_points' y_new_points' z_new_points'];
    %         else
    %             obj.EdgeSegmentationInfo.Edge3.Points = [x_new_points' y_new_points' z_new_points'];
    %         end
    %         obj.EdgeSegmentationInfo.done(edg) = 1;
    new_points = [new_points;current_generated_points];
    %     else %% condition to check if its already divided
    %         Points_for_triangulation = [Points_for_triangulation;obj.EdgeSegmentationInfo.Edge1.Points];
    %     end
end
end