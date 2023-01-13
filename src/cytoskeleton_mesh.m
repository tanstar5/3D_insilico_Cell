classdef cytoskeleton_mesh
    properties
        %% Basic Properties
        ID
        Points
        triangle
        neighbour_IDs
        neighbour_TRs
        %% Derived properties  
        area_of_mesh
        quality
        angles_of_triangle
        %% Dependent properties based on its neighbors
        EdgeSegmentationInfo
        SurfaceSegmentationInfo
        local_membrane_meshConnectivity
        
        %% New segmentation variables
        edge_list
        edge_segmemntation_noter
        
        edge12_indices
        edge23_indices
        edge31_indices
        
        edge12_indices_coors
        edge23_indices_coors
        edge31_indices_coors
        
        surface_points_indices
        surface_points_coors
        
        
        %% Biological and physical properties
        membrane_mesh
        boundary
        
    end
    
    methods
        %% Initializing functions
        function obj = cytoskeleton_mesh(ID)
            if nargin ~= 0
                obj.ID = ID;
                segmented_noter.done = [0;0;0];
                obj.EdgeSegmentationInfo = segmented_noter;
            end
        end
        
        function obj = assign_points(obj,Points)
            obj.Points = Points;
        end
        
        function obj = assign_triangle(obj,tr)
            obj.triangle = tr;
        end
        
        function obj = assign_neighbours(obj,neighbour_IDs)
            obj.neighbour_IDs = neighbour_IDs;
        end
        
        function obj = assign_neighbour_triangles(obj,neighbour_TRs)
            obj.neighbour_TRs = neighbour_TRs;
        end
        
        function obj = assign_edges(obj,tr)
            obj.edge_list = [tr(1) tr(2);
                tr(2) tr(3);
                tr(1) tr(3)];
            obj.edge_list = sort(obj.edge_list,2);
            obj.edge_segmemntation_noter = [0,0,0]';
        end
        
        
        %% Derived property functions
        function [common_TR] = find_opposite_triangles_to_edge(obj,edge)
            %             triangle_of_object = obj.triangle;
            neighbour_of_object = obj.neighbour_IDs;
            neighbourTriangles_of_object = obj.neighbour_TRs;
            logical_mask = ismember(neighbourTriangles_of_object,edge);
            triangle_indicator = sum(double(logical_mask),2);
            logical_triangle_ID = triangle_indicator==2;
            common_TR(1) = neighbour_of_object(logical_triangle_ID);
            common_TR(2) = obj.ID;
        end
        
        function [triangle_info,obj] = triangle_quality_for_segmentation(obj)
            %% This function returns quality value 1 if the triangle is acute and 0 if obtuse
            triangle_coors = obj.Points;
            %% checking size criterion
            edge_vecs = ...
                [ triangle_coors(1,:)-triangle_coors(2,:);...
                triangle_coors(2,:)-triangle_coors(3,:);...
                triangle_coors(3,:)-triangle_coors(1,:)];
            edge_lengths = vecnorm(edge_vecs,2,2);
            semi_perimeter = sum(edge_lengths)/2;
            area = (semi_perimeter*( semi_perimeter-edge_lengths(1) )*( semi_perimeter-edge_lengths(2) )*...
                ( semi_perimeter-edge_lengths(3) ))^(.5);
            %% Calculating angles
            angle_e1_e2 = acos(dot( edge_vecs(1,:),-edge_vecs(2,:) )/edge_lengths(1)/edge_lengths(2))/pi*180;
            angle_e2_e3 = acos(dot( edge_vecs(2,:),-edge_vecs(3,:) )/edge_lengths(2)/edge_lengths(3))/pi*180;
            angle_e3_e1 = acos(dot( edge_vecs(3,:),-edge_vecs(1,:) )/edge_lengths(3)/edge_lengths(1))/pi*180;
            %% Calculating the characteristics of the plane
            [PlaneCharacteristics] = equation_of_plane(triangle_coors(1,:),triangle_coors(2,:),triangle_coors(3,:));
            %% Writing values
            if max([angle_e1_e2,angle_e2_e3,angle_e3_e1])>90
                triangle_info.quality = 0;
            else
                triangle_info.quality = 1;
            end
            triangle_info.edge_vecs = edge_vecs;
            triangle_info.edge_lengths = edge_lengths;
            angles.angles_btw_e1_e2 = angle_e1_e2;
            angles.angles_btw_e2_e3 = angle_e2_e3;
            angles.angles_btw_e3_e1 = angle_e3_e1;
            triangle_info.angles = angles;
            triangle_info.area = area;
            triangle_info.PlaneCharacteristics = PlaneCharacteristics;
            obj.quality = triangle_info.quality;
            obj.angles_of_triangle = triangle_info.angles;
            obj.area_of_mesh = area;
        end
        
        
        %% Dependent property functions
        function [obj_list,new_points_list] = segment_all_edge(obj_list,unit_segment_length)
            new_points_list = [];
            total_points_on_edges = 0;
            for obj_id = 1:size(obj_list,1)
                obj = obj_list(obj_id);
                triangle_current = obj.triangle;
                triangle_points = obj.Points;
                edges = [triangle_current(1),triangle_current(2);...
                    triangle_current(2),triangle_current(3);...
                    triangle_current(3),triangle_current(1)];
                Points_for_triangulation = triangle_points;
                for edg = 1:3
                    edge_current = edges(edg,:);
                    edge_current = sort(edge_current,2);
                    edge_current_vect = triangle_points(triangle_current==edge_current(2),:) - triangle_points(triangle_current==edge_current(1),:);
                    edge_points = [triangle_points(triangle_current==edge_current(2),:);triangle_points(triangle_current==edge_current(1),:)];
                    unit_edge_current_vect = edge_current_vect./vecnorm(edge_current_vect,2,2);
                    segmented_noter_old = obj.EdgeSegmentationInfo.done(edg);
                    
                    if segmented_noter_old == 0
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
                        %                         total_points_on_edges = total_points_on_edges
                        new_index = total_points_on_edges+(1:length(x_new_points))';
                        total_points_on_edges = total_points_on_edges + length(new_index);
                        %                         edge_indexing = [edge_indexing;new_index];
                        if edg == 1
                            obj.EdgeSegmentationInfo.Edge1.Points = [x_new_points' y_new_points' z_new_points'];
                            obj.EdgeSegmentationInfo.Edge1.Points_idx = new_index;
                        elseif edg == 2
                            obj.EdgeSegmentationInfo.Edge2.Points = [x_new_points' y_new_points' z_new_points'];
                            obj.EdgeSegmentationInfo.Edge2.Points_idx = new_index;
                        else
                            obj.EdgeSegmentationInfo.Edge3.Points = [x_new_points' y_new_points' z_new_points'];
                            obj.EdgeSegmentationInfo.Edge3.Points_idx = new_index;
                        end
                        
                        obj.EdgeSegmentationInfo.done(edg) = 1;
                        %% Setting the values for neighbouring triangles the segemented edge information
                        [common_TR] = find_opposite_triangles_to_edge(obj,edge_current);
                        ID_neighbour = common_TR(1);
                        obj_neighbour = obj_list(ID_neighbour);
                        triangle_current_neighbour = obj_neighbour.triangle;
                        %                         triangle_points = obj.Points;
                        edges_neighbour = [triangle_current_neighbour(1),triangle_current_neighbour(2);...
                            triangle_current_neighbour(2),triangle_current_neighbour(3);...
                            triangle_current_neighbour(3),triangle_current_neighbour(1)];
                        edges_neighbour_sorted = sort(edges_neighbour,2);
                        index_edge = (1:3)';
                        index_mask = ismember(edges_neighbour_sorted,edge_current,'rows');
                        edge_neighbor_index = index_edge(index_mask);
                        if edge_neighbor_index == 1
                            obj_neighbour.EdgeSegmentationInfo.Edge1.Points = [x_new_points' y_new_points' z_new_points'];
                            obj_neighbour.EdgeSegmentationInfo.Edge1.Points_idx = new_index;
                        elseif edge_neighbor_index == 2
                            obj_neighbour.EdgeSegmentationInfo.Edge2.Points = [x_new_points' y_new_points' z_new_points'];
                            obj_neighbour.EdgeSegmentationInfo.Edge2.Points_idx = new_index;
                        else
                            obj_neighbour.EdgeSegmentationInfo.Edge3.Points = [x_new_points' y_new_points' z_new_points'];
                            obj_neighbour.EdgeSegmentationInfo.Edge3.Points_idx = new_index;
                        end
                        obj_neighbour.EdgeSegmentationInfo.done(edge_neighbor_index) = 1;
                        obj_list(ID_neighbour) = obj_neighbour;
                        new_points_list = [new_points_list;[x_new_points' y_new_points' z_new_points']];
                    else %% condition to check if its already divided
                        Points_for_triangulation = [Points_for_triangulation;obj.EdgeSegmentationInfo.Edge1.Points];
                    end
                end
                ID_current = obj.ID;
                obj_list(ID_current) = obj;
            end
            %[common_TR] = find_opposite_triangles_to_edge(obj,edge);
        end
        
        function [obj_list,new_points_list] = insert_surface_node_points(obj_list,unit_segment_length)
            new_points_list = [];
            total_points_on_edges = 0;
            for obj_id = 1:size(obj_list,1)
                obj = obj_list(obj_id);
                %% Calculate Angles, Area,
                %% Put points on scaled down triangles on its edges based on unit_segment_length
                %% Calculating angles
                %% Determine if obtuse or acute
                %% Put points on the median line for obtuse triangles based on unit_segment_length.
                [~,obj] = triangle_quality_for_segmentation(obj);
                area = obj.area_of_mesh;
                area_thresh = 0.5*unit_segment_length^2*(3)^(.5)/2;
                large_area = area>area_thresh;
                acute_or_obtuse = obj.quality; %%if acute then = 1
                %% Possibility1: Area>AreaThresh and is acute >> put nodes on scaled triangles
                if large_area==1
                    %                     area_new = area;
                    Points_triangle = obj.Points;
                    %                     scales_triangular_points = Points_triangle;
                    new_triangles = [];
                    on_edge_points = [];
                    iteration = 1;
                    [~,incenter_rad] = scale_triangle_based_on_barycentric(Points_triangle,iteration*unit_segment_length);
                    while((iteration)*unit_segment_length < .95*incenter_rad)
                        [scales_triangular_points,~] = scale_triangle_based_on_barycentric(Points_triangle,iteration*unit_segment_length);
                        iteration = iteration + 1;
                        %                         [area_new] = area_triangle(scales_triangular_points);
                        new_triangles = [new_triangles;scales_triangular_points];
                        on_edge_points = [on_edge_points;segment_edge(scales_triangular_points,unit_segment_length)];
                    end
                    %% Possibility2: Area>AreaThresh and is obtuse >> put nodes on median line
                    
                    %% Possibility3: Area<AreaThresh and is obtuse >> Check for super obtuse if so then insert a single incenter or centrod point
                end
                ID_current = obj.ID;               
%                 [Points_new] = cluster_points(unique([new_triangles;on_edge_points],'rows'),unit_segment_length);
%                 [Points_new] = cluster_points(unique(Points_new,'rows'),unit_segment_length);
%                 [Points_new] = cluster_points([on_edge_points],unit_segment_length);
%                 obj.SurfaceSegmentationInfo.Points = Points_new;
                obj.SurfaceSegmentationInfo.Points = unique([new_triangles;on_edge_points],'rows');               
                new_points_list = [new_points_list; obj.SurfaceSegmentationInfo.Points];
                new_points_created = size(obj.SurfaceSegmentationInfo.Points,1);
                new_index = total_points_on_edges+(1:new_points_created)';
                obj.SurfaceSegmentationInfo.Points_idx = new_index;
                total_points_on_edges = total_points_on_edges + length(new_index);
                obj_list(ID_current) = obj;
                fprintf('%d CYTO TRIANGLE meshed\n',ID_current);
            end
        end
        
        function [obj_list,new_points_list] = make_membrane_patch(obj_list,unit_segment_length)
            new_points_list = [];
            total_points_on_edges = 0;
            %% Take advantage of already initially numbered points to index the new points and rename the triangulations accordingly
            for obj_id = 1:size(obj_list,1)
                obj = obj_list(obj_id);
                
                points_for_triangulation = [obj.Points;...
                    obj.EdgeSegmentationInfo.Edge1.Points;...
                    obj.EdgeSegmentationInfo.Edge2.Points;...
                    obj.EdgeSegmentationInfo.Edge3.Points;
                    obj.SurfaceSegmentationInfo.Points];
                if isempty(obj.SurfaceSegmentationInfo.Points)==1
                    size(obj.SurfaceSegmentationInfo.Points,1)
                end
                internal_nodal_point_start_ind = size(points_for_triangulation,1) - size(obj.SurfaceSegmentationInfo.Points,1)+1;
                
                %% Transforming points to 2D plane
                [output] = equation_of_plane(obj.Points(1,:),obj.Points(2,:),obj.Points(3,:));
                origin_coor = (obj.Points(1,:)+obj.Points(2,:)+obj.Points(3,:))/3;
                [transformed_coor_in_plane] = transform_coor_to_plane(points_for_triangulation,output.BasisVectors,origin_coor);
                transformed_coor_in_plane2D = transformed_coor_in_plane(:,1:2);
                DT = delaunay(transformed_coor_in_plane2D);
                
                
                
                %% Deleting triangles with zero area due to error in delaunay because of points on same line
                triangles_global_mesh = DT;
                all_points = transformed_coor_in_plane2D;
                all_points(:,3) = zeros(size(all_points,1),1);
                area_zero_noter_logical = zeros(size(triangles_global_mesh,1),1);
                for tri = 1:size(triangles_global_mesh,1)
                    tri_ind = triangles_global_mesh(tri,:);
                    Points_current = all_points(tri_ind,:);
                    
                    
                    
                    if abs(real(area_triangle(Points_current)))<=min(10e-6,unit_segment_length*unit_segment_length*.5/10000)
                        area_zero_noter_logical(tri) = 1;
                        fprintf('Zero area triangle detected with area %d in triangle %d\n',area_triangle(Points_current), tri );
                    else
                        area_zero_noter_logical(tri) = 0;
                    end
                end
                triangles_global_mesh(area_zero_noter_logical==1,:)=[];
                

                new_points_list = [new_points_list; obj.SurfaceSegmentationInfo.Points];
                if (isempty(obj.SurfaceSegmentationInfo.Points)==1)
                    new_points_created = 0;
                else
                    new_points_created = size(obj.SurfaceSegmentationInfo.Points,1); 
                end
                
                new_index = total_points_on_edges+(1:new_points_created)';
                obj.SurfaceSegmentationInfo.Points_idx = new_index;
                total_points_on_edges = total_points_on_edges + length(new_index);
                
                %%
%                 DT = triangulation(triangles_global_mesh,all_points);
                membrane_mesh_triangulation = triangulation(triangles_global_mesh,points_for_triangulation);
                [membrane_mesh_triangulation] = allign_triangulation_CounterClock(membrane_mesh_triangulation,[0,0,0]);
                obj.membrane_mesh = membrane_mesh_triangulation;
                %%
                
                ID_current = obj.ID;
                obj_list(ID_current) = obj;
                fprintf('%d Mesh Triangulated\n',ID_current);
            end
        end
        
        function [obj_list] = refine_fine_membrane_mesh(obj_list)
            for obj_id = 1:size(obj_list,1)
                obj = obj_list(obj_id);
                if obj.quality == 0
                    membrane_mesh_current = obj.membrane_mesh;
                    ConnectivityListCurrent = membrane_mesh_current.ConnectivityList;
                    ConnectivityListCurrent_sorted = sort(ConnectivityListCurrent,2);
                    edges_from_triangulation = [ ConnectivityListCurrent_sorted(:,1),ConnectivityListCurrent_sorted(:,2);...
                        ConnectivityListCurrent_sorted(:,2),ConnectivityListCurrent_sorted(:,3);...
                        ConnectivityListCurrent_sorted(:,3),ConnectivityListCurrent_sorted(:,1);];
                    
                    
                    points_current = membrane_mesh_current.Points;
                    %% Making the constraint mask for points on edges which are not allowed to move for preserving mesh boundary condition
                    barbed_mesh_points = obj.Points;
                    num_barbed_mesh_points = size(barbed_mesh_points,1);
                    
                    barbed_mesh_edge_points = [obj.EdgeSegmentationInfo.Edge1.Points ;...
                        obj.EdgeSegmentationInfo.Edge2.Points;...
                        obj.EdgeSegmentationInfo.Edge3.Points];
                    num_barbed_mesh_edge_points = size(barbed_mesh_edge_points,1);
                    
                    surface_points = obj.SurfaceSegmentationInfo.Points;
                    num_surface_points = size(surface_points,1);
                    
                    constraint_mask = zeros(num_barbed_mesh_points+num_barbed_mesh_edge_points+num_surface_points,1);
                    constraint_mask(1:num_barbed_mesh_points+num_barbed_mesh_edge_points) = 1;
                    point_indices = (1:num_barbed_mesh_points+num_barbed_mesh_edge_points+num_surface_points)';
                    indices_to_modify = point_indices(constraint_mask==0);
                    
                    [output] = equation_of_plane(obj.Points(1,:),obj.Points(2,:),obj.Points(3,:));
                    origin_coor = (obj.Points(1,:)+obj.Points(2,:)+obj.Points(3,:))/3;
                    [transformed_coor_in_plane] = transform_coor_to_plane(points_current,output.BasisVectors,origin_coor);
                    transformed_coor_in_plane_new = transformed_coor_in_plane;
                    for to_mod_points = 1:length(indices_to_modify)
                        ind_to_mod = indices_to_modify(to_mod_points);
                        point_to_update = transformed_coor_in_plane(ind_to_mod,:);
                        row_mask_col1 = ismember(edges_from_triangulation(:,1),ind_to_mod);
                        row_mask_col2 = ismember(edges_from_triangulation(:,2),ind_to_mod);
                        row_mask = row_mask_col1|row_mask_col2;
                        
                        imp_rows = edges_from_triangulation(row_mask,:);
                        neigbours = imp_rows(  not(ismember(imp_rows,ind_to_mod))   );
                        
                        point_neighbours_inplane = transformed_coor_in_plane(neigbours,:);
                        transformed_coor_in_plane_new(ind_to_mod,:) = mean(point_neighbours_inplane,1);
                        
                        
                        
                    end
                    [points_modified] = inverse_transform_to_XYZ_coor(transformed_coor_in_plane_new,output.BasisVectors,origin_coor);
                    obj.membrane_mesh = triangulation(ConnectivityListCurrent,points_modified);
                    obj.SurfaceSegmentationInfo.Points = points_modified(indices_to_modify,:);
                    ID_current = obj.ID;
                    obj_list(ID_current) = obj;
                    fprintf('%d Mesh Refined\n',ID_current);
                end
            end
        end
        
        function [global_membrane_mesh,obj_list] =  ...
                concatenate_membrane_patches(obj_list,initial_mesh_points,new_edgePoints_inserted_points,new_surfacePoints_inserted_points)
            ConnectivityList_all = [];
            number_of_initial_mesh_points = size(initial_mesh_points,1);
            number_of_points_onEdges = size(new_edgePoints_inserted_points,1);
            for ob = 1:size(obj_list,1)
                obj = obj_list(ob);
                ConnectivityList = obj.membrane_mesh.ConnectivityList;
                point_ind_map = [obj.triangle';...
                    number_of_initial_mesh_points + obj.EdgeSegmentationInfo.Edge1.Points_idx;...
                    number_of_initial_mesh_points + obj.EdgeSegmentationInfo.Edge2.Points_idx;...
                    number_of_initial_mesh_points + obj.EdgeSegmentationInfo.Edge3.Points_idx;...
                    number_of_initial_mesh_points + number_of_points_onEdges + obj.SurfaceSegmentationInfo.Points_idx];
                point_ind_map = point_ind_map';
                
                if isempty(obj.SurfaceSegmentationInfo.Points_idx) ==1 && isempty(obj.EdgeSegmentationInfo.Edge1.Points_idx)==1 ...
                    && isempty(obj.EdgeSegmentationInfo.Edge2.Points_idx)==1 ...
                    && isempty(obj.EdgeSegmentationInfo.Edge3.Points_idx)==1
                    disp('Error');
                end
%                 ConnectivityList_mapped = NaN(size(ConnectivityList));
                ConnectivityList_mapped = point_ind_map(ConnectivityList);
                obj.local_membrane_meshConnectivity = ConnectivityList_mapped;
                ID_current = obj.ID;
                obj_list(ID_current) = obj;
                ConnectivityList_all = [ConnectivityList_all;ConnectivityList_mapped];
            end
            Points_all = [initial_mesh_points;new_edgePoints_inserted_points;new_surfacePoints_inserted_points];
            global_membrane_mesh = triangulation(ConnectivityList_all,Points_all);
        end
        
        %% New modified efficient functions for mesh generation from cytoskeleton objects
        function [internal_nodal_points,new_points_mask] = insert_face_nodes(obj,unit_segment_length)
            %% Check acuteness
            triangle_current = obj.triangle;
            points_triangle = obj.Points;
            [angles_current] = calculate_angles_in_triangle(points_triangle);
            angle_min = min(angles_current);
            angle_max = max(angles_current);
            
            acute_indicator = angle_max*180/pi<= 90;
            
            %% check area
            edge_12_length = vecnorm(points_triangle(1,:)-points_triangle(2,:));
            edge_23_length = vecnorm(points_triangle(2,:)-points_triangle(3,:));
            edge_31_length = vecnorm(points_triangle(3,:)-points_triangle(1,:));
            
            area_current = area_triangle(points_triangle);
            area_big_indicator = area_current>=3^(0.5)/4* (edge_12_length + edge_23_length + edge_31_length)/3;
            
            %% find nodes
            on_edge_points = [];
            [~,incenter_rad] = scale_triangle_based_on_barycentric(points_triangle,1*unit_segment_length);
            if incenter_rad>.7*unit_segment_length && area_big_indicator && acute_indicator
                num_segments = ceil(incenter_rad/unit_segment_length);
                segment_length = incenter_rad/num_segments;
                inward_lengths = segment_length:segment_length:incenter_rad-segment_length;
                new_triangles = [];
                
                for ind_lgth_id = 1:length(inward_lengths)                    
                    [scales_triangular_points,~] = scale_triangle_based_on_barycentric(points_triangle,inward_lengths(ind_lgth_id));                    
                    new_triangles = [new_triangles;scales_triangular_points];
                    on_edge_points = [on_edge_points;segment_edge(scales_triangular_points,unit_segment_length)];
                end
            end          
        internal_nodal_points = on_edge_points;  
        end
            

        function [membrane_mesh,obj_list] = create_global_membrane_mesh(obj_list,initial_mesh_points,unit_membrane_mesh_length)
            %% Create unique points for edge segment points with proper address
            points_list = initial_mesh_points;
            triangle_global = [];
            for obj_id = 1:size(obj_list,1)       
                
                obj = obj_list(obj_id);
                triangle_current = obj.triangle;
                points_current = obj.Points;
                edges_current_triangle = obj.edge_list;
                point_indices_current = [triangle_current(1,:)];
                triangulation_points = [points_current];
                for edge_id = 1:3
                   num_points_global = size(points_list,1);
                   edge_current = sort(edges_current_triangle(edge_id,:));
                   [common_TR] = find_opposite_triangles_to_edge(obj,edge_current);
                   
                   if obj.edge_segmemntation_noter(edge_id) == 0
                       ID_neighbour = common_TR(1);
                       % insert new points in the edge
                       point_start =  points_list(edge_current(1),:);
                       point_end =  points_list(edge_current(2),:);
                       
                       mag_distance = vecnorm( (point_end - point_start) );
                       direct_vect_unit = (point_end - point_start)/mag_distance;
                       
                       if mag_distance >= 1.5*unit_membrane_mesh_length
                            num_divs = ceil(mag_distance/unit_membrane_mesh_length);
                            dr_mag = mag_distance/num_divs;
                            new_edge_points_mag = dr_mag:dr_mag:mag_distance-dr_mag;
                            new_edge_points_x = point_start(1) + new_edge_points_mag*direct_vect_unit(1);
                            new_edge_points_y = point_start(2) + new_edge_points_mag*direct_vect_unit(2);
                            new_edge_points_z = point_start(3) + new_edge_points_mag*direct_vect_unit(3);
                            
                            new_edge_points = [new_edge_points_x',new_edge_points_y',new_edge_points_z'];
                            num_new_points = size(new_edge_points,1);
                            new_point_indices = num_points_global + (1:num_new_points);
                            if edge_id == 1
                                obj.edge12_indices = new_point_indices;
                                obj.edge12_indices_coors = new_edge_points;
                            elseif edge_id == 2
                                obj.edge23_indices = new_point_indices;
                                obj.edge23_indices_coors = new_edge_points;
                            else
                                obj.edge31_indices = new_point_indices;
                                obj.edge31_indices_coors = new_edge_points;
                            end
                            obj.edge_segmemntation_noter(edge_id) = 1;                     
                            
                            
                            % Find the corresponding neighbour edge in the
                            % neighbour triangle
                            obj_neighbour = obj_list(ID_neighbour);
                            edges_current_neighbour_triangle = obj_neighbour.edge_list;
                            edge_no_mask = ismember( sort(edges_current_neighbour_triangle,2),edge_current,'rows' );
                            obj_neighbour.edge_segmemntation_noter(edge_no_mask) = 1;
                            edge_index = 1:3;
                            edge_id_neighbour = edge_index(edge_no_mask);
                            if edge_id_neighbour == 1
                                obj_neighbour.edge12_indices = new_point_indices;
                                obj_neighbour.edge12_indices_coors = new_edge_points;
                            elseif edge_id_neighbour == 2
                                obj_neighbour.edge23_indices = new_point_indices;
                                obj_neighbour.edge23_indices_coors = new_edge_points;
                            else
                                obj_neighbour.edge31_indices = new_point_indices;
                                obj_neighbour.edge31_indices_coors = new_edge_points;
                            end                            
                            obj_list(ID_neighbour) = obj_neighbour;                            
                            points_list = [points_list;new_edge_points];
                            point_indices_current = [point_indices_current,new_point_indices];
                            triangulation_points = [triangulation_points;new_edge_points];
                       else
                           obj.edge_segmemntation_noter(edge_id) = 1;
                       end
                       
                   else
                       if edge_id == 1
                            existing_point_indices = obj.edge12_indices;
                            existing_edge_points = obj.edge12_indices_coors;
                       elseif edge_id == 2
                            existing_point_indices = obj.edge23_indices;
                            existing_edge_points = obj.edge23_indices_coors;
                       else
                            existing_point_indices = obj.edge31_indices;
                            existing_edge_points = obj.edge31_indices_coors;
                       end
                       point_indices_current = [point_indices_current,existing_point_indices];
                       triangulation_points = [triangulation_points;existing_edge_points];
                       
                   end             
                end
                % Surface segmemntation
                [internal_nodal_points] = insert_face_nodes(obj,unit_membrane_mesh_length);
                num_new_face_points = size(internal_nodal_points,1);
                num_points_global = size(points_list,1);
                points_list = [points_list;internal_nodal_points];
%                 num_points_global = size(points_list,1);
                new_point_indices = num_points_global + (1:num_new_face_points);
                triangulation_points = [triangulation_points;internal_nodal_points];
                point_indices_current = [point_indices_current,new_point_indices];
                
                % Triangulating
                %% Transforming points to 2D plane
                [output] = equation_of_plane(obj.Points(1,:),obj.Points(2,:),obj.Points(3,:));
                origin_coor = (obj.Points(1,:)+obj.Points(2,:)+obj.Points(3,:))/3;
                [transformed_coor_in_plane] = transform_coor_to_plane(triangulation_points,output.BasisVectors,origin_coor);
                transformed_coor_in_plane2D = transformed_coor_in_plane(:,1:2);
                DT = delaunay(transformed_coor_in_plane2D);
                DT_global = point_indices_current(DT);
                if max(DT_global(:))>size(points_list,1)
                   disp('ERROR') 
                end
                triangle_global = [triangle_global;DT_global];
                %%%
                
                obj_list(obj_id) = obj;
            end
            membrane_mesh = triangulation(triangle_global,points_list);
          % delete thin triangles
          membrane_mesh_zero_triangle_noter = NaN(size(triangle_global,1),1);
          parfor tri_id = 1:size(triangle_global,1)
              points_tri_current = points_list(triangle_global(tri_id,:),:);
              if area_triangle( points_tri_current)<=10e-6
                  membrane_mesh_zero_triangle_noter(tri_id) = 1;
              else
                  membrane_mesh_zero_triangle_noter(tri_id) = 0;
              end              
          end
          triangle_global( membrane_mesh_zero_triangle_noter==1,:) = [];
           membrane_mesh = triangulation(triangle_global,points_list);
            
            
        end
        
        %% Need to refine global membrane mesh
    end
    
end