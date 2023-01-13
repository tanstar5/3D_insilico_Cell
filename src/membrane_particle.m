classdef membrane_particle
    
    properties
        % Basic identifiers
        ID
        Pos
        type
        spontaneous_curvature
        bending_rigidity
        gaussian_rigidity
        
        % flags
        need_to_update
        constraint
        
        % about immediate neighbour/topology
        neighbours
        neighbours_coor;
        neighbours_normal_old;
        faces
        faceNormals
        faceArea
        Av_vertex
        Av_projected_vertex
        Normal_vect
        Householder_mat
        Projection_operator
        
        % Derived Topological characteristics
        facePairs
        dihedralAngle
        Re
        He
        Ne
        Se_operators
        We
        Sv_operator
        DarbouxFrame
        
        % Energy calculations
        free_energy
        curvature_energy
        local_surface_tension_energy
        
        
        energy_stored
        work_volume
        work_area
        hard_ball_potential
    end
    
    methods (Access = public)
        %% Initializing functions
        function obj = membrane_particle(ID)
            if nargin ~= 0
                obj.ID = ID;
            end
        end
        
        function [obj_list] = membrane_particles_list(first_obj,num_of_particles)
            obj_list(1,1) =  first_obj;
            for particle = 2:num_of_particles
                obj = membrane_particle(particle);
                obj_list(particle,1) =  obj;
            end
        end
        
        function [obj_list,list_id_logical] = assign_type_poles(obj_list,pole_rad_min,pole_rad_max,type_num)
            list_id_logical = zeros(size(obj_list,1),1);
            parfor obj_id = 1:size(obj_list,1)
                obj = obj_list(obj_id);
                pos = obj.Pos;
                if( (pos(1)^2 + pos(2)^2)^.5<=pole_rad_max && (pos(1)^2 + pos(2)^2)^.5>=pole_rad_min )
                    obj.type = type_num;
                    list_id_logical(obj_id) = 1;
                    obj_list(obj_id) = obj;
                end
            end
        end
        
        function [obj_list] = assign_type_randomly_conserved(obj_list,type_list)
            num_objs = size(obj_list,1);
            num_types = length(type_list);
            for obj_ind = 1:num_objs
                randomize_type_ind = randperm(num_types);
                obj_list(obj_ind).type = type_list(randomize_type_ind(1));
            end
            
        end
        
        function [obj_list] = assign_type_properties(obj_list,type_IDs,spontaneous_curvatures,bending_rigidity,gaussian_rigidity)
            all_type_extract = @(ind) obj_list(ind).type;
            all_type = arrayfun(all_type_extract,[1:size(obj_list,1)]);
            interesting_objs_logical = ismember(all_type,type_IDs);
            obj_idx = 1:size(obj_list,1);
            interesting_objs_idx =obj_idx(interesting_objs_logical==1);
            parfor obj_interesting_id = interesting_objs_idx
                obj_current = obj_list(obj_interesting_id);
                obj_current.spontaneous_curvature = spontaneous_curvatures(ismember(type_IDs,obj_current.type));
                obj_current.bending_rigidity = bending_rigidity(ismember(type_IDs,obj_current.type));
                obj_current.gaussian_rigidity = gaussian_rigidity(ismember(type_IDs,obj_current.type));
                obj_list(obj_interesting_id)=obj_current;
            end
            fprintf('   >>>PROPERTIES assigned to total %d types \n',length(type_IDs));
        end
        
        %% Identifying functions
        function [obj,ind_in_list_arranged,corresponding_IDs_arranged] = findobjID(obj_list,ID)
            all_ID_extract = @(ind) obj_list(ind).ID;
            all_ID = arrayfun(all_ID_extract,[1:size(obj_list,1)]);
            all_index = (1:size(obj_list,1));
            if length(ID)<=1
                if length(all_ID(all_ID==ID))>1
                    fprintf('  >>>MULTIPLE objects found \n');
                    obj = NaN;
                    ind_in_list = NaN;
                elseif length(all_ID(all_ID==ID))<1
                    fprintf('  >>>NO objects found \n');
                    obj = NaN;
                    ind_in_list = NaN;
                else
                    obj = obj_list(all_ID==ID);
                    ind_in_list = all_index(all_ID==ID);
                end
            else
                [obj_list_tmp] = membrane_particles_list(membrane_particle(1),length(ID));
                %                 obj_list_tmp = obj_list(ismember(all_ID,ID));
                %                 ind_in_list = all_index(ismember(all_ID,ID));
                %                 obj = obj_list_tmp;
                all_IDs = extract_IDs(obj_list);
                all_pos_in_obj_list = 1:size(obj_list,1);
                all_quered_pos_ID_mask = ismember(all_IDs,ID);
                all_quered_pos_ID_in_list = all_pos_in_obj_list(all_quered_pos_ID_mask);
                
                extracted_obj_ID = all_IDs(all_quered_pos_ID_mask);
                extracted_objs =  obj_list(all_quered_pos_ID_in_list);
                
                obj = extracted_objs;
                ind_in_list = all_quered_pos_ID_in_list;
                corresponding_IDs = all_IDs(all_quered_pos_ID_in_list);
                
                [~,locb] = ismember(ID,corresponding_IDs);
                corresponding_IDs_arranged = corresponding_IDs(locb);
                ind_in_list_arranged = ind_in_list(locb);
                obj = obj(locb);
                
                
                
                
            end
            
        end
        
        %% Topological functions
        function [obj_list] = assign_local_info_to_objects_fromMesh(obj_list,membrane_mesh)
            
            triangles = membrane_mesh.ConnectivityList;
            face_normal_vects_all = faceNormal(membrane_mesh);
            all_points = membrane_mesh.Points;
            normal_vects_at_each_points = NaN(size(all_points)) ;
            % Correcting the orientation of triangles on anticlockwise order seen from an inside point in a cell
            triangles_arranged = NaN(size(triangles));
            face_normal_vects_all_arranged = NaN(size(face_normal_vects_all));
            for tri_id = 1:size(triangles,1)
                tri_current = triangles(tri_id,:);
                tri_points = all_points(tri_current,:);
                tri_centroid = mean(tri_points,1);
                if sign( dot(face_normal_vects_all(tri_id,:),tri_centroid-[0,0,0] ))<0
                    triangles_arranged(tri_id,:) = fliplr(tri_current);
                    face_normal_vects_all_arranged(tri_id,:) = -face_normal_vects_all(tri_id,:);
                    disp('   >>>TRIANGULATION ANTICLOCKED\n')
                else
                    triangles_arranged(tri_id,:) = tri_current;
                    face_normal_vects_all_arranged(tri_id,:) = face_normal_vects_all(tri_id,:);
                end
            end
            triangles = triangles_arranged;
            face_normal_vects_all = face_normal_vects_all_arranged;
            parfor particle = 1:size(obj_list,1)
                obj = obj_list(particle);
                obj.Pos = membrane_mesh.Points(particle,:);
                % Finding neigbours
                triangle_mask = ismember(triangles,particle);
                triangle_mask_rows = triangle_mask(:,1)|triangle_mask(:,2)|triangle_mask(:,3);
                interesting_triangles = triangles(triangle_mask_rows,:);
                obj.neighbours = unique(interesting_triangles(not(ismember(interesting_triangles,particle))));
                obj.neighbours_coor = membrane_mesh.Points(obj.neighbours,:);
                fprintf('membrane point %d created \n',particle);
                % Finding Normal vect at current vertex and Darboux frame
                faces_of_current_vertex = triangles(triangle_mask_rows,:);
                face_area = NaN(size(faces_of_current_vertex,1),1);
                for face_obj = 1:size(faces_of_current_vertex,1)
                    face_current = faces_of_current_vertex(face_obj,:);
                    Points_vertices = membrane_mesh.Points(face_current,:);
                    face_area(face_obj) = abs(real(area_triangle(Points_vertices)));
                end
                face_normal_vect_current = face_normal_vects_all(triangle_mask_rows,:);
                edges_faces_all = [faces_of_current_vertex(:,1),faces_of_current_vertex(:,2);...
                    faces_of_current_vertex(:,2),faces_of_current_vertex(:,3);...
                    faces_of_current_vertex(:,3),faces_of_current_vertex(:,1)];
                commonVertex = obj.ID;
                common_vertex_edge_logical = ismember(edges_faces_all,commonVertex);
                common_vertex_edge_logical = common_vertex_edge_logical(:,1)|common_vertex_edge_logical(:,2);
                edges_on_circumference = edges_faces_all(common_vertex_edge_logical==0,:);
                edge_ordered_indicator = zeros(size(edges_on_circumference,1),1);
                loop_no = 1;
                id_face = [];
                while not(isempty( edge_ordered_indicator(edge_ordered_indicator==0)))==1
                    if loop_no == 1
                        edge_current = edges_on_circumference(loop_no,:);
                    else
                        if particle == 3526
                            disp('')
                        end
                        vertex_of_interest = edge_current(2);
                        edge_list_to_scan_logical = not(ismember(edges_on_circumference,edge_current,'rows'));
                        edge_list_to_scan = edges_on_circumference(edge_list_to_scan_logical,:);
                        edge_current_mask = ismember(edge_list_to_scan,vertex_of_interest);
                        edge_current_mask = edge_current_mask(:,1)|edge_current_mask(:,2);
                        edge_current = edge_list_to_scan(edge_current_mask,:);
                    end
                    corresponding_face_find_mask = ismember(faces_of_current_vertex,edge_current);
                    corresponding_face_find_mask = double(corresponding_face_find_mask);
                    corresponding_face_find_mask = sum(corresponding_face_find_mask,2);
                    id_face = [id_face find(corresponding_face_find_mask==2)];
                    edge_ordered_indicator(corresponding_face_find_mask==2) = 1;
                    loop_no = loop_no + 1;
                end
                % Checking the right direction test just by looking into the order of the ordered first two face.
                obj.faces = faces_of_current_vertex(id_face,:);
                obj.faceNormals = face_normal_vect_current(id_face,:);
                obj.faceArea = face_area(id_face,:);
                % if cross(face2_normal,face1normal).(common_vertex_vect_fromcenter_to_vertex) >0
                % Then in correct order. Else reverese it
                normal_cross_vect = cross(obj.faceNormals(2,:),obj.faceNormals(1,:));
                common_edge = intersect(obj.faces(1,:),obj.faces(2,:));
                common_neigbourhood_vertex = common_edge( not(ismember(common_edge,commonVertex)) );
                edge_outward_vect = membrane_mesh.Points(common_neigbourhood_vertex,:) - membrane_mesh.Points(commonVertex,:);
                decision = dot(edge_outward_vect, normal_cross_vect);
                if decision<=0
                    obj.faces = flipud(obj.faces);
                    obj.faceNormals = flipud(obj.faceNormals);
                    obj.faceArea = flipud(obj.faceArea);
                end
                weights_based_area = obj.faceArea./sum( obj.faceArea );
                obj.Normal_vect = [ sum(weights_based_area.*obj.faceNormals(:,1)),...
                    sum(weights_based_area.*obj.faceNormals(:,2)),...
                    sum(weights_based_area.*obj.faceNormals(:,3))];
                obj.Normal_vect = obj.Normal_vect./vecnorm(obj.Normal_vect);
                obj.Projection_operator = diag([1,1,1]) - obj.Normal_vect'*obj.Normal_vect;
                normal_vects_at_each_points(particle,:) =  obj.Normal_vect;
                if vecnorm([0,0,1]-obj.Normal_vect)>vecnorm([0,0,1]+obj.Normal_vect)
                    W_for_householder = [0,0,1]-obj.Normal_vect;
                else
                    W_for_householder = [0,0,1]+obj.Normal_vect;
                end
                W_for_householder = W_for_householder/vecnorm(W_for_householder);
                obj.Householder_mat = diag([1,1,1]) - 2*(W_for_householder'*W_for_householder);
                
                % Calculating dihedral angle between each faces.
                faces_current = obj.faces;
                total_faces = size(faces_current,1);
                face_indices = [1:total_faces,1];
                face_pairs = [face_indices(1:end-1)',face_indices(2:end)'];
                dihedral_angle = NaN(size(face_pairs,1),1);
                Re_outward_edgeVect = NaN(size(face_pairs,1),3);
                He_curvature = NaN(size(face_pairs,1),1);
                Ne_bisector  = NaN(size(face_pairs,1),3);
                Se_operator_mat = NaN(3,3,size(face_pairs,1));
                We_weighting = NaN(size(face_pairs,1),1);
                Sv_update = zeros(3,3);
                for pair = 1:size(face_pairs,1)
                    face_pair_current = face_pairs(pair,:);
                    face1 = obj.faces(face_pair_current(1),:);
                    face2 = obj.faces(face_pair_current(2),:);
                    face1_current_normal = obj.faceNormals(face_pair_current(1),:);
                    face2_current_normal = obj.faceNormals(face_pair_current(2),:);
                    Ne_bisector(pair,:) = (face1_current_normal + face2_current_normal);
                    Ne_bisector(pair,:) = Ne_bisector(pair,:)/vecnorm(Ne_bisector(pair,:));
                    if sign(dot(Ne_bisector(pair,:),obj.Normal_vect) )<0
                        disp('error');
                        particle
                    end
                    common_edge = intersect(face1,face2);
                    common_neighbour_vertex = common_edge(not(ismember(common_edge,commonVertex)));
                    Re_outward_edgeVect(pair,:) = membrane_mesh.Points(common_neighbour_vertex,:) -...
                        membrane_mesh.Points(commonVertex,:);
                    dihedral_angle(pair) = -acos(dot(face1_current_normal,face2_current_normal))+pi;
                    He_curvature(pair) = (real(2*vecnorm(Re_outward_edgeVect(pair,:))*cos(dihedral_angle(pair)/2)));
                    Se_operator_mat(:,:,pair) = He_curvature(pair)*[cross( Re_outward_edgeVect(pair,:)...
                        /vecnorm(Re_outward_edgeVect(pair,:)),Ne_bisector(pair,:) )]'*...
                        [cross( Re_outward_edgeVect(pair,:)/vecnorm(Re_outward_edgeVect(pair,:)),Ne_bisector(pair,:) )];
                    We_weighting(pair,:) = dot(obj.Normal_vect, Ne_bisector(pair,:));
                    Sv_update = Sv_update + We_weighting(pair,:)*obj.Projection_operator'*Se_operator_mat(:,:,pair)*...
                        obj.Projection_operator;
                end
                obj.facePairs = face_pairs;
                obj.dihedralAngle = dihedral_angle;
                obj.Re = Re_outward_edgeVect;
                obj.He = He_curvature;
                obj.Ne = Ne_bisector;
                obj.Se_operators = Se_operator_mat;
                obj.We = We_weighting;
                Av = sum(obj.faceArea)/3;
                obj.Av_vertex = Av;
                obj.Sv_operator = Sv_update;
                [Darboux_frame_eigVects,Darboux_frame_eigVals] = eig(obj.Sv_operator,'vector');
                
                % Arranging the darboux frame
                Darboux_frame_eigVals = Darboux_frame_eigVals';
                idx_ortho_noter = zeros(3,1);
                idx = (1:3);
                Darboux_frame_eigVects_arranged = NaN(3,3);
                Darboux_frame_eigVals_arranged = NaN(1,3);
                if particle == 1002
                    disp('error');
                end
                dot_product_noter = NaN(1,3);
                for eig_vec = 1:3
                    eig_vector = Darboux_frame_eigVects(:,eig_vec);
                    eig_vector_unit = eig_vector/vecnorm(eig_vector);
                    dot_product_noter(eig_vec) = abs(dot(eig_vector_unit,obj.Normal_vect));
                    if abs(dot(eig_vector_unit,obj.Normal_vect))>=0.9
                        idx_ortho_noter(eig_vec) = 1;
                    end
                end
                if isempty( idx_ortho_noter(idx_ortho_noter==1)  )==0
                    Darboux_frame_eigVects_arranged(:,3) = Darboux_frame_eigVects(:,idx_ortho_noter==1);
                    Darboux_frame_eigVals_arranged(:,3) = Darboux_frame_eigVals(:,idx_ortho_noter==1);
                    idx_remaining = idx(idx_ortho_noter==0);
                    [in_plane_eigenVals,id] = sort(Darboux_frame_eigVals(idx_remaining));
                    Darboux_frame_eigVals_arranged(:,1:2) = in_plane_eigenVals;
                    Darboux_frame_eigVects_remaining = Darboux_frame_eigVects(:,idx_remaining);
                    Darboux_frame_eigVects_arranged(:,1:2) = Darboux_frame_eigVects_remaining(:,id);
                    obj.DarbouxFrame{1,1} = Darboux_frame_eigVects_arranged;
                    obj.DarbouxFrame{2,1} = Darboux_frame_eigVals_arranged;
                else
                    fprintf('Anomaly Detected\n');
                    Darboux_frame_eigVects
                    Darboux_frame_eigVals
                    Darboux_frame_eigVects_arranged(:,3) = obj.Normal_vect  ;
                    Darboux_frame_eigVals_arranged(:,3) = 0;
                    %                     idx_remaining = idx(idx_ortho_noter==0);
                    %                     [in_plane_eigenVals,id] = sort(Darboux_frame_eigVals(idx_remaining));
                    Darboux_frame_eigVals_arranged(:,1:2) = [0,0];
                    %                     Darboux_frame_eigVects_remaining = Darboux_frame_eigVects(:,idx_remaining);
                    Darboux_frame_eigVects_arranged(:,1:2) = NaN(3,2);
                    obj.DarbouxFrame{1,1} = Darboux_frame_eigVects_arranged;
                    obj.DarbouxFrame{2,1} = Darboux_frame_eigVals_arranged;
                    
                end
                %                 Darboux_frame_eigVects_arranged(:,3) = Darboux_frame_eigVects(:,idx_ortho_noter==1);
                %                 Darboux_frame_eigVals_arranged(:,3) = Darboux_frame_eigVals(:,idx_ortho_noter==1);
                %                 idx_remaining = idx(idx_ortho_noter==0);
                %                 [in_plane_eigenVals,id] = sort(Darboux_frame_eigVals(idx_remaining));
                
                %                 Darboux_frame_eigVals_arranged(:,1:2) = in_plane_eigenVals;
                %                 Darboux_frame_eigVects_remaining = Darboux_frame_eigVects(:,idx_remaining);
                %                 Darboux_frame_eigVects_arranged(:,1:2) = Darboux_frame_eigVects_remaining(:,id);
                %                 obj.DarbouxFrame{1,1} = Darboux_frame_eigVects_arranged;
                %                 obj.DarbouxFrame{2,1} = Darboux_frame_eigVals_arranged;
                % Energy at vertex
                obj.energy_stored = 0.5*( (mean(obj.DarbouxFrame{2,1}(1:2))-obj.spontaneous_curvature)^2   );
                
                pressure_displacement =  mean(measure_distance_from_plane(obj.Pos,obj.Normal_vect,obj.neighbours_coor));
                projected_face_area = sum(calculate_projected_area_of_faces(obj.Normal_vect,obj.faceNormals,obj.faceArea))/3;
                obj.Av_projected_vertex = projected_face_area;
                volume_of_push = 1/3*projected_face_area*pressure_displacement;
                obj.work_volume = volume_of_push*sign(mean(Darboux_frame_eigVals_arranged));
                obj.work_area = +Av - projected_face_area;
                obj_list(particle) = obj;
                
            end
            % Assigning the normal vects to the neighbours
            parfor particle = 1:size(obj_list,1)
                obj = obj_list(particle);
                obj.neighbours_normal_old = normal_vects_at_each_points(obj.neighbours,:);
                obj_list(particle) = obj;
            end
            
        end
        
        
        %% Derived functions
        function [obj_list,total_energy] = calculate_energy_stored(obj_list,hardball_radius,overlapping_energy)
            energy_counter = NaN(size(obj_list,1));
            
            posx_all_extract = @(id) obj_list(id).Pos(1);
            posy_all_extract = @(id) obj_list(id).Pos(2);
            posz_all_extract = @(id) obj_list(id).Pos(3);
            
            posx_all = (arrayfun(posx_all_extract,1:size(obj_list,1)))';
            posy_all = (arrayfun(posy_all_extract,1:size(obj_list,1)))';
            posz_all = (arrayfun(posz_all_extract,1:size(obj_list,1)))';
            
            for obj_id = 1:size(obj_list,1)
                % Energy due to bending
                obj_current = obj_list(obj_id);
                principle_curvatures = obj_current.DarbouxFrame{2,1}(1:2);
                obj_current.energy_stored = 0.5*obj_current.bending_rigidity*(mean(principle_curvatures)-obj_current.spontaneous_curvature);
                obj_current.energy_stored = 0.5*obj_current.Av_vertex*obj_current.bending_rigidity*(principle_curvatures(1)-obj_current.spontaneous_curvature)^2+...
                    0.5*obj_current.Av_vertex*obj_current.bending_rigidity*(principle_curvatures(2)-obj_current.spontaneous_curvature)^2;
                
                energy_counter(obj_id) = obj_current.energy_stored;
                % Energy due to hard potential
                neigbours_current = obj_current.neighbours;
                neighbours_pos = [ posx_all(neigbours_current) posy_all(neigbours_current) posz_all(neigbours_current)];
                edge_vectors = [obj_current.Pos(1)-neighbours_pos(:,1),...
                    obj_current.Pos(2)-neighbours_pos(:,2),...
                    obj_current.Pos(3)-neighbours_pos(:,3)];
                edge_vect_mag = vecnorm(edge_vectors,2,2);
                if length(edge_vect_mag(edge_vect_mag<=2*hardball_radius))>0
                    obj_current.hard_ball_potential = length(edge_vect_mag(edge_vect_mag<=2*hardball_radius))*overlapping_energy;
                else
                    obj_current.hard_ball_potential = 0;
                end
                
                obj_list(obj_id) = obj_current;
            end
            total_energy = sum(energy_counter,'omitnan');
            fprintf('   >>>ENERGY calculated \n');
        end
        
        
        %% Diffusion functions
        % update change calculation functions -----------------------------
        function [updated_obj] = ...
                calculate_topological_characteristics(obj,hardball_radius,overlapping_energy)
            current_pos = obj.Pos;
            current_ID = obj.ID;
            neighbour_ID_all = obj.neighbours;
            all_neighbour_points = obj.neighbours_coor;
            all_points_insequence = [current_pos;all_neighbour_points];
            
            current_all_id_list = [current_ID,neighbour_ID_all']';
            current_obj_faces = obj.faces;
            current_obj_faces_maskedSequenced = mask_mat_to_sequence(current_obj_faces,current_all_id_list);
            
            [obj.facePairs,~] = update_facePair(obj,obj.faces);
            current_face_pairs = obj.facePairs;
            
            Re_all_edge_vects = NaN(size(current_face_pairs,1),3);
            Ne_all_edge_vects = NaN(size(current_face_pairs,1),3);
            triangulation_current_mesh = triangulation(current_obj_faces_maskedSequenced,all_points_insequence);
            face_normals = faceNormal(triangulation_current_mesh);
            face_Area = NaN(size(current_obj_faces_maskedSequenced,1),1);
            
            for faces_id = 1:size(current_obj_faces_maskedSequenced,1)
                face_current = current_obj_faces_maskedSequenced(faces_id,:);
                current_face_points = all_points_insequence(face_current,:);
                face_Area(faces_id,1) = abs(real(area_triangle(current_face_points)));
                
            end
            weights_based_area = face_Area/sum( face_Area );
            Normal_vect_obj = [ sum(weights_based_area.*face_normals(:,1)),...
                sum(weights_based_area.*face_normals(:,2)),...
                sum(weights_based_area.*face_normals(:,3))];
            Normal_vect_obj  = Normal_vect_obj ./vecnorm(Normal_vect_obj );
            Projection_operator_obj = diag([1,1,1]) - Normal_vect_obj'*Normal_vect_obj;
            
            
            dihedral_angle = NaN(size(current_face_pairs,1),1);
            He_curvature = NaN(size(current_face_pairs,1),1);
            Se_operator_mat = NaN(3,3,size(current_face_pairs,1));
            We_weighting = NaN(size(current_face_pairs,1),1);
            Sv_update = zeros(3,3);
            for edge_id = 1:size(Re_all_edge_vects,1)
                corresponding_face_pair = current_face_pairs(edge_id,:);
                faceid1 = corresponding_face_pair(1);
                faceid2 = corresponding_face_pair(2);
                
                face1 = current_obj_faces_maskedSequenced(faceid1,:);
                face2 = current_obj_faces_maskedSequenced(faceid2,:);
                
                face1_normal = face_normals(faceid1,:);
                face2_normal = face_normals(faceid2,:);
                Ne_all_edge_vects(edge_id,:) = (face1_normal + face2_normal);
                Ne_all_edge_vects(edge_id,:) = Ne_all_edge_vects(edge_id,:)/vecnorm(Ne_all_edge_vects(edge_id,:));
                
                common_edge = intersect(face1,face2);
                common_neighbourhood_vertex = common_edge(not( ismember(common_edge,1) ));
                Re_all_edge_vects(edge_id,:) = all_points_insequence(common_neighbourhood_vertex,:)...
                    - all_points_insequence(1,:);
                %                 dihedral_angle(edge_id) = -acos(dot(face1_normal,face2_normal))+pi;
                inward_or_outward_deformation_sign = sign( dot( cross(face1_normal,face2_normal), Re_all_edge_vects(edge_id,:) ) );
                dihedral_angle(edge_id) = inward_or_outward_deformation_sign*acos(dot(face1_normal,face2_normal))+pi;
                He_curvature(edge_id) = (real(2*vecnorm(Re_all_edge_vects(edge_id,:))*cos(dihedral_angle(edge_id)/2)));
                
                Se_operator_mat(:,:,edge_id) = He_curvature(edge_id)*[cross( Re_all_edge_vects(edge_id,:)...
                    /vecnorm(Re_all_edge_vects(edge_id,:)),Ne_all_edge_vects(edge_id,:) )]'*...
                    [cross( Re_all_edge_vects(edge_id,:)/vecnorm(Re_all_edge_vects(edge_id,:)),Ne_all_edge_vects(edge_id,:) )];
                
                We_weighting(edge_id,:) = dot(Normal_vect_obj , Ne_all_edge_vects(edge_id,:));
                Sv_update = Sv_update + We_weighting(edge_id,:)*Projection_operator_obj'*Se_operator_mat(:,:,edge_id)*...
                    Projection_operator_obj;
            end
            Av = sum(face_Area)/3;
            [Darboux_frame_eigVects,Darboux_frame_eigVals] = eig(Sv_update,'vector');
            % Arranging the darboux frame
            Darboux_frame_eigVals = Darboux_frame_eigVals';
            idx_ortho_noter = zeros(3,1);
            idx = (1:3);
            Darboux_frame_eigVects_arranged = NaN(3,3);
            Darboux_frame_eigVals_arranged = NaN(1,3);
            dot_product_noter = NaN(1,3);
            dot_pro = NaN(1,3);
            for eig_vec = 1:3
                eig_vector = Darboux_frame_eigVects(:,eig_vec);
                eig_vector_unit = eig_vector/vecnorm(eig_vector);
                dot_product_noter(eig_vec) = abs(dot(eig_vector_unit,Normal_vect_obj));
                dot_pro(eig_vec) =  dot(eig_vector_unit,Normal_vect_obj);
                if abs(dot(eig_vector_unit,Normal_vect_obj))>=0.9
                    idx_ortho_noter(eig_vec) = 1;
                end
            end
            %%
            if length( idx_ortho_noter(idx_ortho_noter==1)  )==1
                Darboux_frame_eigVects_arranged(:,3) = Darboux_frame_eigVects(:,idx_ortho_noter==1);
                Darboux_frame_eigVals_arranged(:,3) = Darboux_frame_eigVals(:,idx_ortho_noter==1);
                idx_remaining = idx(idx_ortho_noter==0);
                [in_plane_eigenVals,id] = sort(Darboux_frame_eigVals(idx_remaining));
                Darboux_frame_eigVals_arranged(:,1:2) = in_plane_eigenVals;
                Darboux_frame_eigVects_remaining = Darboux_frame_eigVects(:,idx_remaining);
                Darboux_frame_eigVects_arranged(:,1:2) = Darboux_frame_eigVects_remaining(:,id);
                %                 obj.DarbouxFrame{1,1} = Darboux_frame_eigVects_arranged;
                %                 obj.DarbouxFrame{2,1} = Darboux_frame_eigVals_arranged;
            else
                fprintf('Anomaly Detected\n');
                Darboux_frame_eigVects
                Darboux_frame_eigVals
                dot_pro
                Darboux_frame_eigVects_arranged(:,3) = obj.Normal_vect  ;
                Darboux_frame_eigVals_arranged(:,3) = 0;
                Darboux_frame_eigVals_arranged(:,1:2) = [0,0];
                Darboux_frame_eigVects_arranged(:,1:2) = NaN(3,2);
            end
            principle_curvatures = Darboux_frame_eigVals_arranged(:,1:2);
            new_energy_stored = 0.5*Av*obj.bending_rigidity*( (principle_curvatures(1)+principle_curvatures(2))/2-obj.spontaneous_curvature)^2;
            pressure_displacement =  mean(measure_distance_from_plane(obj.Pos,Normal_vect_obj,obj.neighbours_coor));
            projected_face_area = sum(calculate_projected_area_of_faces(Normal_vect_obj,face_normals,face_Area))/3;
            volume_of_push = 1/3*projected_face_area*pressure_displacement;
            edge_vect_mag = vecnorm(Re_all_edge_vects,2,2);
            if isempty(edge_vect_mag(edge_vect_mag<=2*hardball_radius))==0
                hard_ball_potential_obj = length(edge_vect_mag(edge_vect_mag<=2*hardball_radius))*overlapping_energy;
            else
                hard_ball_potential_obj = 0;
            end
            total_energy_new = hard_ball_potential_obj + new_energy_stored;
            % updating the local properties to updated obj
            updated_obj = obj;
            updated_obj.faceNormals = face_normals;
            updated_obj.faceArea = face_Area;
            updated_obj.Av_vertex = Av;
            updated_obj.Av_projected_vertex = projected_face_area;
            updated_obj.Normal_vect = Normal_vect_obj;
            updated_obj.Householder_mat;
            updated_obj.Projection_operator = Projection_operator_obj;
            
            updated_obj.dihedralAngle = dihedral_angle;
            updated_obj.Re = Re_all_edge_vects;
            updated_obj.He = He_curvature;
            updated_obj.Ne = Ne_all_edge_vects;
            updated_obj.Se_operators = Se_operator_mat;
            updated_obj.We = We_weighting;
            updated_obj.Sv_operator = Sv_update;
            updated_obj.DarbouxFrame{1,1} = Darboux_frame_eigVects_arranged;
            updated_obj.DarbouxFrame{2,1} = Darboux_frame_eigVals_arranged;
            
            %             updated_obj.energy_stored = new_energy_stored;
            %             updated_obj.work_volume = volume_of_push*sign(mean(Darboux_frame_eigVals_arranged));
            %             updated_obj.work_area = +Av - projected_face_area;
            %             updated_obj.hard_ball_potential = hard_ball_potential_obj;
            
            %             del_energy = total_energy_new - obj.hard_ball_potential-obj.energy_stored;
        end
        
        function [free_energy_total,energy_curvature,energy_surface_tension_normalized,energy_stored_by_pressure_outward,updated_obj] = ...
                calculate_energy_from_geometry(obj,surface_tension,hardball_radius,pressure)
            % Energy due to bending rigidity
            % <0.5*bending_rigidity*(H-H0)^2>
            H_mean = mean(obj.DarbouxFrame{2,1}(1:2));
            principle_curv =  obj.DarbouxFrame{2,1}(1:2);
            bending_rigidity_obj = obj.bending_rigidity;
            gauusian_rigidity_obj = obj.gaussian_rigidity;
            H_spontaneous = obj.spontaneous_curvature;
            energy_curvature = 0.5*obj.Av_vertex*bending_rigidity_obj*(2*H_mean - H_spontaneous)^2 + ...
                0.5*gauusian_rigidity_obj*( obj.DarbouxFrame{2,1}(1)*obj.DarbouxFrame{2,1}(2) - (0*H_spontaneous/2).^2 ).^2;
            %             energy_curvature = 0.5*obj.Av_vertex*bending_rigidity_obj*((principle_curv(1) - H_spontaneous)^2+(principle_curv(2) - H_spontaneous)^2);
            faces_objs = obj.faces;
            ID_obj = obj.ID;
            neighbours_obj = obj.neighbours;
            neighbours_coor_obj = obj.neighbours_coor;
            
            ID_list = [ID_obj,neighbours_obj'];
            ID_coor_list = [obj.Pos;neighbours_coor_obj];
            %             local_id_list = 1:size(ID_coor_list,1);
            
            faces_objs_local =  mask_mat_to_sequence(faces_objs,ID_list);
            
            energy_surface_tension_normalized = surface_tension*obj.Av_vertex - 0*hardball_radius;
            % free energy due to pressure
            distance_of_push = mean(measure_distance_from_plane(obj.Pos,obj.Normal_vect,obj.neighbours_coor));
            projected_face_area = sum(abs(calculate_projected_area_of_faces(obj.Normal_vect,obj.faceNormals,obj.faceArea)))/3;
            energy_stored_by_pressure_outward = pressure*distance_of_push*abs(projected_face_area);
            free_energy_total = -energy_stored_by_pressure_outward + energy_surface_tension_normalized + energy_curvature;
            obj.free_energy = free_energy_total;
            obj.curvature_energy = energy_curvature;
            obj.local_surface_tension_energy = energy_surface_tension_normalized;
            updated_obj = obj;
            
        end
        
        function [free_energy_total,energy_curvature,energy_surface_tension_normalized,energy_stored_by_pressure_outward,updated_obj] = ...
                calculate_energy_from_geometry_obj_coor_list_updated(obj,all_points,all_IDs,surface_tension,hardball_radius,pressure)
            neighbour_IDs = obj.neighbours;
            [~,locb] = ismember(neighbour_IDs,all_IDs);
            obj.neighbours_coor = all_points(locb,:);
            [updated_obj] =  calculate_topological_characteristics(obj,hardball_radius,100);
            [free_energy_total,energy_curvature,energy_surface_tension_normalized,energy_stored_by_pressure_outward,updated_obj] = ...
                calculate_energy_from_geometry(updated_obj,surface_tension,hardball_radius,pressure);
        end
        
        function [updated_obj_list] = update_neighbours_normals(obj_list)
            [all_IDs] = extract_IDs(obj_list);
            [normal_vects_all] = extract_Normal_vects(obj_list);
            updated_obj_list = obj_list;
            parfor obj_id = 1:size(obj_list,1)
                obj_current = obj_list(obj_id,1);
                obj_current_neighbours = obj_current.neighbours;
                obj_current_neighbours_normal = NaN(size(obj_current.neighbours_coor));
                for neighbour_ind = 1:length(obj_current_neighbours)
                    obj_current_neighbours_normal(neighbour_ind,:) = ...
                        normal_vects_all(all_IDs==obj_current_neighbours(neighbour_ind),:);
                end
                obj_current.neighbours_normal_old = obj_current_neighbours_normal;
                updated_obj_list(obj_id) = obj_current;
                fprintf('    >>>Neighbour normals updated %d\n',obj_current.ID);
            end
        end
        
        
        function [dihedral_angle_constraint, hardball_constraint,fold_mem_constraint] =...
                constraint_check(updated_obj,old_obj,min_dihedral_angles,max_dihedral_angles,hardball_rad)
            % Check Constraints--------------------------------------------
            % Check face constraint
            dihedral_angles = updated_obj.dihedralAngle;
            dihedral_angles_old = old_obj.dihedralAngle;
            dihedral_angle_constraint = min(abs(dihedral_angles))>=min_dihedral_angles & max(abs(dihedral_angles))<=max_dihedral_angles;
            %             dihedral_angle_constraint = dihedral_angle_constraint || ( min(abs(dihedral_angles))> min(abs(dihedral_angles_old)) && max(abs(dihedral_angles))< max(abs(dihedral_angles_old)) );
            
            
            % Check hardball/overlapping constarint
            neighbour_vects_mag = vecnorm(updated_obj.Re,2,2);
            neighbour_vects_mag_old = vecnorm(old_obj.Re,2,2);
            
            hardball_constraint = min(neighbour_vects_mag)>= min([2*hardball_rad,min(neighbour_vects_mag_old) ]) ...
                && max(neighbour_vects_mag)<= max([3^(.5)*2*hardball_rad max(neighbour_vects_mag_old)]);
            
            
            %             if min(neighbour_vects_mag)< 2*hardball_rad && max(neighbour_vects_mag)<3^(.5)*2*hardball_rad
            %                 hardball_constraint = (min(neighbour_vects_mag)>min(neighbour_vects_mag_old));
            %             end
            %
            %             if max(neighbour_vects_mag)>3^(.5)*2*hardball_rad && min(neighbour_vects_mag)< 2*hardball_rad
            %                 hardball_constraint = (max(neighbour_vects_mag)<max(neighbour_vects_mag_old));
            %             end
            
            normal_neighbours_old = updated_obj.neighbours_normal_old;
            updated_face_normals = updated_obj.faceNormals;
            updated_faces = updated_obj.faces;
            all_point_IDs = [updated_obj.ID,updated_obj.neighbours'];
            updated_faces_local = mask_mat_to_sequence(updated_faces,all_point_IDs);
            all_point_normal_vects = [updated_obj.Normal_vect;updated_obj.neighbours_normal_old];
            face_flipped_indicator = NaN(size(updated_faces,1),1);
            for face_id = 1:size(updated_faces,1)
                current_face = updated_faces_local(face_id,:);
                [~,idx] = find(current_face==1);
                current_face_arranged = circshift(current_face,1-idx);
                edge_vertex_normal_vect_mean = mean( all_point_normal_vects(current_face_arranged(:),:),1 );
                edge_vertex_normal_vect_mean = edge_vertex_normal_vect_mean./vecnorm(edge_vertex_normal_vect_mean);
                face_flipped_indicator(face_id) = (dot(all_point_normal_vects(1,:),edge_vertex_normal_vect_mean))<=0.2;
            end
            fold_mem_constraint = not(sum(face_flipped_indicator)>0);
            if fold_mem_constraint == 0
                disp('MEMBRANE FLIPPING EVENT FOUND\n\n');
                fprintf('DOT ALLpoints = %d\n\n\n',dot(all_point_normal_vects(1,:),edge_vertex_normal_vect_mean));
            end
        end
        
        
        
        %% new constraints
        %% local constraint
        function [area_constaint_factor] =...
                constraint_check_local(updated_obj,old_obj,area_reduced,k_a,temperature)
            %% New method
            triangles_updated_area = updated_obj.faceArea;
            area_constraint_updated = k_a*3*area_reduced*(triangles_updated_area/(3*area_reduced) - 1).^2;
            area_constraint_updated_sum = sum(area_constraint_updated);           
            
            triangles_old_area = old_obj.faceArea;
            area_constraint_old = k_a*3*area_reduced*(triangles_old_area/(3*area_reduced) - 1).^2;
            area_constraint_old_sum = sum(area_constraint_old);
            
            area_constaint_factor = exp(-(area_constraint_updated_sum-area_constraint_old_sum)/temperature);
            fprintf('   >>> FUNCTION:<constraint_check_local> area_constraint_factor = %d and hamiltonian_diff = %d\n',...
                area_constaint_factor,(area_constraint_updated_sum-area_constraint_old_sum));           
        end
        
        function [area_constaint_factor] =...
                constraint_check_local_leonard_jones(updated_obj,old_obj,area_reduced,k_a,temperature)
            %% New method
            triangles_updated_area = updated_obj.faceArea;
            area_constraint_updated = k_a*3*area_reduced*(triangles_updated_area/(3*area_reduced) - 1).^2;
            area_constraint_updated_sum = sum(area_constraint_updated);           
            
            triangles_old_area = old_obj.faceArea;
            area_constraint_old = k_a*3*area_reduced*(triangles_old_area/(3*area_reduced) - 1).^2;
            area_constraint_old_sum = sum(area_constraint_old);
            
            area_constaint_factor = exp(-(area_constraint_updated_sum-area_constraint_old_sum)/temperature);
            fprintf('   >>> FUNCTION:<constraint_check_local> area_constraint_factor = %d and hamiltonian_diff = %d\n',...
                area_constaint_factor,(area_constraint_updated_sum-area_constraint_old_sum));           
        end
        
        %% global constarint
        function [global_constaint_factor] =...
                constraint_check_global(updated_obj_list,old_obj_list,surface_tension,volume_reduced_total,K_V,temperature)
            
            [energy_curvature_sum_new,totalVolume_new,totalArea_new] = extract_energy_obj_list(updated_obj_list,surface_tension,1,1);
            [energy_curvature_sum_old,totalVolume_old,totalArea_old] = extract_energy_obj_list(old_obj_list,surface_tension,1,1);
            
            global_constraint_hamiltonian_new = surface_tension*totalArea_new + K_V*(totalVolume_new/volume_reduced_total-1)^2;
            global_constraint_hamiltonian_old = surface_tension*totalArea_old + K_V*(totalVolume_old/volume_reduced_total-1)^2;
            
            global_constaint_factor = min([0.9 ...
                exp(-(energy_curvature_sum_new-energy_curvature_sum_old)/temperature/size(old_obj_list,1))* ...
                exp(-(global_constraint_hamiltonian_new-global_constraint_hamiltonian_old)/temperature/size(old_obj_list,1))]);
            
            fprintf('   >>> FUNCTION:constraint_check_global global_constaint_factor = %d min(0.9 %d %d) \n',global_constaint_factor, ...
               exp(-(energy_curvature_sum_new-energy_curvature_sum_old)/temperature/size(old_obj_list,1)),...
               exp(-(global_constraint_hamiltonian_new-global_constraint_hamiltonian_old)/temperature/size(old_obj_list,1)));
            
        end
        
        %%
        function [updated_obj,object_updated_indicator,del_free_energy] = ...
                check_metropolis_MCSweep_new(updated_obj,Current_obj,surface_tension,hardball_radius,pressure,temperature,area_reduced,k_a)
            % Calculate Energy of previous and new object -----------------
            [~,energy_curvature_initial,~,~] = ...
                calculate_energy_from_geometry(Current_obj,surface_tension,hardball_radius,pressure);
            [~,energy_curvature_final,~,~] = ...
                calculate_energy_from_geometry(updated_obj,surface_tension,hardball_radius,pressure);
            
%             [dihedral_angle_constraint, hardball_constraint,fold_mem_constraint] = constraint_check(updated_obj,Current_obj,min_dihedral_angles,max_dihedral_angles,hardball_radius);
            [area_constaint_factor] =...
                constraint_check_local(updated_obj,Current_obj,area_reduced,k_a,temperature);
            % Look into the probability or acceptance rate
            del_free_energy = (energy_curvature_final) - (energy_curvature_initial); %% Since free energy must be always >=0 according to defination;
            acceptance_rate = min([1,exp(-del_free_energy/temperature)*area_constaint_factor]);
            
            object_updated_indicator = acceptance_rate>=1 ;
            if object_updated_indicator == 1
                %                 updated_obj = updated_obj;
                fprintf('    >>> FUNCTION: check_metropolis_MCSweep_new  ObjectID =%d del_curvature_constraint = %d && constraint_surface_tension = %d \n',...
                    updated_obj.ID,exp(-del_free_energy/temperature),area_constaint_factor);
            else
                updated_obj = Current_obj;
                fprintf('    >>>FUNCTION: check_metropolis_MCSweep_new ObjectID =%d NOT UPDATED del_curvature_constraint = %d && constraint_surface_tension = %d\n',...
                    updated_obj.ID,exp(-del_free_energy/temperature),area_constaint_factor);
            end
        end      
        
        function [updated_obj_list,object_updated_indicator,del_free_energy] = ...
                check_metropolis_MCSweep_obj_list(updated_obj_list,Current_obj_list,surface_tension,hardball_radius,pressure,min_dihedral_angles,max_dihedral_angles,temperature)
            % Calculate Energy of previous and new object -----------------
            del_free_energy = 0;
            dihedral_angle_constraint = 1;
            hardball_constraint = 1;
            fold_mem_constraint = 1;
            for ob_ind = 1:size(Current_obj_list,1)
                Current_obj = Current_obj_list(ob_ind,1);
                updated_obj = updated_obj_list(ob_ind,1);
                [free_energy_total_initial,energy_curvature_initial,energy_surface_tension_normalized_initial,energy_stored_by_pressure_outward_initial] = ...
                    calculate_energy_from_geometry(Current_obj,surface_tension,hardball_radius,pressure);
                [free_energy_total_final,energy_curvature_final,energy_surface_tension_normalized_final,energy_stored_by_pressure_outward_final] = ...
                    calculate_energy_from_geometry(updated_obj,surface_tension,hardball_radius,pressure);
                
                [dihedral_angle_constraint_current, hardball_constraint_current,fold_mem_constraint_current] = constraint_check(updated_obj,Current_obj,min_dihedral_angles,max_dihedral_angles,hardball_radius);
                
                % Look into the probability or acceptance rate
                del_free_energy = del_free_energy +((energy_curvature_final) - (energy_curvature_initial)); %% Since free energy must be always >=0 according to defination;
                %                 acceptance_rate = min(1,exp(-del_free_energy/temperature));
                dihedral_angle_constraint = dihedral_angle_constraint*dihedral_angle_constraint_current;
                hardball_constraint = hardball_constraint*hardball_constraint_current;
                fold_mem_constraint = fold_mem_constraint_current*fold_mem_constraint;
                
            end
            acceptance_rate = min(1,exp(-del_free_energy/temperature));
            object_updated_indicator = acceptance_rate>=1 & (dihedral_angle_constraint==1 & hardball_constraint==1 & fold_mem_constraint==1 );
            if object_updated_indicator == 1
                %                 updated_obj = updated_obj;
                fprintf('       >>>METROPOLIS Satisfied along with Constraints object list updated\n       >>>[delEnergy(%d) AcceptanceRate(%d)]\n',...
                    del_free_energy,acceptance_rate);
            else
                updated_obj_list = Current_obj_list;
                %                                 fprintf('       >>>METROPOLIS or CONSTRAINTS NOT Satisfied object %d  Not updated\n       >>>[delEnergy(%d) AcceptanceRate(%d) dihedralConstarint(%d) HardballConstarint(%d)]\n',...
                %                                     updated_obj.ID,del_free_energy,acceptance_rate,dihedral_angle_constraint,hardball_constraint);
            end
            
            
        end
        
        
        function [updated_obj_list,object_updated_indicator] = ...
                check_metropolis_MCSweep_obj_list_new(updated_obj_list,Current_obj_list,surface_tension,hardball_radius,pressure,temperature,area_reduced,k_a)
            % Calculate Energy of previous and new object -----------------
%             del_free_energy = 0;
%             dihedral_angle_constraint = 1;
%             hardball_constraint = 1;
%             fold_mem_constraint = 1;
            
            constraint_factor = 1;
            energy_metripolis_factor = 1;
            for ob_ind = 1:size(Current_obj_list,1)
                Current_obj = Current_obj_list(ob_ind,1);
                updated_obj = updated_obj_list(ob_ind,1);
                [~,energy_curvature_initial,~,~,Current_obj] = ...
                    calculate_energy_from_geometry(Current_obj,surface_tension,hardball_radius,pressure);
                [~,energy_curvature_final,~,~,updated_obj] = ...
                    calculate_energy_from_geometry(updated_obj,surface_tension,hardball_radius,pressure);
                
%                 [dihedral_angle_constraint_current, hardball_constraint_current,fold_mem_constraint_current] = constraint_check(updated_obj,Current_obj,min_dihedral_angles,max_dihedral_angles,hardball_radius);
                [area_constaint_factor] =   constraint_check_local(updated_obj,Current_obj,area_reduced,k_a,temperature);
                constraint_factor = constraint_factor*area_constaint_factor;
                energy_metripolis_factor = energy_metripolis_factor*exp(-(energy_curvature_final - energy_curvature_initial)/temperature);
                
                
%                 % Look into the probability or acceptance rate
%                 del_free_energy = del_free_energy +((energy_curvature_final) - (energy_curvature_initial)); %% Since free energy must be always >=0 according to defination;
%                 %                 acceptance_rate = min(1,exp(-del_free_energy/temperature));
%                 dihedral_angle_constraint = dihedral_angle_constraint*dihedral_angle_constraint_current;
%                 hardball_constraint = hardball_constraint*hardball_constraint_current;
%                 fold_mem_constraint = fold_mem_constraint_current*fold_mem_constraint;
                
            end
%             acceptance_rate = min(1,exp(-del_free_energy/temperature));
            acceptance_rate = min(1,energy_metripolis_factor*constraint_factor);
            object_updated_indicator = acceptance_rate;
            if object_updated_indicator == 1
                %                 updated_obj = updated_obj;
                fprintf('   >>> FUNCTION: check_metropolis_MCSweep_obj_list_new Object_list UPDATED based on accepetance score from energy %d and surface_constraint %d\n',...
                    energy_metripolis_factor,constraint_factor);
            else
                updated_obj_list = Current_obj_list;
                fprintf('   >>> FUNCTION: check_metropolis_MCSweep_obj_list_new Object_list NOT UPDATED based on accepetance score from energy %d and surface_constraint %d\n',...
                    energy_metripolis_factor,constraint_factor);
            end
            
            
        end
        
        
        % Degree of freedom: < LOCAL > Vertex move ----------       
        function [obj_list,acceptance_rate] = diffuse_fluctuation_orthogonal_instantaneous(obj_list,hardball_radius,kick_displacement,overlapping_energy,temperature_current,pressure,surface_tension,area_reduced,k_a)
            % Initialize
            fprintf('Initializing all IDS for diffusion\n');
            all_ID_extract_func = @(index) obj_list(index).ID;
            all_index = 1:size(obj_list,1);
            all_IDs = arrayfun(all_ID_extract_func,all_index);
            %---------------------------------------------------------
            % Randomize
            fprintf('Randomizing all IDS for diffusion\n');
            randomized_index = randperm(size(obj_list,1));
            randomized_IDs = all_IDs(randomized_index);
            %---------------------------------------------------------
            % Assign IDS for 3 degrees of freedom for mesh diffusion
            fprintf('Assigning all IDS for 3 Degs of diffusion\n');
            displacement_IDs = randomized_IDs(1:floor(end));
            %             vertex_flip_IDs = randomized_IDs(floor(end/3)+1:floor(2*end/3));
            %             Edge_flip_IDs = randomized_IDs(floor(2*end/3)+1:end);
            %---------------------------------------------------------
            % Update obj_list due to vertex displacement.
            fprintf('Diffusing\n');
            [displacement_IDs_objs,ind_in_list_displacement] = findobjID(obj_list,displacement_IDs);
            all_IDs = extract_IDs(displacement_IDs_objs);
            all_points = extract_coors_points(displacement_IDs_objs);
            extract_neighbours_of_current_obj = @(id) (obj_list(id).neighbours)';
            
            diffused_id_list_mask = zeros(size(displacement_IDs_objs,1),1);
            for id = 1:size(displacement_IDs_objs,1)
                current_ID = displacement_IDs(id);
                Current_obj = displacement_IDs_objs(id);
%                 displacement_vect = (rand(1,3)-0.5)*kick_displacement;
                displacement_vect = (rand(1,1)-0.5)*kick_displacement*Current_obj.Normal_vect;
                displacement_vect_mag = vecnorm(displacement_vect,2,2);
                updated_obj = Current_obj;
                updated_obj.Pos = updated_obj.Pos + displacement_vect;               
                
                % Checking metropolis
                [free_energy_total,energy_curvature,energy_surface_tension_normalized,energy_stored_by_pressure_outward,Current_obj] = ...
                    calculate_energy_from_geometry_obj_coor_list_updated(Current_obj,all_points,all_IDs,surface_tension,hardball_radius,pressure);
                [free_energy_total,energy_curvature,energy_surface_tension_normalized,energy_stored_by_pressure_outward,updated_obj] = ...
                    calculate_energy_from_geometry_obj_coor_list_updated(updated_obj,all_points,all_IDs,surface_tension,hardball_radius,pressure);                
                
                [updated_obj,can_move] = ...
                    check_metropolis_MCSweep_new(updated_obj,Current_obj,surface_tension,hardball_radius,pressure,temperature_current,area_reduced,k_a);
                
                
                
                if can_move
                    displacement_IDs_objs(id) = updated_obj;
                    all_points(id,:) = updated_obj.Pos;
                    diffused_id_list_mask(id) = 1;
                else
                    displacement_IDs_objs(id) = Current_obj;
                end
            end
            
            
            % Updating objects list
            obj_list(ind_in_list_displacement(:)) = displacement_IDs_objs;
            acceptance_rate = sum(diffused_id_list_mask)/length(diffused_id_list_mask)*100;
            fprintf('Acceptance Rate = %d percent\n',sum(diffused_id_list_mask)/length(diffused_id_list_mask)*100);
            
            
            
            %% updating neighbourhood normals
            [obj_list] = update_neighbours_normals(obj_list);
        end
        
        % Degree of freedom: < LOCAL > Kawashaki move ----------
        function [updated_objs] = vertex_flip(objs)
            %             old_obj1 = objs(1);
            %             old_obj2 = objs(2);
            updated_objs = objs;
            %  interchange ID
            %             updated_objs(1).ID = objs(2).ID;
            %             updated_objs(2).ID = objs(1).ID;
            
            %  interchange type
            updated_objs(1).type = objs(2).type;
            updated_objs(2).type = objs(1).type;
            
            
            %  interchange spontaneous_curvature
            updated_objs(1).spontaneous_curvature = objs(2).spontaneous_curvature;
            updated_objs(2).spontaneous_curvature = objs(1).spontaneous_curvature;
            
            %  interchange bending_rigidity
            updated_objs(1).bending_rigidity = objs(2).bending_rigidity;
            updated_objs(2).bending_rigidity = objs(1).bending_rigidity;            
            
            %
            
            % Calculating the energy
            principle_curvatures1 = updated_objs(1).DarbouxFrame{2,1};
            new_energy_stored1 = 0.5*updated_objs(1).Av_vertex*updated_objs(1).bending_rigidity*(principle_curvatures1(1)-updated_objs(1).spontaneous_curvature)^2+...
                0.5*updated_objs(1).Av_vertex*updated_objs(1).bending_rigidity*(principle_curvatures1(2)-updated_objs(1).spontaneous_curvature)^2;
            updated_objs(1).energy_stored = new_energy_stored1;
            
            principle_curvatures2 = updated_objs(2).DarbouxFrame{2,1};
            new_energy_stored2 = 0.5*updated_objs(2).Av_vertex*updated_objs(2).bending_rigidity*(principle_curvatures2(2)-updated_objs(2).spontaneous_curvature)^2+...
                0.5*updated_objs(2).Av_vertex*updated_objs(2).bending_rigidity*(principle_curvatures2(2)-updated_objs(2).spontaneous_curvature)^2;
            updated_objs(2).energy_stored = new_energy_stored2;
            
        end
        
        function [obj_list,acceptance_rate] = diffuse_vertex_flip(obj_list,temperature_current,hardball_radius,pressure,surface_tension,area_reduced,k_a)
            % Randomly selecting edges------------------------------------
            fprintf('Randomly selecting edges for vertex flip\n');
            [all_IDs] = extract_IDs(obj_list);
            randomized_index = randperm(size(obj_list,1));
            randomized_IDs = all_IDs(randomized_index);
            vertex_flip_obj_list = findobjID(obj_list,randomized_IDs(1:(size(obj_list,1))));
            [edge_to_be_flipped] = randomized_Partner_from_Neighbour(vertex_flip_obj_list);
            fprintf('Making obj list pair\n');
            
            
            
            [objs_for_vertex_flip_from_edgesCol1,edge_to_be_flipped_pos_in_obj_listcol1] = findobjID(obj_list,edge_to_be_flipped(:,1));
            [objs_for_vertex_flip_from_edgesCol1_IDs] = extract_IDs(objs_for_vertex_flip_from_edgesCol1);
            bi_pair_indices_inobj_list = 1:size(objs_for_vertex_flip_from_edgesCol1,1);
            find_indices_in_array = @(ID) bi_pair_indices_inobj_list( objs_for_vertex_flip_from_edgesCol1_IDs==ID  );
            pos_in_bi_pair_row = arrayfun(find_indices_in_array,edge_to_be_flipped(:,1));
            objs_for_vertex_flip_from_edgesCol1 = objs_for_vertex_flip_from_edgesCol1(pos_in_bi_pair_row);
            edge_to_be_flipped_pos_in_obj_listcol1 = edge_to_be_flipped_pos_in_obj_listcol1(pos_in_bi_pair_row);
            
            [objs_for_vertex_flip_from_edgesCol2,edge_to_be_flipped_pos_in_obj_listcol2] = findobjID(obj_list,edge_to_be_flipped(:,2));
            [objs_for_vertex_flip_from_edgesCol2_IDs] = extract_IDs(objs_for_vertex_flip_from_edgesCol2);
            bi_pair_indices_inobj_list = 1:size(objs_for_vertex_flip_from_edgesCol2,1);
            find_indices_in_array = @(ID) bi_pair_indices_inobj_list( objs_for_vertex_flip_from_edgesCol2_IDs==ID  );
            pos_in_bi_pair_row = arrayfun(find_indices_in_array,edge_to_be_flipped(:,2));
            objs_for_vertex_flip_from_edgesCol2 = objs_for_vertex_flip_from_edgesCol2(pos_in_bi_pair_row);
            edge_to_be_flipped_pos_in_obj_listcol2 = edge_to_be_flipped_pos_in_obj_listcol2(pos_in_bi_pair_row);
            
            
            obj_list_pair = [objs_for_vertex_flip_from_edgesCol1,objs_for_vertex_flip_from_edgesCol2];
            [indices_col1] = edge_to_be_flipped_pos_in_obj_listcol1;%indices_in_objlist(obj_list,edge_to_be_flipped(:,1));
            [indices_col2] = edge_to_be_flipped_pos_in_obj_listcol2;%indices_in_objlist(obj_list,edge_to_be_flipped(:,2));
            % Do vertex flip ---------------------------------------------
            pairs_flipped_mask = zeros(size(obj_list_pair,1),1);
            fprintf('Started FLIPPING\n');
            %             flipping_rate_mask = zeros(size(obj_list_pair,1),1);
            parfor obj_pair = 1:size(obj_list_pair,1)
                current_obj_pair = obj_list_pair(obj_pair,:);
                [updated_objs_pair] = vertex_flip(current_obj_pair);
                
                [updated_objs_pair,can_move_obj_list] = ...
                    check_metropolis_MCSweep_obj_list_new(updated_objs_pair,current_obj_pair,surface_tension,hardball_radius,pressure,temperature_current,area_reduced,k_a);           
                
                if can_move_obj_list == 1
                    %                     updated_objs_pair = [updated_obj1,updated_obj2];
                    obj_list_pair(obj_pair,:) = updated_objs_pair;
                    fprintf('   >>>VERTEX FLIPPED ID %d  with ID %d with acceptance score %d \n',...
                        updated_objs_pair(1).ID,updated_objs_pair(2).ID,can_move_obj_list);
                    pairs_flipped_mask(obj_pair) = 1;
                    %                 else
                    %                     fprintf('%d\n',acceptence_probability);
                end
            end
            acceptance_rate = sum(pairs_flipped_mask )/length(pairs_flipped_mask)*100;
            fprintf('Flipping rate = %d percent\n', acceptance_rate);
            % update list -------------------------------------------
            fprintf('FLIPPING completed and Updating\n');
            obj_list( indices_col1(:) ) = obj_list_pair(:,1);
            obj_list( indices_col2(:) ) = obj_list_pair(:,2);
            [obj_list] = update_neighbours_normals(obj_list);
        end
        
        % Degree of freedom: < GLOBAL > Volume move (scaling move) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
        function [updated_obj_list] = ...
                enhance_structures_based_on_curvature...
                (obj_list,normal_displacement,surface_tension,hardball_radius,overlapping_energy,pressure,Temperature,volume_constraint_coefficient,V_initial)
            
            updated_obj_list = obj_list;
            free_energy_curvature = NaN(size(obj_list,1),1);
            %%update positions according to normal of all objects
            parfor obj_no = 1:size(obj_list,1)
                current_obj = obj_list(obj_no);
                normal_vect_current = current_obj.Normal_vect;
                current_obj_neighbours_coor = current_obj.neighbours_coor;
                current_obj_neighbours_normal = current_obj.neighbours_normal_old;
                mean_curvature = mean(current_obj.DarbouxFrame{2,1}(1:2));
                normal_displacement_rand = rand(1)*normal_displacement;
                if sign(mean_curvature)>=0
                    displacement_vect = normal_displacement_rand*normal_vect_current;
                else
                    displacement_vect = -normal_displacement_rand*normal_vect_current;
                end
                current_obj.Pos = current_obj.Pos + displacement_vect;
                updated_obj_list(obj_no) = current_obj;
                fprintf('   >>> Function (enhance_structures_based_on_curvature) %d particle displaced according to curvature\n',obj_no);
            end
            
            
            %%update neighbours coor of all objects
            [updated_obj_list] = update_neighbour_info(updated_obj_list);
            
            %%calculate energy of the whole object list
            [updated_obj_list] = update_geometry_obj_list(updated_obj_list,hardball_radius,overlapping_energy);
            %             [energy_curvature_list] = extract_energy_obj_list(updated_obj_list,surface_tension,hardball_radius,pressure);
            [updated_obj_list,accpetance_probability] = ...
                enhance_structure_metropolis(updated_obj_list,obj_list,surface_tension,hardball_radius,pressure,Temperature,volume_constraint_coefficient,V_initial);
        end
        
        function [updated_obj_list,global_constaint_factor] = ...
                enhance_structure_metropolis(updated_obj_list,obj_list,surface_tension,hardball_radius,pressure,Temperature,K_V,volume_reduced_total)
%             [energy_curvature_sum_old,totalVolume_old,totalArea_old] = extract_energy_obj_list(obj_list,surface_tension,hardball_radius,pressure);
%             [energy_curvature_sum_new,totalVolume_new,totalArea_new] = extract_energy_obj_list(updated_obj_list,surface_tension,hardball_radius,pressure);
%             
%             free_energy_updated = energy_curvature_sum_new + volume_constraint_coefficient*(totalVolume_new/V_initial -1)^2 + surface_tension*totalArea_new;
%             free_energy_old = energy_curvature_sum_old + volume_constraint_coefficient*(totalVolume_old/V_initial -1)^2 + surface_tension*totalArea_old;
%             del_hamiltonian =  free_energy_updated - free_energy_old;
            
            [global_constaint_factor] = ...
                constraint_check_global(updated_obj_list,obj_list,surface_tension,volume_reduced_total,K_V,Temperature);
            
%             accpetance_probability = min(.72, exp(-del_hamiltonian/Temperature/size(obj_list,1)));
            fprintf('   >>> Function (enhance_structure_metropolis) VOLUME MOVE metropolis performed with acceptance = %d\n',global_constaint_factor);
            if global_constaint_factor<.9
                updated_obj_list = obj_list;
            end
        end
        
        
        %% Degree of freedom: <Local> Alexander move or edge flip ----------
        function [randomized_tetra_pairsID] = randomize_tetraPair(obj_list)
            all_IDs = extract_IDs(obj_list);
            %             all_ID_index_pos_in_obj_list = 1:size(obj_list,1);
            % Randomly select not connected tetra pairs
            randomized_id_index = randperm(size(obj_list,1));
            randomly_selected_objs_ID = all_IDs(randomized_id_index( 1:floor(size(obj_list,1)/8) ));
            [randomly_selected_objs,~] = findobjID(obj_list,randomly_selected_objs_ID);
            [opposite_pairs,edge_selected] = opposite_vertices_of_random_edge(randomly_selected_objs);
            
            tetra_pair = [edge_selected,opposite_pairs];
            % exclude overlapping tetra pairs or adjacent tetra pairs
            tetra_pair_unique = unique(tetra_pair,'rows');
            unique_vertex_list = unique(tetra_pair_unique(:));
            counts_each_unique_vertex =  histc(tetra_pair_unique(:),unique_vertex_list);
            %             adjacent_tetra_vertices = unique_vertex_list(counts_each_unique_vertex>1);
            tetra_pair_unique_valency = NaN(size(tetra_pair_unique));
            parfor tetra_pair_ind = 1:size(tetra_pair_unique,1)
                current_tetra_pair = sort(tetra_pair_unique(tetra_pair_ind,:));
                valencied_tetra_pair = counts_each_unique_vertex(ismember(unique_vertex_list,current_tetra_pair));
                tetra_pair_unique_valency(tetra_pair_ind,:) = valencied_tetra_pair;
            end
            % refine tetra_pair_unique for valency 4
            valency_4_containing_pairs_mask = tetra_pair_unique_valency(:,1)>=4|tetra_pair_unique_valency(:,2)>=4 ...
                |tetra_pair_unique_valency(:,3)>=4 | tetra_pair_unique_valency(:,4)>=4;
            tetra_pair_unique( valency_4_containing_pairs_mask ,: ) = [];
            
            unique_vertex_list = unique(tetra_pair_unique(:));
            counts_each_unique_vertex =  histc(tetra_pair_unique(:),unique_vertex_list);
            %             adjacent_tetra_vertices = unique_vertex_list(counts_each_unique_vertex>1);
            tetra_pair_unique_valency = NaN(size(tetra_pair_unique));
            parfor tetra_pair_ind = 1:size(tetra_pair_unique,1)
                current_tetra_pair = sort(tetra_pair_unique(tetra_pair_ind,:));
                valencied_tetra_pair = counts_each_unique_vertex(ismember(unique_vertex_list,current_tetra_pair));
                tetra_pair_unique_valency(tetra_pair_ind,:) = valencied_tetra_pair;
            end
            
            % refine tetra_pair_unique for valency 3
            valency_3_containing_pairs_mask = tetra_pair_unique_valency(:,1)==3|tetra_pair_unique_valency(:,2)==3 ...
                |tetra_pair_unique_valency(:,3)==3 | tetra_pair_unique_valency(:,4)==3;
            tetra_pair_unique( valency_3_containing_pairs_mask ,: ) = [];
            
            unique_vertex_list = unique(tetra_pair_unique(:));
            counts_each_unique_vertex =  histc(tetra_pair_unique(:),unique_vertex_list);
            %             adjacent_tetra_vertices = unique_vertex_list(counts_each_unique_vertex>1);
            tetra_pair_unique_valency = NaN(size(tetra_pair_unique));
            parfor tetra_pair_ind = 1:size(tetra_pair_unique,1)
                current_tetra_pair = sort(tetra_pair_unique(tetra_pair_ind,:));
                valencied_tetra_pair = counts_each_unique_vertex(ismember(unique_vertex_list,current_tetra_pair));
                tetra_pair_unique_valency(tetra_pair_ind,:) = valencied_tetra_pair;
            end
            
            % refine tetra_pair_unique for valency 2
            valency_3_containing_pairs_mask = tetra_pair_unique_valency(:,1)==2|tetra_pair_unique_valency(:,2)==2 ...
                |tetra_pair_unique_valency(:,3)==2 | tetra_pair_unique_valency(:,4)==2;
            tetra_pair_unique( valency_3_containing_pairs_mask ,: ) = [];
            randomized_tetra_pairsID = tetra_pair_unique;
        end
        
        function [decision_to_update] = angle_constraint_tetrapair(tetra_pair_obj)
            edge_to_remove_objs = tetra_pair_obj(1:2);
            edge_to_form_objs = tetra_pair_obj(3:4);
            % check valency to remove pair
            edge_to_rm_obj1 = edge_to_remove_objs(1);
            edge_to_rm_obj2 = edge_to_remove_objs(2);
            
            edge_to_rm_obj1_neighbours = edge_to_rm_obj1.neighbours;
            edge_to_rm_obj2_neighbours = edge_to_rm_obj2.neighbours;
            valencies_to_rm_objs = [length(edge_to_rm_obj1_neighbours);length(edge_to_rm_obj2_neighbours)];
            
            edge_to_form_obj1_neighbours = edge_to_form_objs(1).neighbours;
            edge_to_form_obj2_neighbours = edge_to_form_objs(2).neighbours;
            valencies_to_form_objs = [length(edge_to_form_obj1_neighbours);length(edge_to_form_obj2_neighbours)];
            
            if isempty(valencies_to_rm_objs(valencies_to_rm_objs<=3))==1 && isempty(valencies_to_form_objs(valencies_to_form_objs>=7))==1
                %% find the interesting plane (angle bisector of the normal vects of the remove points)
                normal_vect_plane = edge_to_form_objs(1).Normal_vect + edge_to_form_objs(2).Normal_vect +...
                    edge_to_remove_objs(1).Normal_vect + edge_to_remove_objs(2).Normal_vect;
                normal_vect_plane = normal_vect_plane/vecnorm(normal_vect_plane);
                
                
                
                %% Check projected angle constraint
                interseting_new_triangles = [ edge_to_form_objs(1).ID edge_to_rm_obj1.ID edge_to_rm_obj2.ID; ...
                    edge_to_form_objs(2).ID edge_to_rm_obj2.ID edge_to_rm_obj1.ID;];
                
                interesting_tri1_points =  [edge_to_form_objs(1).Pos;...
                    edge_to_rm_obj1.Pos;...
                    edge_to_rm_obj2.Pos];
                
                interesting_tri2_points =  [edge_to_form_objs(2).Pos;...
                    edge_to_rm_obj2.Pos;...
                    edge_to_rm_obj1.Pos];
                
                interesting_edge_vects_tri1 = [interesting_tri1_points(2,:)-interesting_tri1_points(1,:);...
                    interesting_tri1_points(3,:)-interesting_tri1_points(1,:)];
                component_along_normal =  diag(interesting_edge_vects_tri1*(normal_vect_plane'*[1,1]));
                interesting_edge_vects_tri1_transverse = interesting_edge_vects_tri1 - ...
                    [component_along_normal(1)*normal_vect_plane;...
                    component_along_normal(2)*normal_vect_plane];
                interesting_edge_vects_tri1_transverse_mag = vecnorm(interesting_edge_vects_tri1_transverse,2,2);
                interesting_edge_vects_tri1_transverse_unit = [interesting_edge_vects_tri1_transverse(1,:)/interesting_edge_vects_tri1_transverse_mag(1);...
                    interesting_edge_vects_tri1_transverse(2,:)/interesting_edge_vects_tri1_transverse_mag(2)];
                
                
                interesting_edge_vects_tri2 = [interesting_tri2_points(2,:)-interesting_tri2_points(1,:);...
                    interesting_tri2_points(3,:)-interesting_tri2_points(1,:)];
                component_along_normal =  diag(interesting_edge_vects_tri2*(normal_vect_plane'*[1,1]));
                interesting_edge_vects_tri2_transverse = interesting_edge_vects_tri2 - ...
                    [component_along_normal(1)*normal_vect_plane;...
                    component_along_normal(2)*normal_vect_plane];
                interesting_edge_vects_tri2_transverse_mag = vecnorm(interesting_edge_vects_tri2_transverse,2,2);
                interesting_edge_vects_tri2_transverse_unit = [interesting_edge_vects_tri2_transverse(1,:)/interesting_edge_vects_tri2_transverse_mag(1);...
                    interesting_edge_vects_tri2_transverse(2,:)/interesting_edge_vects_tri2_transverse_mag(2)];
                
                % calculate angle with new edge formed with edges in
                % projected plane
                angle_bisector_tri1 = interesting_edge_vects_tri1_transverse_unit(1,:)+interesting_edge_vects_tri1_transverse_unit(2,:);
                angle_bisector_tri1_unit = angle_bisector_tri1/vecnorm(angle_bisector_tri1);
                angle_bisector_tri2 = interesting_edge_vects_tri2_transverse_unit(1,:)+interesting_edge_vects_tri2_transverse_unit(2,:);
                angle_bisector_tri2_unit = angle_bisector_tri2/vecnorm(angle_bisector_tri2);
                
                angle_tri1 = acos(dot(interesting_edge_vects_tri1_transverse_unit(1,:), interesting_edge_vects_tri1_transverse_unit(2,:)))*180/pi;
                angle_tri2 = acos(dot(interesting_edge_vects_tri2_transverse_unit(1,:), interesting_edge_vects_tri2_transverse_unit(2,:)))*180/pi;
                
                new_edge_vect = edge_to_form_objs(2).Pos - edge_to_form_objs(1).Pos;
                new_edge_vect_projected = new_edge_vect - dot(new_edge_vect,normal_vect_plane)*normal_vect_plane;
                new_edge_vect_projected_unit = new_edge_vect_projected/vecnorm(new_edge_vect_projected);
                
                new_edge_vect_unit_fromtri1 =  new_edge_vect_projected_unit;
                angle_with_bisector_tri1 = acos(dot( new_edge_vect_unit_fromtri1,angle_bisector_tri1_unit ))*180/pi;
                angle_with_bisector_tri2 = acos(dot( -new_edge_vect_unit_fromtri1,angle_bisector_tri2_unit ))*180/pi;
                
                if angle_with_bisector_tri1<angle_tri1/3 && angle_with_bisector_tri2<angle_tri2/3
                    decision_to_update = 1;
                else
                    decision_to_update = 0;
                end
                
                
            else
                decision_to_update = 0;
                fprintf('Valency 3 tetra pair found!\n');
            end
            
            edge_to_rm_obj1_ = edge_to_rm_obj1.neighbours;
            edge_to_rm_obj2_neighbours = edge_to_rm_obj2.neighbours;
            
            
            
            
        end
        
        function [updated_tetra_pair_obj_row,updated_indicator] = update_tetra_pair(tetra_pair_obj,hardball_radius,overlapping_energy)
            % randomized_tetra_pairsID is arranged in [minEDGElive
            % maxEDGElive minEDGEformed maxEDGEformed] and each tetra pair
            % is mutually exclusive
            %             updated_tetra_pair_obj_row = tetra_pair_obj;
            edge_to_remove_objs = tetra_pair_obj(1:2);
            edge_to_form_objs = tetra_pair_obj(3:4);
            length_of_deleted_edge = vecnorm(edge_to_remove_objs(1).Pos - edge_to_remove_objs(2).Pos);
            length_of_form_edge = vecnorm(edge_to_form_objs(1).Pos - edge_to_form_objs(2).Pos);
            %             decision_to_update_length = length_of_form_edge<=length_of_deleted_edge;
            
            decision_to_update_length = length_of_form_edge<=max([3^(.5)*2*hardball_radius length_of_deleted_edge])...
                && length_of_form_edge>=min([2*hardball_radius length_of_deleted_edge]);
            % check before flipping if in new edge either of the old edge
            % fall
            edge_to_remove_objs1_pos = edge_to_remove_objs(1).Pos;
            edge_to_remove_objs2_pos = edge_to_remove_objs(2).Pos;
            
            edge_to_del_vect =  edge_to_remove_objs2_pos-edge_to_remove_objs1_pos;
            edge_to_del_vect_unit = edge_to_del_vect/vecnorm(edge_to_del_vect);
            
            edge_to_form_objs1_pos = edge_to_form_objs(1).Pos;
            edge_to_form_objs2_pos = edge_to_form_objs(2).Pos;
            
            edge_to_form_vect =  edge_to_form_objs2_pos-edge_to_form_objs1_pos;
            edge_to_form_vect_unit = edge_to_form_vect./vecnorm(edge_to_form_vect);
            
            edge_vect_triangle1 = edge_to_remove_objs1_pos - edge_to_form_objs1_pos;
            edge_vect_triangle1 = edge_vect_triangle1./vecnorm(edge_vect_triangle1);
            
            edge_vect_triangle2 = edge_to_remove_objs2_pos - edge_to_form_objs1_pos;
            edge_vect_triangle2 = edge_vect_triangle2./vecnorm(edge_vect_triangle2);
            
            edge_vect_triangle3 = edge_to_remove_objs1_pos - edge_to_form_objs2_pos;
            edge_vect_triangle3 = edge_vect_triangle3./vecnorm(edge_vect_triangle3);
            
            edge_vect_triangle4 = edge_to_remove_objs2_pos - edge_to_form_objs2_pos;
            edge_vect_triangle4 = edge_vect_triangle4./vecnorm(edge_vect_triangle4);
            
            vects_to_check_dot_product = [edge_vect_triangle1;...
                edge_vect_triangle2;...
                edge_vect_triangle3;...
                edge_vect_triangle4];
            idx_scalerProductequal1_indicator = NaN(4,1);
            for vect_to_check_ind =1:4
                idx_scalerProductequal1_indicator(vect_to_check_ind) = abs(dot(edge_to_form_vect_unit,vects_to_check_dot_product(vect_to_check_ind,:)));
            end
            idx_scalerProductequal1_indicator(idx_scalerProductequal1_indicator<0.95) = 0;
            decision_to_update_scaler = isempty(idx_scalerProductequal1_indicator(idx_scalerProductequal1_indicator~=0));
            
            
            % check if the angle between existing edges is more twice then the
            % angle between the angle bisector and the new edge
            angle_between_existing_edge_tri1 = dot(edge_vect_triangle1,edge_vect_triangle2);
            angle_bisector_vect1 = (edge_vect_triangle1+edge_vect_triangle2)/1;
            angle_bisector_vect1 = angle_bisector_vect1/vecnorm(angle_bisector_vect1);
            angle_tri1 = acos(angle_between_existing_edge_tri1);
            angle_with_edge_bisector_tri1 = acos((dot(angle_bisector_vect1,edge_to_form_vect_unit)));
            
            angle_between_existing_edge_tri2 = dot(edge_vect_triangle3,edge_vect_triangle4);
            angle_bisector_vect2 = (edge_vect_triangle3 + edge_vect_triangle4)/1;
            angle_bisector_vect2 = angle_bisector_vect2/vecnorm(angle_bisector_vect2);
            angle_tri2 = acos(angle_between_existing_edge_tri2);
            angle_with_edge_bisector_tri2 = acos((dot(angle_bisector_vect2,-edge_to_form_vect_unit)));
            
            %             decision_to_update = ...
            %                 angle_with_edge_bisector_tri1<angle_tri1/2*2/3 & angle_with_edge_bisector_tri2<angle_tri2/2*2/3 & length_of_form_edge <1.5*length_of_deleted_edge;
            [decision_to_update] = angle_constraint_tetrapair(tetra_pair_obj);
            
            %% check the projected surface angle constraints
            
            
            decision_to_update = decision_to_update & decision_to_update_length;
            if decision_to_update==1
                
                % Updating minLIVE --------------------------------------------
                edge_to_remove_objs1 = edge_to_remove_objs(1);
                % update the neighbours and neighbours coor
                old_neighbours_ob1 = edge_to_remove_objs1.neighbours;
                old_neighbours_ob1_coor = edge_to_remove_objs1.neighbours_coor;
                neighbourID_toDel = edge_to_remove_objs(2).ID;
                updated_neighbours_ob1 = old_neighbours_ob1;
                updated_neighbours_ob1(old_neighbours_ob1==neighbourID_toDel) = [];
                updated_neighbours_ob1_coor = old_neighbours_ob1_coor;
                updated_neighbours_ob1_coor(old_neighbours_ob1==neighbourID_toDel,:) = [];
                
                % update faces
                old_faces_ob1 = edge_to_remove_objs1.faces;
                edge_to_del = [edge_to_remove_objs(1).ID,edge_to_remove_objs(2).ID];
                faces_containg_edge_toDel_mask = sum(double(ismember(old_faces_ob1,edge_to_del)),2)==2;
                faces_containg_edge_toDel = old_faces_ob1(faces_containg_edge_toDel_mask,:);
                new_edge = unique(faces_containg_edge_toDel(not(ismember(faces_containg_edge_toDel,edge_to_del))));
                new_face = [edge_to_remove_objs1.ID,new_edge'];
                
                % check face order if alligned in the direction of normal vect
                new_face_vertex_coor = [edge_to_remove_objs1.Pos;...
                    updated_neighbours_ob1_coor(updated_neighbours_ob1==new_face(2),:);...
                    updated_neighbours_ob1_coor(updated_neighbours_ob1==new_face(3),:)];
                new_face_edge_vect21 = new_face_vertex_coor(2,:)-new_face_vertex_coor(1,:);
                new_face_edge_vect32 = new_face_vertex_coor(3,:)-new_face_vertex_coor(2,:);
                face_normal_newface = cross(new_face_edge_vect21,new_face_edge_vect32);
                if sign(dot(edge_to_remove_objs1.Normal_vect,face_normal_newface))<0
                    new_face = fliplr(new_face);
                end
                updated_faces = old_faces_ob1;
                updated_faces(faces_containg_edge_toDel_mask,:) = [];
                updated_faces(end+1,:) = new_face;
                %                 [updated_face_pairs,anomaly1] = update_facePair(edge_to_remove_objs1,updated_faces);
                
                
                updated_edge_to_remove_objs1 = edge_to_remove_objs1;
                updated_edge_to_remove_objs1.neighbours = updated_neighbours_ob1;
                updated_edge_to_remove_objs1.neighbours_coor = updated_neighbours_ob1_coor;
                updated_edge_to_remove_objs1.faces = updated_faces;
                [updated_face_pairs,anomaly1] = update_facePair(updated_edge_to_remove_objs1,updated_faces);
                updated_edge_to_remove_objs1.facePairs = updated_face_pairs;
                % -------------------------------------------------------------
                
                
                
                % Updating maxLIVE --------------------------------------------
                edge_to_remove_objs2 = edge_to_remove_objs(2);
                % update the neighbours and neighbours coor
                old_neighbours_ob2 = edge_to_remove_objs2.neighbours;
                old_neighbours_ob2_coor = edge_to_remove_objs2.neighbours_coor;
                neighbourID_toDel = edge_to_remove_objs(1).ID;
                updated_neighbours_ob2 = old_neighbours_ob2;
                updated_neighbours_ob2(old_neighbours_ob2==neighbourID_toDel) = [];
                updated_neighbours_ob2_coor = old_neighbours_ob2_coor;
                updated_neighbours_ob2_coor(old_neighbours_ob2==neighbourID_toDel,:) = [];
                
                % update faces
                old_faces_ob2 = edge_to_remove_objs2.faces;
                edge_to_del = [edge_to_remove_objs(2).ID,edge_to_remove_objs(1).ID];
                faces_containg_edge_toDel_mask = sum(double(ismember(old_faces_ob2,edge_to_del)),2)==2;
                faces_containg_edge_toDel = old_faces_ob2(faces_containg_edge_toDel_mask,:);
                new_edge = unique(faces_containg_edge_toDel(not(ismember(faces_containg_edge_toDel,edge_to_del))));
                new_face = [edge_to_remove_objs2.ID,new_edge'];
                
                % check face order if alligned in the direction of normal vect
                new_face_vertex_coor = [edge_to_remove_objs2.Pos;...
                    updated_neighbours_ob2_coor(updated_neighbours_ob2==new_face(2),:);...
                    updated_neighbours_ob2_coor(updated_neighbours_ob2==new_face(3),:)];
                new_face_edge_vect21 = new_face_vertex_coor(2,:)-new_face_vertex_coor(1,:);
                new_face_edge_vect32 = new_face_vertex_coor(3,:)-new_face_vertex_coor(2,:);
                face_normal_newface = cross(new_face_edge_vect21,new_face_edge_vect32);
                if sign(dot(edge_to_remove_objs2.Normal_vect,face_normal_newface))<0
                    new_face = fliplr(new_face);
                end
                updated_faces = old_faces_ob2;
                updated_faces(faces_containg_edge_toDel_mask,:) = [];
                updated_faces(end+1,:) = new_face;
                %                 [updated_face_pairs,anomaly2] = update_facePair(edge_to_remove_objs2,updated_faces);
                
                
                updated_edge_to_remove_objs2 = edge_to_remove_objs2;
                updated_edge_to_remove_objs2.neighbours = updated_neighbours_ob2;
                updated_edge_to_remove_objs2.neighbours_coor = updated_neighbours_ob2_coor;
                updated_edge_to_remove_objs2.faces = updated_faces;
                [updated_face_pairs,anomaly2] = update_facePair(updated_edge_to_remove_objs2,updated_faces);
                updated_edge_to_remove_objs2.facePairs = updated_face_pairs;
                % -------------------------------------------------------------
                
                
                % Updating minFORM --------------------------------------------
                edge_to_form_objs1 = edge_to_form_objs(1);
                % update the neighbours and neighbours coor
                old_neighbours_ob1 = edge_to_form_objs1.neighbours;
                old_neighbours_ob1_coor = edge_to_form_objs1.neighbours_coor;
                neighbourID_toConnect = edge_to_form_objs(2).ID;
                updated_neighbours_ob1 = old_neighbours_ob1;
                updated_neighbours_ob1(end+1) = neighbourID_toConnect;
                updated_neighbours_ob1_coor = old_neighbours_ob1_coor;
                updated_neighbours_ob1_coor(end + 1,:) = edge_to_form_objs(2).Pos;
                
                % update faces
                old_faces_ob1 = edge_to_form_objs1.faces;
                edge_to_del = [edge_to_remove_objs(1).ID,edge_to_remove_objs(2).ID];
                faces_containg_edge_toDel_mask = sum(double(ismember(old_faces_ob1,edge_to_del)),2)==2;
                faces_containg_edge_toDel = old_faces_ob1(faces_containg_edge_toDel_mask,:);
                edges_of_the_facetoDel = [faces_containg_edge_toDel(1),faces_containg_edge_toDel(2);...
                    faces_containg_edge_toDel(1),faces_containg_edge_toDel(3);...
                    faces_containg_edge_toDel(2),faces_containg_edge_toDel(3)];
                edges_of_face_from_current_vertex_mask = sum(double(ismember(edges_of_the_facetoDel,[edge_to_del])),2)==1;
                edges_of_face_from_current_vertex = edges_of_the_facetoDel( edges_of_face_from_current_vertex_mask,:);
                new_face1 = [edges_of_face_from_current_vertex(1,:),neighbourID_toConnect];
                new_face2 = [edges_of_face_from_current_vertex(2,:),neighbourID_toConnect];
                [~,new_face1_current_vertex_ind] = find(new_face1==edge_to_form_objs1.ID);
                new_face1_circshift = circshift(new_face1,1-new_face1_current_vertex_ind);
                new_face1 = new_face1_circshift;
                
                [~,new_face2_current_vertex_ind] = find(new_face2==edge_to_form_objs1.ID);
                new_face2_circshift = circshift(new_face2,1-new_face2_current_vertex_ind);
                new_face2 = new_face2_circshift;
                
                % check face order if alligned in the direction of normal vect
                new_face_vertex_coor = [edge_to_form_objs1.Pos;...
                    updated_neighbours_ob1_coor(updated_neighbours_ob1==new_face1(2),:);...
                    updated_neighbours_ob1_coor(updated_neighbours_ob1==new_face1(3),:)];
                
                new_face_edge_vect21 = new_face_vertex_coor(2,:)-new_face_vertex_coor(1,:);
                new_face_edge_vect32 = new_face_vertex_coor(3,:)-new_face_vertex_coor(2,:);
                face_normal_newface = cross(new_face_edge_vect21,new_face_edge_vect32);
                if sign(dot(edge_to_form_objs1.Normal_vect,face_normal_newface))<0
                    new_face1 = fliplr(new_face1);
                end
                
                new_face_vertex_coor = [edge_to_form_objs1.Pos;...
                    updated_neighbours_ob1_coor(updated_neighbours_ob1==new_face2(2),:);...
                    updated_neighbours_ob1_coor(updated_neighbours_ob1==new_face2(3),:)];
                
                new_face_edge_vect21 = new_face_vertex_coor(2,:)-new_face_vertex_coor(1,:);
                new_face_edge_vect32 = new_face_vertex_coor(3,:)-new_face_vertex_coor(2,:);
                face_normal_newface = cross(new_face_edge_vect21,new_face_edge_vect32);
                if sign(dot(edge_to_form_objs1.Normal_vect,face_normal_newface))<0
                    new_face2 = fliplr(new_face2);
                end
                updated_faces = old_faces_ob1;
                updated_faces(faces_containg_edge_toDel_mask,:) = [];
                updated_faces = [updated_faces;new_face1;new_face2];
                %                 [updated_face_pairs,anomaly3] = update_facePair(edge_to_form_objs1,updated_faces);
                
                updated_edge_to_form_objs1 = edge_to_form_objs1;
                updated_edge_to_form_objs1.neighbours = updated_neighbours_ob1;
                updated_edge_to_form_objs1.neighbours_coor = updated_neighbours_ob1_coor;
                updated_edge_to_form_objs1.faces = updated_faces;
                [updated_face_pairs,anomaly3] = update_facePair(updated_edge_to_form_objs1,updated_faces);
                updated_edge_to_form_objs1.facePairs = updated_face_pairs;
                % -------------------------------------------------------------
                
                % Updating maxFORM --------------------------------------------
                edge_to_form_objs2 = edge_to_form_objs(2);
                % update the neighbours and neighbours coor
                old_neighbours_ob2 = edge_to_form_objs2.neighbours;
                old_neighbours_ob2_coor = edge_to_form_objs2.neighbours_coor;
                neighbourID_toConnect = edge_to_form_objs(1).ID;
                updated_neighbours_ob2 = old_neighbours_ob2;
                updated_neighbours_ob2(end+1) = neighbourID_toConnect;
                updated_neighbours_ob2_coor = old_neighbours_ob2_coor;
                updated_neighbours_ob2_coor(end + 1,:) = edge_to_form_objs(1).Pos;
                
                % update faces
                old_faces_ob2 = edge_to_form_objs2.faces;
                edge_to_del = [edge_to_remove_objs(1).ID,edge_to_remove_objs(2).ID];
                faces_containg_edge_toDel_mask = sum(double(ismember(old_faces_ob2,edge_to_del)),2)==2;
                faces_containg_edge_toDel = old_faces_ob2(faces_containg_edge_toDel_mask,:);
                edges_of_the_facetoDel = [faces_containg_edge_toDel(1),faces_containg_edge_toDel(2);...
                    faces_containg_edge_toDel(1),faces_containg_edge_toDel(3);...
                    faces_containg_edge_toDel(2),faces_containg_edge_toDel(3)];
                edges_of_face_from_current_vertex_mask = sum(double(ismember(edges_of_the_facetoDel,[edge_to_del])),2)==1;
                edges_of_face_from_current_vertex = edges_of_the_facetoDel( edges_of_face_from_current_vertex_mask,:);
                new_face1 = [edges_of_face_from_current_vertex(1,:),neighbourID_toConnect];
                new_face2 = [edges_of_face_from_current_vertex(2,:),neighbourID_toConnect];
                [~,new_face1_current_vertex_ind] = find(new_face1==edge_to_form_objs2.ID);
                new_face1_circshift = circshift(new_face1,1-new_face1_current_vertex_ind);
                new_face1 = new_face1_circshift;
                
                [~,new_face2_current_vertex_ind] = find(new_face2==edge_to_form_objs2.ID);
                new_face2_circshift = circshift(new_face2,1-new_face2_current_vertex_ind);
                new_face2 = new_face2_circshift;
                
                % check face order if alligned in the direction of normal vect
                new_face_vertex_coor = [edge_to_form_objs2.Pos;...
                    updated_neighbours_ob2_coor(updated_neighbours_ob2==new_face1(2),:);...
                    updated_neighbours_ob2_coor(updated_neighbours_ob2==new_face1(3),:)];
                
                new_face_edge_vect21 = new_face_vertex_coor(2,:)-new_face_vertex_coor(1,:);
                new_face_edge_vect32 = new_face_vertex_coor(3,:)-new_face_vertex_coor(2,:);
                face_normal_newface = cross(new_face_edge_vect21,new_face_edge_vect32);
                if sign(dot(edge_to_form_objs2.Normal_vect,face_normal_newface))<0
                    new_face1 = fliplr(new_face1);
                end
                
                new_face_vertex_coor = [edge_to_form_objs2.Pos;...
                    updated_neighbours_ob2_coor(updated_neighbours_ob2==new_face2(2),:);...
                    updated_neighbours_ob2_coor(updated_neighbours_ob2==new_face2(3),:)];
                
                new_face_edge_vect21 = new_face_vertex_coor(2,:)-new_face_vertex_coor(1,:);
                new_face_edge_vect32 = new_face_vertex_coor(3,:)-new_face_vertex_coor(2,:);
                face_normal_newface = cross(new_face_edge_vect21,new_face_edge_vect32);
                if sign(dot(edge_to_form_objs2.Normal_vect,face_normal_newface))<0
                    new_face2 = fliplr(new_face2);
                end
                updated_faces = old_faces_ob2;
                updated_faces(faces_containg_edge_toDel_mask,:) = [];
                updated_faces = [updated_faces;new_face1;new_face2];
                %                 [updated_face_pairs,anomaly4] = update_facePair(edge_to_form_objs2,updated_faces);
                
                updated_edge_to_form_objs2 = edge_to_form_objs2;
                updated_edge_to_form_objs2.neighbours = updated_neighbours_ob2;
                updated_edge_to_form_objs2.neighbours_coor = updated_neighbours_ob2_coor;
                updated_edge_to_form_objs2.faces = updated_faces;
                [updated_face_pairs,anomaly4] = update_facePair(updated_edge_to_form_objs2,updated_faces);
                updated_edge_to_form_objs2.facePairs = updated_face_pairs;
                % -------------------------------------------------------------
                if not(anomaly1 | anomaly2 | anomaly3 | anomaly4)
                    [updated_edge_to_remove_objs1] = calculate_topological_characteristics(updated_edge_to_remove_objs1,hardball_radius,overlapping_energy);
                    [updated_edge_to_remove_objs2] = calculate_topological_characteristics(updated_edge_to_remove_objs2,hardball_radius,overlapping_energy);
                    [updated_edge_to_form_objs1] = calculate_topological_characteristics(updated_edge_to_form_objs1,hardball_radius,overlapping_energy);
                    [updated_edge_to_form_objs2] = calculate_topological_characteristics(updated_edge_to_form_objs2,hardball_radius,overlapping_energy);
                    
                    updated_tetra_pair_obj_row = [updated_edge_to_remove_objs1,updated_edge_to_remove_objs2,updated_edge_to_form_objs1,updated_edge_to_form_objs2];
                    %                     fprintf('%f<%f   %f<%f\n',...
                    %                         angle_with_edge_bisector_tri1*180/pi,angle_tri1*180/pi/2, angle_with_edge_bisector_tri2*180/pi,angle_tri2*180/pi/2);
                    
                    %                     idx_scalerProductequal1_indicator
                    %                     fprintf('updated tetra pair \n');
                    updated_indicator = 1;
                else
                    %                     [updated_edge_to_remove_objs1] = edge_to_remove_objs1;
                    %                     [updated_edge_to_remove_objs2] = edge_to_remove_objs2;
                    %                     [updated_edge_to_form_objs1] = edge_to_form_objs1;
                    %                     [updated_edge_to_form_objs2] = edge_to_form_objs2;
                    fprintf('ANOMALY Detected \n');
                    
                    updated_tetra_pair_obj_row = tetra_pair_obj;
                    updated_indicator = 0;
                end
            else
                updated_tetra_pair_obj_row = tetra_pair_obj;
                
                %                 fprintf('angle_with_edge_bisector_tri1(%d)<angle_tri1(%d)/2 & angle_with_edge_bisector_tri2(%d)<angle_tri2(%d)/2 & decision_scaler %d\n',...
                %                     angle_with_edge_bisector_tri1*180/pi,angle_tri1*180/pi, angle_with_edge_bisector_tri2*180/pi,angle_tri2*180/pi,decision_to_update_scaler);
                %
                %                 fprintf('%f<%f   %f<%f\n',...
                %                     angle_with_edge_bisector_tri1*180/pi,angle_tri1*180/pi/2, angle_with_edge_bisector_tri2*180/pi,angle_tri2*180/pi/2);
                %
                %                 idx_scalerProductequal1_indicator
                %                 fprintf('Same tetra pair due to geometric constraint\n');
                updated_indicator = 0;
            end
            
            %             if length_of_form_edge <1.5*length_of_deleted_edge
            %                 updated_indicator = 0;
            %             end
            
            
            
        end
        
        function [tetra_pair_obj_list_updated,updated_obj_list,Edge_flip_rate] = ...
                edge_flip_tetra_pair_obj_list(obj_list,randomized_tetra_pairsID,hardball_radius,overlapping_energy,temperature_current,pressure,surface_tension)
            % randomized_tetra_pairsID is arranged in [minEDGElive
            % maxEDGElive minEDGEformed maxEDGEformed] and each tetra pair
            % is mutually exclusive
            [minEDGELive_objs,pos_in_obj_list_minEDGElive] = findobjID(obj_list,randomized_tetra_pairsID(:,1));
            [minEDGELive_objs_IDs] = extract_IDs(minEDGELive_objs);
            tetra_pair_indices_inobj_list = 1:size(minEDGELive_objs,1);
            find_indices_in_array = @(ID) tetra_pair_indices_inobj_list( minEDGELive_objs_IDs==ID  );
            pos_in_tetra_pair_row = arrayfun(find_indices_in_array,randomized_tetra_pairsID(:,1));
            minEDGELive_objs = minEDGELive_objs(pos_in_tetra_pair_row);
            pos_in_obj_list_minEDGElive = pos_in_obj_list_minEDGElive(pos_in_tetra_pair_row);
            
            [maxEDGELive_objs,pos_in_obj_list_maxEDGElive] = findobjID(obj_list,randomized_tetra_pairsID(:,2));
            [maxEDGELive_objs_IDs] = extract_IDs(maxEDGELive_objs);
            tetra_pair_indices_inobj_list = 1:size(maxEDGELive_objs,1);
            find_indices_in_array = @(ID) tetra_pair_indices_inobj_list( maxEDGELive_objs_IDs==ID  );
            pos_in_tetra_pair_row = arrayfun(find_indices_in_array,randomized_tetra_pairsID(:,2));
            maxEDGELive_objs = maxEDGELive_objs(pos_in_tetra_pair_row);
            pos_in_obj_list_maxEDGElive = pos_in_obj_list_maxEDGElive(pos_in_tetra_pair_row);
            
            [minEDGEForm_objs,pos_in_obj_list_minEDGEForm] = findobjID(obj_list,randomized_tetra_pairsID(:,3));
            [minEDGEForm_objs_IDs] = extract_IDs(minEDGEForm_objs);
            tetra_pair_indices_inobj_list = 1:size(minEDGEForm_objs,1);
            find_indices_in_array = @(ID) tetra_pair_indices_inobj_list( minEDGEForm_objs_IDs==ID  );
            pos_in_tetra_pair_row = arrayfun(find_indices_in_array,randomized_tetra_pairsID(:,3));
            minEDGEForm_objs = minEDGEForm_objs(pos_in_tetra_pair_row);
            pos_in_obj_list_minEDGEForm = pos_in_obj_list_minEDGEForm(pos_in_tetra_pair_row);
            
            [maxEDGEForm_objs,pos_in_obj_list_maxEDGEForm] = findobjID(obj_list,randomized_tetra_pairsID(:,4));
            [maxEDGEForm_objs_IDs] = extract_IDs(maxEDGEForm_objs);
            tetra_pair_indices_inobj_list = 1:size(maxEDGEForm_objs,1);
            find_indices_in_array = @(ID) tetra_pair_indices_inobj_list( maxEDGEForm_objs_IDs==ID  );
            pos_in_tetra_pair_row = arrayfun(find_indices_in_array,randomized_tetra_pairsID(:,4));
            maxEDGEForm_objs = maxEDGEForm_objs(pos_in_tetra_pair_row);
            pos_in_obj_list_maxEDGEForm = pos_in_obj_list_maxEDGEForm(pos_in_tetra_pair_row);
            
            tetra_pair_obj_list = [minEDGELive_objs,maxEDGELive_objs,minEDGEForm_objs,maxEDGEForm_objs];
            tetra_pair_pos_in_obj_list = ...
                [pos_in_obj_list_minEDGElive',pos_in_obj_list_maxEDGElive',pos_in_obj_list_minEDGEForm',pos_in_obj_list_maxEDGEForm'];
            % Doing metropolis
            %             tic
            tetra_pair_pos_in_obj_list_toChangeMarker = zeros(size(tetra_pair_pos_in_obj_list));
            %             tetra_pair_obj_list_updated = NaN(size(randomized_tetra_pairsID));
            tetra_pair_obj_list_updated = [membrane_particles_list(membrane_particle(1),size(randomized_tetra_pairsID,1)),...
                membrane_particles_list(membrane_particle(1),size(randomized_tetra_pairsID,1)),...
                membrane_particles_list(membrane_particle(1),size(randomized_tetra_pairsID,1)),...
                membrane_particles_list(membrane_particle(1),size(randomized_tetra_pairsID,1))];
            updated_indicator_all = zeros(size(tetra_pair_obj_list,1),1);
            parfor tetra_pair_row = 1:size(tetra_pair_obj_list,1)
                current_tetra_pair = tetra_pair_obj_list(tetra_pair_row,:);
                [updated_current_tetra_pair,updated_indicator] = update_tetra_pair(current_tetra_pair,hardball_radius,overlapping_energy);
                [updated_obj_list,can_move,del_energy] = ...
                    check_metropolis_MCSweep_obj_list(updated_current_tetra_pair,current_tetra_pair,surface_tension,hardball_radius,pressure,30*pi/180,(360-30)*pi/180,temperature_current);
                
                if can_move && (updated_indicator==1) ==1
                    fprintf('   >>> EDGE %d_%d broken and EDGE %d_%d formed with del energy %d\n\n\n',...
                        current_tetra_pair(1).ID,current_tetra_pair(2).ID,updated_current_tetra_pair(3).ID,updated_current_tetra_pair(4).ID,del_energy);
                    tetra_pair_pos_in_obj_list_toChangeMarker(tetra_pair_row,:) = [1,1,1,1];
                    tetra_pair_obj_list_updated(tetra_pair_row,:) = updated_obj_list;
                    updated_indicator_all(tetra_pair_row) = 1;
                else
                    tetra_pair_obj_list_updated(tetra_pair_row,:) = current_tetra_pair;
                end
                
            end
            
            Edge_flip_rate = sum(updated_indicator_all)/length(updated_indicator_all)*100;
            fprintf('Edge flipping Rate = %d\n',Edge_flip_rate)
            updated_obj_list = obj_list;
            updated_obj_list( tetra_pair_pos_in_obj_list(:,1) )  = ...
                tetra_pair_obj_list_updated(:,1);
            updated_obj_list( tetra_pair_pos_in_obj_list(:,2) )  = ...
                tetra_pair_obj_list_updated(:,2);
            updated_obj_list( tetra_pair_pos_in_obj_list(:,3) )  = ...
                tetra_pair_obj_list_updated(:,3);
            updated_obj_list( tetra_pair_pos_in_obj_list(: ,4) )  = ...
                tetra_pair_obj_list_updated(: ,4);
        end
        
        function [obj_list,acceptance_rate] = ...
                diffuse_edge_flip(obj_list,hardball_radius,overlapping_energy,temperature_current,pressure,surface_tension)
            % Randomly select edges and their corresponding vertices----
            fprintf('Randomly selecting edges for edge flip\n');
            [randomized_tetra_pairsID] = randomize_tetraPair(obj_list);
            
            %             [all_IDs] = extract_IDs(obj_list);
            %             [edge_to_del_objs_list,new_edge_to_connect_objs_list] = randomized_vertex_tetrapairs_for_edgeflip(obj_list);
            
            fprintf('Making obj list tetra pair\n');
            [~,obj_list,Edge_flip_rate] = edge_flip_tetra_pair_obj_list(obj_list,randomized_tetra_pairsID,hardball_radius,overlapping_energy,temperature_current,pressure,surface_tension);
            acceptance_rate = Edge_flip_rate;
            [obj_list] = update_neighbours_normals(obj_list);
        end
        

        %% Nearest neighbourhood check %%%%%%%%%%%%%%%
        
        
        %% Prediction functions
        function [free_energy,curvature_energy] = calculate_energy_ofType(obj_list,type_list,pressure,surface_tension,temperature,hardball_rad,toplot,figure_handle)
            fprintf('EXTRACTING Energy info of different types \n')
            tic
            total_types = length(type_list);
            
            [type_all] = extract_type(obj_list);
            obj_ind = 1:size(obj_list,1);
            
            free_energy = cell(1,total_types);
            curvature_energy = cell(1,total_types);
            corresponding_mean_curvature = cell(1,total_types);
            corresponding_curvature = cell(1,total_types);
            
            type_rigidity = NaN(1,total_types);
            type_Ho = NaN(1,total_types);
            
            
            for type_ind = 1:total_types
                current_type = type_list(type_ind);
                obj_inds_current_type = obj_ind(ismember(type_all,current_type));
                objs_current_type = obj_list(obj_inds_current_type,:);
                
                
                type_rigidity(1,type_ind) = objs_current_type(1).bending_rigidity;
                type_Ho(1,type_ind) = objs_current_type(1).spontaneous_curvature;
                free_energy_current = NaN(length(obj_inds_current_type),1);
                curvature_energy_current = NaN(length(obj_inds_current_type),1);
                parfor obj_id = 1:size(objs_current_type,1)
                    obj = objs_current_type(obj_id,1);
                    [free_energy_current(obj_id),curvature_energy_current(obj_id),energy_surface_tension_normalized,energy_stored_by_pressure_outward] = ...
                        calculate_energy_from_geometry(obj,surface_tension,hardball_rad,pressure);
                end
                free_energy{1,type_ind} = free_energy_current;
                curvature_energy{1,type_ind} = curvature_energy_current;
                
                corresponding_mean_curvature{type_ind} = mean(extract_curvature(objs_current_type),2);
                corresponding_curvature{type_ind} = extract_curvature(objs_current_type);
            end
            
            toc
            % Plotting in figure
            curvature_median_plot = mean(type_Ho);
            curvature_min_plot = min(type_Ho)-3*(curvature_median_plot - min(type_Ho));
            curvature_max_plot = max(type_Ho)+3*(-curvature_median_plot + max(type_Ho));
            
            if toplot == 1
                figure(figure_handle);
                h2_sub = subplot(1,1,1);
                for type_ind = 1:total_types
                    type_bending_rigidity = type_rigidity(type_ind);
                    Ho = type_Ho(type_ind);
                    Hmin_all = corresponding_curvature{type_ind}(:,1);
                    Hmax_all = corresponding_curvature{type_ind}(:,2);
                    
                    %                     energy_profile_curve = @(Hmean) 0.5*type_bending_rigidity*(Hmean-Ho)^2;
                    %                     curvature_line_space = curvature_min_plot:(curvature_max_plot-curvature_min_plot)/100:curvature_max_plot;
                    %                     energy_curvature_profile = arrayfun(energy_profile_curve,curvature_line_space);
                    energy_profile_curve2D = @(Hmin,Hmax) 0.5*type_bending_rigidity*((Hmin-Ho)^2+(Hmax-Ho)^2);
                    [H_min_space,H_max_space] = meshgrid((curvature_min_plot:(curvature_max_plot-curvature_min_plot)/100:curvature_max_plot),...
                        (curvature_min_plot:(curvature_max_plot-curvature_min_plot)/200:curvature_max_plot));
                    energy_curvature_profile = zeros(size(H_min_space));
                    energy_curvature_profile(:) = arrayfun(energy_profile_curve2D,H_min_space,H_max_space);
                    if type_ind ~=1
                        hold(h2_sub, 'on');
                    end
                    %                     plot(h2_sub,curvature_line_space,energy_curvature_profile);hold(h2_sub, 'off');
                    %                     hold(h2_sub, 'on');plot(h2_sub,corresponding_mean_curvature{type_ind},curvature_energy{1,type_ind},'o');hold(h2_sub, 'off');
                    surf(h2_sub,H_min_space,H_max_space,energy_curvature_profile,'FaceAlpha',0.4);%view(h2_sub,2);
                    hold(h2_sub,'on'); scatter3(h2_sub,Hmin_all,Hmax_all,curvature_energy{1,type_ind},'filled');hold(h2_sub, 'off');daspect(h2_sub,[1,1,1]);
                end
            end
            hold(h2_sub, 'on'); surf(h2_sub,H_min_space,H_max_space,temperature*ones(size(H_min_space)),'FaceAlpha',0.2);
            zlim([0,2*temperature]);view(h2_sub,[-15.3599,4.5027]);
            hold(h2_sub, 'off');
            
        end
        
        function compare_energies(obj_to_plot,kick_displacement,pressure,surface_tension,temperature,hardball_radius,figure_number)
            current_fig = figure(figure_number);
            current_fig_sub1 = subplot(1,2,1);
            
            % plot the energy distribution of obj with varying curvature
            [current_obj_star_mesh] = extract_obj_star(obj_to_plot);
            Points_shifted_to_origin = current_obj_star_mesh.Points - ...
                (current_obj_star_mesh.Points(1,:)'*ones(1,size(current_obj_star_mesh.Points ,1)))';
            current_obj_star_mesh = triangulation(current_obj_star_mesh.ConnectivityList,Points_shifted_to_origin);
            trisurf(current_obj_star_mesh,'FaceColor','red','FaceAlpha',0.3); hold(current_fig_sub1,'off');
            hold(current_fig_sub1,'on');scatter3(current_fig_sub1,0,0,0,20,'Filled');hold(current_fig_sub1,'off');
            hold(current_fig_sub1,'on');quiver3(current_fig_sub1,0,0,0,...
                obj_to_plot.Normal_vect(1),obj_to_plot.Normal_vect(2),obj_to_plot.Normal_vect(3));hold(current_fig_sub1,'off');daspect(current_fig_sub1,[1,1,1]);
            
            for nei_id = 1:length(obj_to_plot.neighbours)
                hold(current_fig_sub1,'on');quiver3(current_fig_sub1,Points_shifted_to_origin(nei_id+1,1),Points_shifted_to_origin(nei_id+1,2),Points_shifted_to_origin(nei_id+1,3),...
                    obj_to_plot.neighbours_normal_old(nei_id,1),obj_to_plot.neighbours_normal_old(nei_id,2),obj_to_plot.neighbours_normal_old(nei_id,3));hold(current_fig_sub1,'off');daspect(current_fig_sub1,[1,1,1]);
            end
            vecnorm(obj_to_plot.Normal_vect(3))
            el_m = asin( -obj_to_plot.Normal_vect(3)/vecnorm(obj_to_plot.Normal_vect) )*180/pi + 90;
            az_m = acos( dot([0,1],-obj_to_plot.Normal_vect(1:2)./vecnorm(obj_to_plot.Normal_vect(1:2))   ) )*180/pi;
            %             azimuth_for_cam_placement = acos( obj_to_plot.Normal_vect(1)/vecnorm(obj_to_plot.Normal_vect(1:2)) )*180/pi
            view([az_m,el_m]);
            
            updated_obj_to_plot = obj_to_plot;
            updated_obj_to_plot.Pos = obj_to_plot.Pos + kick_displacement/1*obj_to_plot.Normal_vect;
            current_obj_star_mesh_points_positive = current_obj_star_mesh.Points;
            current_obj_star_mesh_points_positive(1,:) = current_obj_star_mesh_points_positive(1,:) + kick_displacement/1*updated_obj_to_plot.Normal_vect;
            current_obj_star_mesh_positive = triangulation(current_obj_star_mesh.ConnectivityList,current_obj_star_mesh_points_positive);
            updated_obj_to_plot = calculate_topological_characteristics(updated_obj_to_plot,hardball_radius,0);
            [free_energy_total_positive,energy_curvature_positive,energy_surface_tension_normalized_positive,energy_stored_by_pressure_positive] = ...
                calculate_energy_from_geometry(updated_obj_to_plot,surface_tension,hardball_radius,pressure);
            hold(current_fig_sub1,'on');trisurf(current_obj_star_mesh_positive,'FaceColor','red','FaceAlpha',0.1); hold(current_fig_sub1,'off');
            
            updated_obj_to_plot = obj_to_plot;
            updated_obj_to_plot.Pos = obj_to_plot.Pos - kick_displacement/1*obj_to_plot.Normal_vect;
            current_obj_star_mesh_points_negative = current_obj_star_mesh.Points;
            current_obj_star_mesh_points_negative(1,:) = current_obj_star_mesh_points_negative(1,:) - kick_displacement/1*updated_obj_to_plot.Normal_vect;
            current_obj_star_mesh_negative = triangulation(current_obj_star_mesh.ConnectivityList,current_obj_star_mesh_points_negative);
            updated_obj_to_plot = calculate_topological_characteristics(updated_obj_to_plot,0,0);
            [free_energy_total_negative,energy_curvature_negative,energy_surface_tension_normalized_negative,energy_stored_by_pressure_negative] = ...
                calculate_energy_from_geometry(updated_obj_to_plot,surface_tension,hardball_radius,pressure);
            hold(current_fig_sub1,'on');trisurf(current_obj_star_mesh_negative,'FaceColor','red','FaceAlpha',0.1); hold(current_fig_sub1,'off');
            
            % distribution about quantal steps in the cube (kick displament^3)
            deformations_along_normal_direction = -kick_displacement:2*kick_displacement/100:kick_displacement;
            curvature = NaN(size(deformations_along_normal_direction));
            free_energy_total = NaN(size(deformations_along_normal_direction));
            energy_curvature = NaN(size(deformations_along_normal_direction));
            energy_surface_tension = NaN(size(deformations_along_normal_direction));
            energy_stored_by_pressure = NaN(size(deformations_along_normal_direction));
            pos_displaced = NaN(length(deformations_along_normal_direction),3);
            pos_displaced_normal = NaN(length(deformations_along_normal_direction),1);
            parfor deform_ind = 1:length(deformations_along_normal_direction)
                deform_current = deformations_along_normal_direction(deform_ind);
                updated_obj_to_plot = obj_to_plot;
                updated_obj_to_plot.Pos = obj_to_plot.Pos + deform_current*obj_to_plot.Normal_vect;
                pos_displaced(deform_ind,:) =  updated_obj_to_plot.Pos-obj_to_plot.Pos;
                pos_displaced_normal(deform_ind) = deform_current;
                [updated_obj_to_plot.facePairs,anomaly] = update_facePair(updated_obj_to_plot,updated_obj_to_plot.faces);
                updated_obj_to_plot = calculate_topological_characteristics(updated_obj_to_plot,0,0);
                [free_energy_total(deform_ind),energy_curvature(deform_ind),energy_surface_tension(deform_ind),energy_stored_by_pressure(deform_ind)] = ...
                    calculate_energy_from_geometry(updated_obj_to_plot,surface_tension,hardball_radius,pressure);
                curvature(deform_ind) = mean(updated_obj_to_plot.DarbouxFrame{2,1}(1:2))
            end
            current_fig = figure(figure_number);
            current_fig_sub2 = subplot(1,2,2);
            plot(current_fig_sub2,deformations_along_normal_direction,energy_curvature,'blue'); hold(current_fig_sub2,'off');
            hold(current_fig_sub2,'on');plot(current_fig_sub2,deformations_along_normal_direction,energy_surface_tension,'r');hold(current_fig_sub2,'off');
            hold(current_fig_sub2,'on');plot(current_fig_sub2,deformations_along_normal_direction,energy_stored_by_pressure,'green');hold(current_fig_sub2,'off');
            hold(current_fig_sub2,'on');plot(current_fig_sub2,deformations_along_normal_direction,temperature*ones(size(deformations_along_normal_direction)),'yellow');hold(current_fig_sub2,'off');
            hold(current_fig_sub2,'on');plot(current_fig_sub2,deformations_along_normal_direction,curvature,'black');hold(current_fig_sub2,'off');
            hold(current_fig_sub1,'on');scatter3(current_fig_sub1,pos_displaced(:,1),pos_displaced(:,2),pos_displaced(:,3),10,pos_displaced_normal,'filled');hold(current_fig_sub1,'off');
            legend(current_fig_sub2,'energy_curvature','surface_tension','pressure' ,'temperature','curvature' );
            % Distribution of edge length
            
            % distribution of curved volume
            
            % distribution of surface tension energy
            
            % distribution of volume energy
        end
        
        
        %% Plotting, query, refining functions
        function [IDs] = extract_IDs(obj_list)
            extract_ID_func = @(id) obj_list(id).ID;
            IDs = arrayfun(extract_ID_func,1:size(obj_list,1));
        end
        
        function [indices] = indices_in_objlist(obj_list,IDs)
            [~,indices] = findobjID(obj_list,IDs);
        end
        
        function [coors] = extract_coors_points(obj_list)
            coors = NaN(size(obj_list,1),3);
            parfor obj_id = 1:size(obj_list,1)
                coors(obj_id,:) = obj_list(obj_id).Pos;
            end
        end
        
        function [normal_vects_all] = extract_Normal_vects(obj_list)
            extract_Normal_func1 = @(id) obj_list(id).Normal_vect(1);
            extract_Normal_func2 = @(id) obj_list(id).Normal_vect(2);
            extract_Normal_func3 = @(id) obj_list(id).Normal_vect(3);
            normal_vects_all = [(arrayfun(extract_Normal_func1,1:size(obj_list,1)))',...
                (arrayfun(extract_Normal_func2,1:size(obj_list,1)))',...
                (arrayfun(extract_Normal_func3,1:size(obj_list,1)))'];
        end
        
        function [coors] = extract_coors_points_ID_specific(obj_list)
            coors = NaN(size(obj_list,1),3);
            ID_extract = @(id)obj_list(id).ID;
            ID_extract_all = arrayfun(ID_extract,1:size(obj_list,1));
            for obj_id = 1:size(obj_list,1)
                ID_current = ID_extract_all(obj_id);
                coors(ID_current,:) = obj_list(obj_id).Pos;
            end
            
        end
        
        function [faces_all_mat_uniqued] = extract_all_faces(obj_list)
            faces_extract = @(id) obj_list(id).faces;
            num_faces_extract = @(id) size(obj_list(id).faces,1);
            num_faces_all = arrayfun(num_faces_extract,1:size(obj_list,1));
            total_faces = sum(num_faces_all);
            faces_all = cell(size(obj_list,1),1);
            
            for fc_id = 1:size(obj_list,1)
                faces_all{fc_id} = faces_extract(fc_id);
            end
            faces_all_mat = cell2mat(faces_all);
            faces_all_mat_uniqued = unique(sort(faces_all_mat,2),'rows');
        end
        
        function [edge_to_be_flipped] = randomized_Partner_from_Neighbour(obj_list)
            obj_ID = NaN(size(obj_list,1),1);
            extracted_partner = NaN(size(obj_list,1),1);
            %%extracted_opposite_vertices = NaN( size(obj_list,1),2 );
            parfor obj_ind = 1:size(obj_list,1)
                obj_ID(obj_ind) = obj_list(obj_ind).ID;
                %                 neighbours_idx = 1:length(obj_list(obj_ind).neighbours);
                neighbours_idx_randomized = randperm(length(obj_list(obj_ind).neighbours));
                neighbours_current = obj_list(obj_ind).neighbours(neighbours_idx_randomized) ;
                extracted_partner(obj_ind) =  neighbours_current(1);
                % finding faces containing the edge for opposite vertices
                %%current_faces = obj_list(obj_ind).faces;
                %%indication_mask = double(ismember(current_faces,[extracted_partner(obj_ind),obj_ID(obj_ind)]));
                %%interesting_faces = current_faces( sum(indication_mask,2)==2,: );
                %%interesting_vertices = interesting_faces( not( ismember(interesting_faces,[extracted_partner(obj_ind),obj_ID(obj_ind)]) ) );
                
            end
            %Removing  vertex flip events of same ID object from two
            %different edges
            edge_common_unique = unique(sort([extracted_partner,obj_ID],2),'rows');
            
            [edge_common_unique_based_oncol1,indicies] = unique(edge_common_unique(:,1));
            edge_common_unique_tmp = [edge_common_unique_based_oncol1,edge_common_unique(indicies,2)];
            edge_common_unique = edge_common_unique_tmp;
            [edge_common_unique_based_oncol2,indicies] = unique(edge_common_unique(:,2));
            edge_common_unique_tmp = [edge_common_unique(indicies,1),edge_common_unique_based_oncol2];
            edge_common_unique = edge_common_unique_tmp;
            
            
            edge_common = intersect(edge_common_unique(:,1),edge_common_unique(:,2));
            common_vertex_mask = ismember(edge_common_unique(:,1),edge_common) | ismember(edge_common_unique(:,2),edge_common);
            edge_common_unique(common_vertex_mask,:) = [];
            edge_to_be_flipped = edge_common_unique;
            rand_ind_seq = randperm(size(edge_to_be_flipped,1));
            edge_to_be_flipped = edge_to_be_flipped(rand_ind_seq,:);
        end
        
        function [updated_face_pairs,anomaly] = update_facePair(obj,updated_faces)
            faces_sorted = sort(updated_faces,2);
            all_edges = [faces_sorted(:,1),faces_sorted(:,2);...
                faces_sorted(:,1),faces_sorted(:,3);...
                faces_sorted(:,2),faces_sorted(:,3)];
            all_edges_unique = unique(all_edges,'rows');
            common_edges_mask = sum(double(ismember(all_edges_unique,obj.ID)),2)==1;
            common_edges = all_edges_unique(common_edges_mask,:);
            
            all_points_ID = [obj.ID,obj.neighbours'];
            all_points_coor = [obj.Pos;obj.neighbours_coor];
            
            % find face pairs for each edge.
            face_ids = 1:size(updated_faces,1);
            updated_face_pairs = NaN(size(common_edges,1),2);
            anomaly = false;
            for ed_pair = 1:size(common_edges,1)
                current_common_edge = common_edges(ed_pair,:);
                faces_containing_edge_mask = sum(double(ismember(updated_faces,current_common_edge)),2)==2;
                if length( face_ids(faces_containing_edge_mask) )~=2
                    fprintf('COMMON EDGE NOT FOUND \n')
                    anomaly = or(anomaly,true);
                else
                    updated_face_pairs(ed_pair,:) = face_ids(faces_containing_edge_mask);
                    % check order if anticlocked
                    face1_id = updated_face_pairs(ed_pair,1);
                    face2_id = updated_face_pairs(ed_pair,2);
                    face1 = updated_faces(face1_id,:);
                    face2 = updated_faces(face2_id,:);
                    faces_pair = [face1;face2];
                    faces_pair_local = mask_mat_to_sequence(faces_pair,all_points_ID);
                    [~,idx] = find(faces_pair_local(1,:)==1);
                    faces_pair_local(1,:) = circshift(faces_pair_local(1,:),1-idx);
                    [~,idx] = find(faces_pair_local(2,:)==1);
                    faces_pair_local(2,:) = circshift(faces_pair_local(2,:),1-idx);
                    
                    face1_centroid = mean( all_points_coor(faces_pair_local(1,:),:),1 );
                    face2_centroid = mean( all_points_coor(faces_pair_local(2,:),:),1 );
                    
                    direction_current = face2_centroid - face1_centroid;
                    direction_edgeNeighbor_face1 = ...
                        all_points_coor(faces_pair_local(1,3),:) - all_points_coor(faces_pair_local(1,2),:);
                    direction_edgeNeighbor_face2 = ...
                        all_points_coor(faces_pair_local(2,3),:) - all_points_coor(faces_pair_local(2,2),:);
                    
                    net_edge_direction = direction_edgeNeighbor_face1 + direction_edgeNeighbor_face2;
                    
                    if sign( dot(direction_current,net_edge_direction))<0
                        updated_face_pairs(ed_pair,:) = fliplr(updated_face_pairs(ed_pair,:));
                    end
                end
            end
        end
        
        function [neighbours_effected_obj_all] = extract_all_neighbours(obj_list)
            extract_neighbours_of_current_obj = @(id) (obj_list(id).neighbours)';
            neighbours_effected_obj_all = ...
                unique(cell2mat(arrayfun(extract_neighbours_of_current_obj,1,'UniformOutput',false)));
        end
        
        function [randomly_selected_neigbour] = random_neighbour(obj_list)
            num_objs = size(obj_list,1);
            random_neighbour_select = @(id) obj_list(id).neighbours(randi(end));
            randomly_selected_neigbour = arrayfun( random_neighbour_select,1:num_objs );
        end
        
        function [opposite_pairs,edge_selected] = opposite_vertices_of_random_edge(obj_list)
            opposite_pairs = NaN(size(obj_list,1),2);
            edge_selected = NaN(size(obj_list,1),2);
            %             obj_ID = NaN(size(obj_list,1),1);
            for ob_ind = 1:size(obj_list,1)
                obj_ID = obj_list(ob_ind).ID;
                corresponding_neighbour_vertex = obj_list(ob_ind).neighbours(randi(end));
                ob_faces = obj_list(ob_ind).faces;
                common_faces_ind_mask = sum(double(ismember(ob_faces,[corresponding_neighbour_vertex,obj_list(ob_ind).ID])),2)==2   ;
                common_faces = ob_faces(common_faces_ind_mask,:);
                opposite_vertices = common_faces(not(ismember(common_faces,[corresponding_neighbour_vertex,obj_list(ob_ind).ID])));
                opposite_vertices = unique(opposite_vertices(:));
                opposite_pairs(ob_ind,:) = opposite_vertices;
                edge_selected(ob_ind,:) = sort([obj_ID,corresponding_neighbour_vertex]);
            end
        end
        
        function [type_all] = extract_type(obj_list)
            extract_type = @(ind) obj_list(ind).type;
            type_all = arrayfun(extract_type,1:size(obj_list,1));
        end
        
        function [type_all_Ho] = extract_spontaneousCurvature(obj_list)
            extract_Ho = @(ind) obj_list(ind).spontaneous_curvature;
            type_all_Ho = arrayfun(extract_Ho,1:size(obj_list,1));
        end
        
        function [obj_list] = allign_faceOrder_based_on_faceNormal(obj_list)
            num_objs = size(obj_list,1);
            new_obj_list_faceordered = obj_list;
            parfor ob_ind = 1:num_objs
                current_obj = obj_list(ob_ind);
                normal_vect_current = current_obj.Normal_vect;
                neighbours_current = current_obj.neighbours;
                all_vertex_IDs = [current_obj.ID,neighbours_current'];
                all_vertex_index = 1:size(all_vertex_IDs,1);
                neighbours_current_coors = current_obj.neighbours_coor;
                all_vertex_coors = [current_obj.Pos;neighbours_current_coors];
                current_faces = current_obj.faces;
                current_obj_faces_maskedSequenced = mask_mat_to_sequence(current_faces,all_vertex_IDs);
                % check normal directionality
                for face_id = 1:size(current_faces,1)
                    current_face = current_obj_faces_maskedSequenced(face_id,:);
                    edge1 = all_vertex_coors(current_face(2),:) - all_vertex_coors(current_face(1),:);
                    edge2 = all_vertex_coors(current_face(3),:) - all_vertex_coors(current_face(2),:);
                    decision = sign( dot(cross(edge1,edge2),normal_vect_current) );
                    if decision<0
                        obj_list(ob_ind).faces(face_id,:) =  fliplr(obj_list(ob_ind).faces(face_id,:));
                        fprintf('   >>>object with ID %d face id %d corrected',ob_ind,face_id);
                    end
                end
            end
        end
        
        function [principle_curvatures] = extract_curvature(obj_list)
            darboux_eigen_val_extract = @(id) obj_list(id).DarbouxFrame{2,1}(1:2);
            principle_curvature_extracted = arrayfun(darboux_eigen_val_extract,1:size(obj_list,1),'UniformOutput',false);
            principle_curvatures = cell2mat(principle_curvature_extracted');
        end
        
        function [surface_mesh] = generate_mesh_from_objlist(obj_list)
            fprintf('Generating Membrane Mesh \n');
            tic
            [coors] = extract_coors_points_ID_specific(obj_list);
            [faces_all_mat_uniqued] = extract_all_faces(obj_list);
            %             [all_coors] = extract_coors_points(obj_list);
            [IDs] = extract_IDs(obj_list);
            surface_mesh = triangulation(faces_all_mat_uniqued,coors);
            all_face_normals = faceNormal(surface_mesh);
            %% ordering according to normal vect
            faces_all_mat_uniqued_ordered = faces_all_mat_uniqued;
            [normal_vects_all] = extract_Normal_vects(obj_list);
            parfor face_no = 1:size(faces_all_mat_uniqued,1)
                face_current = faces_all_mat_uniqued(face_no,:);
                current_normals = normal_vects_all(ismember(IDs,face_current),:);
                current_normals_sum = sum(current_normals,1);
                current_face_normal = all_face_normals(face_no,:);
                if sign(dot(current_face_normal,current_normals_sum))<0
                    faces_all_mat_uniqued_ordered(face_no,:) = fliplr(faces_all_mat_uniqued(face_no,:));
                end
            end
            surface_mesh = triangulation(faces_all_mat_uniqued_ordered,coors);
            toc
            fprintf('Membrane Mesh generated \n');
        end        
        
        function [face_curvature] = generate_curvature_plotSurface(obj_list,surface_mesh)
            %num_faces = size(surface_mesh.Connectivity,1);
            fprintf('Extracting mean curvature for faces for CURVATURE SURFACE PLOT\n')
            [All_IDs] = extract_IDs(obj_list);
            obj_indices = 1:size(obj_list,1);
            find_indices_in_array = @(ID) obj_indices( All_IDs == ID);
            %             find_indices_in_array = @(ID) bi_pair_indices_inobj_list( objs_for_vertex_flip_from_edgesCol1_IDs==ID  );
            col1_IDs = surface_mesh.ConnectivityList(:,1);
            col2_IDs = surface_mesh.ConnectivityList(:,2);
            col3_IDs = surface_mesh.ConnectivityList(:,3);
            pos_of_IDs_col1inObjList = arrayfun(find_indices_in_array,col1_IDs);
            pos_of_IDs_col2inObjList = arrayfun(find_indices_in_array,col2_IDs);
            pos_of_IDs_col3inObjList = arrayfun(find_indices_in_array,col3_IDs);
            
            
            extract_mean_principle_curv = @(ind) mean(obj_list(ind).DarbouxFrame{2,1}(1:2),2);
            principle_curv_mat = [arrayfun(extract_mean_principle_curv,pos_of_IDs_col1inObjList),...
                arrayfun(extract_mean_principle_curv,pos_of_IDs_col2inObjList),...
                arrayfun(extract_mean_principle_curv,pos_of_IDs_col3inObjList)];
            face_curvature = mean(principle_curv_mat,2);
            
        end
        
        function [star_mesh] = extract_obj_star(obj)
            obj_ID = obj.ID;
            obj_neighbours = obj.neighbours;
            all_IDs = [obj_ID,obj_neighbours'];
            
            obj_pos = obj.Pos;
            obj_neighbours_coor = obj.neighbours_coor;
            all_coors = [obj_pos;obj_neighbours_coor];
            
            obj_faces = obj.faces;
            obj_faces_local = mask_mat_to_sequence(obj_faces,all_IDs);
            star_mesh = triangulation(obj_faces_local,all_coors);
        end
        
        function [updated_neighbour_info_obj_list] = update_neighbour_info(obj_list_position_updated)
            [all_coors] = extract_coors_points(obj_list_position_updated);
            [all_IDS] =   extract_IDs(obj_list_position_updated);
            updated_neighbour_info_obj_list = obj_list_position_updated;
            parfor obj_no = 1:size(obj_list_position_updated,1)
                current_obj = obj_list_position_updated(obj_no);
                neighbour_IDS = current_obj.neighbours;
                [~,locb] = ismember(neighbour_IDS,all_IDS);
                current_obj.neighbours_coor = all_coors(locb,:);
                updated_neighbour_info_obj_list(obj_no) = current_obj;
                fprintf('   >>> Neighbour information updated to list %d\n',obj_no);
            end
        end
        
        function [voronoi_area] = extract_avg_voronoi_area_volume(obj)
            obj_faces = obj.faces;
            obj_neighbours = obj.neighbours;
            obj_ID = obj.ID;
            local_vertex_IDs = [obj_ID,obj_neighbours'];
            local_vertex_coors = [obj.Pos;obj.neighbours_coor];
            
            local_faces = mask_mat_to_sequence(obj_faces,local_vertex_IDs);
            %% the area normalization is missing in energy calculation ??????
            current_facePairs = obj.facePairs;
            
            
        end
        
        function [obj_list] = update_geometry_obj_list(obj_list,hardball_radius,overlapping_energy)
            [obj_list] = update_neighbour_info(obj_list);
            parfor ob_no = 1:size(obj_list,1)
                obj = obj_list(ob_no);
                updated_obj = calculate_topological_characteristics(obj,hardball_radius,overlapping_energy);
                obj_list(ob_no) = updated_obj;
                fprintf('   >>>(FUNCTION:update_geometry_obj_list) Object.ID %d topology updated \n',obj.ID);
            end
            [obj_list] = update_neighbours_normals(obj_list);
        end
        
        function [energy_curvature_sum,totalVolume,totalArea] = extract_energy_obj_list(obj_list,surface_tension,hardball_radius,pressure)
            energy_curvature_list = NaN(size(obj_list,1),1);
            parfor ob_no = 1:size(obj_list,1)
                obj = obj_list(ob_no);
                [~,energy_curvature,~,~] = ...
                    calculate_energy_from_geometry(obj,surface_tension,hardball_radius,pressure)
                energy_curvature_list(ob_no) = energy_curvature;
                fprintf('   >>>(FUNCTION:extract_energy_obj_list) Object.ID %d CURVATURE energy extracted \n',obj.ID);
            end
            energy_curvature_sum = sum(energy_curvature_list(:));
            fprintf('>>>(FUNCTION:extract_energy_obj_list) GENERATING MESH for global energy\n');
            [surface_mesh] = generate_mesh_from_objlist(obj_list);
            [totalVolume,totalArea] = area_volume_surface(surface_mesh.Points',surface_mesh.ConnectivityList');
            
        end        
        
        
    end
end