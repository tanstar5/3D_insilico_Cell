classdef membrane_patch
    
    properties
        %% Identity
        ID
        Pos
        Normal_vect
        
        %% Topology
        faces
        neighbours
        neighbours_coor
        neighbours_normal_old;
        neighbours_lipid_compositions_up
        neighbours_lipid_compositions_down
        
        %% u_Sigma_T ensemble properties
        
        chemical_potential_initial % standard chemical potential based on initial area and number of particles
        n_particles_up %total_number of particles in upper leaflet of bilayer
        n_particles_down %total_number of particles in upper leaflet of bilayer
        avg_area_per_lipid %avergage area occupied by a lipid = tot_area/num_lipids
        
        lipid_ratio_up   % 3 X 1 row vector representing molar ratio of +,-,0 curvature lipid in upper leaflet
        lipid_ratio_down   % 3 X 1 row vector representing molar ratio of +,-,0 curvature lipid in upper leaflet
        
        %% Energies and other derived geometrical quantities
        Av_vertex
        c1 %principle curvature1
        c2 %principle curvature2
        principal_directions
        coarse_grained_curvature_energy
        
        %% Signalling properties
        protein_bound_type_ID
    end
    
    
    methods
        %% Membrane patch list Initializing functions
        function obj = membrane_patch(ID)
            if nargin ~= 0
                obj.ID = ID;
            end
        end
        
        function [obj_list] = membrane_patch_list(first_obj,num_of_particles)
            obj_list(1,1) =  first_obj;
            for particle = 2:num_of_particles
                obj = membrane_patch(particle);
                obj_list(particle,1) =  obj;
            end
        end
        
        function [obj_list,total_surface_area] = load_spatial_properties_from_mesh(obj_list,membrane_mesh)
            fprintf('>>> Function: load_spatial_properties_from_mesh; CORRECTING ORIENTATION...\n\n');
            all_points = membrane_mesh.Points;
            triangles = membrane_mesh.ConnectivityList;
            face_normal_vects_all = faceNormal(membrane_mesh);
            %             normal_vects_at_each_points = NaN(size(all_points)) ;
            % Correcting the orientation of triangles on anticlockwise order seen from an inside point in a cell
            triangles_arranged = NaN(size(triangles));
            face_normal_vects_all_arranged = NaN(size(face_normal_vects_all));
            area_noter = NaN(size(triangles,1),1);
            for tri_id = 1:size(triangles,1)
                tri_current = triangles(tri_id,:);
                tri_points = all_points(tri_current,:);
                tri_centroid = mean(tri_points,1);
                area_noter(tri_id) = area_triangle(tri_points);
                if sign( dot(face_normal_vects_all(tri_id,:),tri_centroid-[0,0,0] ))<0
                    %                     triangles_arranged(tri_id,:) = fliplr(tri_current);
                    triangles_arranged(tri_id,:) = (tri_current);
                    face_normal_vects_all_arranged(tri_id,:) = -face_normal_vects_all(tri_id,:);
                    fprintf('Function: load_spatial_properties_from_mesh; %d triangle anticloked \n',tri_id);
                else
                    triangles_arranged(tri_id,:) = tri_current;
                    face_normal_vects_all_arranged(tri_id,:) = face_normal_vects_all(tri_id,:);
                end
            end
            total_surface_area = sum(area_noter);
            
            triangles = triangles_arranged;
            
            fprintf('>>> Function: load_spatial_properties_from_mesh; NEIGHBOUR COORS LOADING...\n\n');
            parfor particle = 1:size(obj_list,1)
                obj = obj_list(particle);
                obj.Pos = all_points(particle,:);
                triangle_mask = ismember(triangles,particle);
                triangle_mask_rows = triangle_mask(:,1)|triangle_mask(:,2)|triangle_mask(:,3);
                interesting_triangles = triangles(triangle_mask_rows,:);
                obj.neighbours = unique(interesting_triangles(not(ismember(interesting_triangles,particle))));
                obj.neighbours_coor = all_points(obj.neighbours,:);
                fprintf('>>> Function: load_spatial_properties_from_mesh; generating face info...\n');
                %                 face_normal_vect_current = face_normal_vects_all(triangle_mask_rows,:);
                faces_of_current_vertex = triangles(triangle_mask_rows,:);
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
                        vertex_of_interest = edge_current(2);
                    else
                        if particle == 3526
                            disp('')
                        end
                        vertex_of_interest = edge_current(edge_current~=vertex_of_interest);
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
                %                 normal_cross_vect = cross(obj.faceNormals(2,:),obj.faceNormals(1,:));
                %                 common_edge = intersect(obj.faces(1,:),obj.faces(2,:));
                %                 common_neigbourhood_vertex = common_edge( not(ismember(common_edge,commonVertex)) );
                %                 edge_outward_vect = all_points(common_neigbourhood_vertex,:) - membrane_mesh.Points(commonVertex,:);
                %                 decision = dot(edge_outward_vect, normal_cross_vect);
                %                 if decision<=0
                %                     obj.faces = flipud(obj.faces);
                % %                     obj.faceNormals = flipud(obj.faceNormals);
                % %                     obj.faceArea = flipud(obj.faceArea);
                %                 end
                obj.protein_bound_type_ID = 0;
                obj_list(particle) = obj;
                fprintf('>>> Function: load_spatial_properties_from_mesh; %d vertex neighbours loaded\n',particle);
            end
        end
        
        %% Derived Geometrical Quantities Functions
        function [updated_obj] = derive_geometrical_quatities(obj)
            fprintf('>>> Function: derive_geometrical_quatities; ID %d object geometry deriving\n',obj.ID);
            [c1_current,c2_current,r,Normal_vect_current,Av,~,principal_directions_current] = determine_principal_curvatures(obj);
            updated_obj = obj;
            updated_obj.c1 = c1_current; %2*principle_curvatures(1);
            updated_obj.c2 = c2_current; %2*principle_curvatures(2);
            updated_obj.principal_directions = principal_directions_current; 
            updated_obj.Normal_vect = Normal_vect_current;%Normal_vect_obj;
            updated_obj.Av_vertex = Av; 
            
        end
        
        function [updated_list] = derive_geometrical_quatities_all(obj_list)
            updated_list = obj_list;
            fprintf('>>> Function: derive_geometrical_quatities_all; deriving Geometrical quantities\n');
            parfor obj_no = 1:size(obj_list,1)
                [updated_obj] = derive_geometrical_quatities(obj_list(obj_no));
                updated_list(obj_no) = updated_obj;
            end
        end
        
        %% Initialize parameters and initial conditions functions
        function [obj_list] = determine_num_particles_per_patch_basedArea(obj_list,avg_molecules_per_triangle)
            fprintf('>>> Function: determine_num_particles_per_patch_basedArea; Calculating particles per patch in each leaflet \n');
            [surface_mesh] = generate_mesh_from_objlist(obj_list);
            num_triangles = size(surface_mesh.ConnectivityList,1);
            [totalVolume,totalArea] = area_volume_surface(surface_mesh.Points',surface_mesh.ConnectivityList');
            total_number_of_molecules = avg_molecules_per_triangle*num_triangles;
            density_of_molecules = total_number_of_molecules/totalArea;
            
            num_patches = size(obj_list,1);
            parfor patch_no = 1:num_patches
                obj = obj_list(patch_no);
                [star_mesh] = extract_obj_star(obj);
                area = 0;
                for tri_no = 1:size(star_mesh.ConnectivityList,1)
                    Points = star_mesh.Points(star_mesh.ConnectivityList(tri_no,:),:);
                    [area] = area + area_triangle(Points);
                end
                area_patch_current = area/3;
                obj.n_particles_up = round(area_patch_current*density_of_molecules);
                obj.n_particles_down = round(area_patch_current*density_of_molecules);
                obj.avg_area_per_lipid = 1/density_of_molecules;
                obj.Av_vertex = area_patch_current;
                fprintf('>>> Function: determine_num_particles_per_patch_basedArea; patch_no = %d; num_particle = %d; area = %d\n',obj.ID,obj.n_particles_up,area_patch_current);
                obj_list(patch_no) = obj;
            end
            
        end
        
        function [obj_list] = determine_num_particles_per_patch_basedAreaCorrected(obj_list,avg_molecules_per_triangle,total_surface_area,temperature)
            fprintf('>>> Function: determine_num_particles_per_patch_basedArea; Calculating particles per patch in each leaflet \n');
            [surface_mesh] = generate_mesh_from_objlist(obj_list);
            num_triangles = size(surface_mesh.ConnectivityList,1);
            [~,totalArea] = area_volume_surface(surface_mesh.Points',surface_mesh.ConnectivityList');
            total_number_of_molecules = avg_molecules_per_triangle*num_triangles;
            density_of_molecules = total_number_of_molecules/total_surface_area;
            
            num_patches = size(obj_list,1);
            for patch_no = 1:num_patches
                obj = obj_list(patch_no);
                [~,face_areas,projected_face_areas] = extract_obj_star(obj);


                total_projected_area = sum(projected_face_areas)/3;
                area_patch_current = total_projected_area;
                obj.n_particles_up = round(area_patch_current*density_of_molecules);
                obj.n_particles_down = round(area_patch_current*density_of_molecules);
                obj.avg_area_per_lipid = 1/density_of_molecules;
                obj.Av_vertex = sum(face_areas)/3;
                

                area_each_lipid = obj.avg_area_per_lipid;
                chemical_potential_current = temperature*log(obj.Av_vertex/area_each_lipid);
                obj.chemical_potential_initial = chemical_potential_current;
                
                fprintf('>>> Function: determine_num_particles_per_patch_basedArea; patch_no = %d; num_particle = %d; area = %d\n',obj.ID,obj.n_particles_up,area_patch_current);
                obj_list(patch_no) = obj;
            end
            
        end
        
        function [obj_list] = randomly_distribute_lipids(obj_list,avg_lipid_proportionality)
            % lipid proportinality determines the average density of lipids
            % [+ 0 -] type curvatue lipids
            avg_lipid_proportionality = avg_lipid_proportionality/sum(avg_lipid_proportionality);
            fprintf('>>> Function: randomly_distribute_lipids; Distributing lipids in each patch \n');
            for obj_no = 1:size(obj_list,1)
                obj = obj_list(obj_no);
                %upper leaflet calculation
                lipid_ratio_up_current = zeros(1,3);
                num_lipids_upper = obj.n_particles_up;
                
                %                 lipid_ratio_up_current = zeros(1,3);
                rand_sampling_pos = rand(num_lipids_upper,1);
                lipid_ratio_up_current(1) = length( rand_sampling_pos( rand_sampling_pos<=avg_lipid_proportionality(1)) );
                
                rand_sampling_pos = rand(num_lipids_upper,1);
                lipid_ratio_up_current(3) = length( rand_sampling_pos( rand_sampling_pos<=avg_lipid_proportionality(3)) );
                
                lipid_ratio_up_current(2) = num_lipids_upper - lipid_ratio_up_current(1) - lipid_ratio_up_current(3);
                obj.lipid_ratio_up = lipid_ratio_up_current;
                
                %lower leaflet calculation
                lipid_ratio_lower_current = zeros(1,3);
                num_lipids_upper = obj.n_particles_up;
                
                rand_sampling_pos = rand(num_lipids_upper,1);
                lipid_ratio_lower_current(1) = length( rand_sampling_pos( rand_sampling_pos<=avg_lipid_proportionality(1)) );
                
                rand_sampling_pos = rand(num_lipids_upper,1);
                lipid_ratio_lower_current(3) = length( rand_sampling_pos( rand_sampling_pos<=avg_lipid_proportionality(3)) );
                
                lipid_ratio_lower_current(2) = num_lipids_upper - lipid_ratio_lower_current(1) - lipid_ratio_lower_current(3);
                obj.lipid_ratio_down = lipid_ratio_lower_current;
                fprintf('>>> Function: randomly_distribute_lipids; patch_ID: %d with lipids upper [ %d %d %d ] and lower [ %d %d %d ]\n',...
                    obj.ID,obj.lipid_ratio_up(1),obj.lipid_ratio_up(2),obj.lipid_ratio_up(3),obj.lipid_ratio_down(1),obj.lipid_ratio_down(2),obj.lipid_ratio_down(3));
                obj_list(obj_no) = obj;
            end
            
            [lipid_comp_upper] = extract_lipids_composition_ID_specific(obj_list,1);
            parfor patch_no = 1:size(obj_list,1)
                current_patch = obj_list(patch_no);
%                 current_patch.Pos = all_points(current_patch.ID,:);
%                 current_patch.neighbours_coor = all_points(current_patch.neighbours,:);
                current_patch.neighbours_lipid_compositions_up = lipid_comp_upper(current_patch.neighbours,:);
                obj_list(patch_no) = current_patch;
                %                 fprintf('>>> Function: vertex_displacement_MC_local_move; Loading the updated pos and neighbours coor to patch ID = %d \n',current_patch.ID);
            end
            
        end
        
        %% Metropolis and statitical energy calculations of ensembles
        % local moves
        function [curvature_energy,entropy_of_mixing,surface_stretching_energy,distortion_energy,obj] =...
                coarse_grain_energy_mem_patch(obj,type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous,quantal_change,each_mole_content)
            % type properties is defined as [ bending_rigidity_row; spontaneous_curv; area_scale_factor ]
            
            [obj] = derive_geometrical_quatities(obj);
            
            avg_lipid_area = obj.avg_area_per_lipid;
            bending_rigidities_array = type_properties(1,1:3);
            spontaneous_curv_array = type_properties(2,1:3);
            area_scale_factor_array = type_properties(3,1:3);
            H_mean = (obj.c1 + obj.c2)/2;
            c1_curvature = obj.c1;
            c2_curvature = obj.c2;
            [distance_from_edge] = calculate_distance_from_edge(obj);
            std_radius = std(distance_from_edge);
            
            neighbourhood_lipid_comp = obj.neighbours_lipid_compositions_up;
            neighbours_area_from_comp = neighbourhood_lipid_comp*area_scale_factor_array';
            %             mean_neighbours_area = mean(neighbours_area_from_comp);
            %             std_lipid_comp_array = std( avg_lipid_area*(neighbourhood_lipid_comp*area_scale_factor_array') );
            
            %             std_vertex_fluctuation = std_lipid_comp_array/(mean_neighbours_area*pi)^(.5)/2;
            %             std_vertex_fluctuation = 1;
            
            area_scale_factor_array_RMS_fact = sqrt(sum( area_scale_factor_array.^2 ));
            area_scale_factor_array_normalized = area_scale_factor_array./area_scale_factor_array_RMS_fact;
            allowed_area_fluctuation = quantal_change*avg_lipid_area*( sqrt( 1 - 1/sqrt(3)*sum(area_scale_factor_array_normalized )   ) );
            %             allowed_area_fluctuation_rad = sqrt(allowed_area_fluctuation/pi);
            r_mean = sqrt( sum( avg_lipid_area*(obj.lipid_ratio_up.*area_scale_factor_array) )/pi );
            allowed_area_fluctuation_rad = allowed_area_fluctuation/(2*pi*r_mean);
            area_inequality_factor = std(avg_lipid_area*neighbours_area_from_comp);
            area_inequality_factor_rad = area_inequality_factor/(2*pi*r_mean);
            
            
            % Upperleaflet calculations
            fprintf('>>> Function: coarse_grain_energy_mem_patch; Calculating coarse grain energies of membrane patch ID = %d\n',obj.ID);
            [~,face_areas,projected_face_areas] = extract_obj_star(obj);
            upper_leaflet_lipid_composition_array = obj.lipid_ratio_up;
            
            total_area = sum(abs(face_areas))/3;
            projected_area_total = sum(abs(projected_face_areas))/3;
            N_particles = sum(obj.lipid_ratio_up);
            
            mole_fraction = upper_leaflet_lipid_composition_array./sum(upper_leaflet_lipid_composition_array);
            spontaneous_curv_array_mean = sum(mole_fraction.*spontaneous_curv_array);
            %             spontaneous_curv_array_mean = ((upper_leaflet_lipid_composition_array/sum(upper_leaflet_lipid_composition_array)).*spontaneous_curv_array);
            %             total_area = sum(avg_lipid_area.*upper_leaflet_lipid_composition_array);
            avg_bending_rigidity = mean(bending_rigidities_array);
%             if obj.ID == 100
%                 N_particles
%                 avg_lipid_area
%                 avg_bending_rigidity
%             end
            curvature_energy_type_array = ... % 1/2*number_of_ith_type*bending_rigidity_ith_type*(2*H_mean - spontaneous_curv_ith_type)^2
                0.5*N_particles*avg_bending_rigidity*( (-H_mean + 2*spontaneous_curv_array_mean).^2 )...
                + 0.5*gaussian_modulus*N_particles*(c1_curvature*c2_curvature);
            %                   0.5*N_particles*avg_lipid_area*avg_bending_rigidity*( (c1_curvature - spontaneous_curv_array_mean).^2 + (c2_curvature - spontaneous_curv_array_mean).^2 );
            %               0.5*N_particles*avg_lipid_area*avg_bending_rigidity*( (2*H_mean - spontaneous_curv_array_mean).^2 );
            %                 0.5*N_particles*avg_lipid_area*avg_bending_rigidity*( (2*H_mean - spontaneous_curv_array_mean).^2 );
            
            %             mole_fraction_array = upper_leaflet_lipid_composition_array./sum(upper_leaflet_lipid_composition_array);
            entropy_of_mixing_type_array = ... % no_of_ith_type*ln(mole_fraction_ith_type)
                -mole_fraction.*log(mole_fraction)*N_particles*temperature/each_mole_content*8.314;
            
            equillibrium_area_patch_array =  area_scale_factor_array.*upper_leaflet_lipid_composition_array*avg_lipid_area;
            
            
            check = abs(total_area-N_particles*avg_lipid_area);
            surface_stretching_energy = surface_modulus*N_particles*(total_area/sum(equillibrium_area_patch_array)-1).^2;
            
            %             surface_stretching_energy = surface_modulus*((total_area-N_particles*avg_lipid_area));
            %             curvature_energy = (curvature_energy_type_array) + 0.5*avg_lipid_area*gaussian_modulus*sum(upper_leaflet_lipid_composition_array)*(obj.c1*obj.c2);
%             curvature_energy = (curvature_energy_type_array);
            [curvature_energy] = calculate_enthalpic_energy_lipids(obj,type_properties(2,:),type_properties(1,:),H_mean);
            %             entropy_of_mixing = sum(entropy_of_mixing_type_array);
%             entropy_of_mixing =  sum( upper_leaflet_lipid_composition_array.*log(mole_fraction) )  ...
%                 -sum( upper_leaflet_lipid_composition_array.*log(upper_leaflet_lipid_composition_array) )  ...
%                 +sum(upper_leaflet_lipid_composition_array)*log( sum(upper_leaflet_lipid_composition_array) );
            entropy_of_mixing = 1*((N_particles*log(N_particles) - 2*sum(upper_leaflet_lipid_composition_array.*log(upper_leaflet_lipid_composition_array))))*temperature;
            
            %             distortion_energy = distortion_modulous*avg_lipid_area*sum(upper_leaflet_lipid_composition_array)*( std_vertex_fluctuation/std_radius -1 )^2;
            distortion_energy = real(distortion_modulous*N_particles*avg_lipid_area*( (allowed_area_fluctuation_rad + area_inequality_factor_rad)/std_radius -1 )^2);
            
        end
        
        function [curvature_energy,distortion_energy,obj] =...
                coarse_grain_energy_mem_patch_vertex_displacement(obj,obj_neighbours,type_properties,distortion_modulous,quantal_change)
            % type properties is defined as [ bending_rigidity_row; spontaneous_curv; area_scale_factor ]
            
            [obj] = derive_geometrical_quatities(obj);
            
            % Calculating curvature energy of the patch
            N_particles = sum(obj.lipid_ratio_up);
            bending_rigidities_array = type_properties(1,1:3);
            avg_bending_rigidity = mean(bending_rigidities_array);
            [principal_curvatures_neighbours,spontaneous_curvature_neigh,N_particles_neigh,~] = ...
                determine_principle_curvatures_at_neighbours_dis_at_center(obj,obj_neighbours,type_properties(2,:));
            
            [adjacent_areas_to_edges,total_area_patch] = calculate_adjacent_areas_to_neighbours(obj);
            
            weights_from_edge_for_curvature = 1/3*adjacent_areas_to_edges/total_area_patch;
            weights_for_center_for_curvature = obj.Av_vertex/total_area_patch;
            
            H_mean_edges = mean(principal_curvatures_neighbours,2);
            curvature_energy_contributionEdges = ...
                0.5*avg_bending_rigidity*(N_particles_neigh.*weights_from_edge_for_curvature.*( (-H_mean_edges + 2*spontaneous_curvature_neigh).^2 ));
            
            upper_leaflet_lipid_composition_array = obj.lipid_ratio_up;
            mole_fraction = upper_leaflet_lipid_composition_array./sum(upper_leaflet_lipid_composition_array);
            spontaneous_curv_array_mean = sum(mole_fraction.*type_properties(2,:));
            H_mean_center = (obj.c1 + obj.c2)/2;
            
            curvature_energy_contributionCenter = ...
                0.5*N_particles*avg_bending_rigidity*weights_for_center_for_curvature*(-H_mean_center + 2*spontaneous_curv_array_mean)^2  ;            
            
            curvature_energy = curvature_energy_contributionCenter + sum(curvature_energy_contributionEdges);       

            avg_lipid_area = obj.avg_area_per_lipid;
            area_scale_factor_array = type_properties(3,1:3);
            [distance_from_edge] = calculate_distance_from_edge(obj);
            std_radius = std(distance_from_edge);
            
            neighbourhood_lipid_comp = obj.neighbours_lipid_compositions_up;
            neighbours_area_from_comp = neighbourhood_lipid_comp*area_scale_factor_array';
            area_scale_factor_array_RMS_fact = sqrt(sum( area_scale_factor_array.^2 ));
            area_scale_factor_array_normalized = area_scale_factor_array./area_scale_factor_array_RMS_fact;
            allowed_area_fluctuation = quantal_change*avg_lipid_area*( sqrt( 1 - 1/sqrt(3)*sum(area_scale_factor_array_normalized )   ) );

            r_mean = sqrt( sum( avg_lipid_area*(obj.lipid_ratio_up.*area_scale_factor_array) )/pi );
            allowed_area_fluctuation_rad = allowed_area_fluctuation/(2*pi*r_mean);
            area_inequality_factor = std(avg_lipid_area*neighbours_area_from_comp);
            area_inequality_factor_rad = area_inequality_factor/(2*pi*r_mean);
            
            fprintf('>>> Function: coarse_grain_energy_mem_patch; Calculating coarse grain energies of membrane patch ID = %d\n',obj.ID);
            N_particles = sum(obj.lipid_ratio_up);       
            distortion_energy = real(distortion_modulous*N_particles*avg_lipid_area*( (allowed_area_fluctuation_rad + area_inequality_factor_rad)/std_radius -1 )^2);
            
        end
        
        
        
        function [decision] = metropolis_local_vertex_displacement(obj_new,obj_old,obj_neighbours,type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous,quantal_change,each_mole_content)
            fprintf('>>> Function: metropolis_local_vertex_displacement; Calculating Hamiltonians of %d old membrane patch \n',obj_old.ID);
            % Calculating energies of old_obj            
            [curvature_energy_old,distortion_energy_old,obj_old_updated] =...
                coarse_grain_energy_mem_patch_vertex_displacement(obj_old,obj_neighbours,type_properties,distortion_modulous,quantal_change);
            
            fprintf('>>> Function: metropolis_local_vertex_displacement; Calculating Hamiltonians of %d new membrane patch \n',obj_new.ID);
            % Calculating energies of new_obj        
            [curvature_energy_new,distortion_energy_new,updated_old_updated] =...
                coarse_grain_energy_mem_patch_vertex_displacement(obj_new,obj_neighbours,type_properties,distortion_modulous,quantal_change);
            
%             [del_entropy] = calculate_change_in_entropy_configurationBased(updated_old_updated,obj_old_updated,each_mole_content);
%             del_entropy = (entropy_of_mixing_new - entropy_of_mixing_old)/temperature;
            
            
            del_internal_energy_changeCurvature = curvature_energy_new - curvature_energy_old;
            
            [surface_tension_new,area_face_new,~] = calculate_changed_surface_tension(updated_old_updated,surface_modulus);
            [surface_tension_old,area_face_old,~] = calculate_changed_surface_tension(obj_old_updated,surface_modulus);            
            del_internal_surface_energy = surface_tension_old*( area_face_new - area_face_old );            
            
            del_internal_energy = del_internal_energy_changeCurvature + del_internal_surface_energy;
            
            decision = ...
                min( [1, exp( -del_internal_energy/temperature )] );
            
            
%             decision = ...
%                 min(1, exp( -(curvature_energy_new-curvature_energy_old)/temperature )*...
%                 exp( (del_entropy)/temperature )*...
%                 exp(-(surface_stretching_energy_new-surface_stretching_energy_old)/temperature )*...
%                 exp(-(distortion_energy_new-distortion_energy_old)/temperature ));
            
            fprintf('>>> Function: metropolis_local_vertex_displacement; decision = %d; curvature_factor = %d; surface_tension_energy = %d; distortion_energy = %d \n',...
                decision, exp( -del_internal_energy_changeCurvature/temperature ),...
                exp(-(del_internal_surface_energy)/temperature ),...
                exp(-(distortion_energy_new-distortion_energy_old)/temperature ))
            
        end
        
        function [decision] = metropolis_lipid_exchange_displacement(patch_pair_new,patch_pair_old,type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous,quantal_change,each_mole_content)
            fprintf('>>> Function: metropolis_lipid_exchange_displacement; Calculating Hamiltonians of [%d,%d] old membrane patch pair \n',patch_pair_old(1).ID,patch_pair_old(2).ID);
            % calculate energies of old pair
            obj_old_1 = patch_pair_old(1);
            obj_old_2 = patch_pair_old(2);
            [curvature_energy_old_1,entropy_of_mixing_old_1,surface_stretching_energy_old_1,distortion_energy_old_1,updated_obj_old_1] =...
                coarse_grain_energy_mem_patch(obj_old_1,type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous,quantal_change,each_mole_content);
            [curvature_energy_old_2,entropy_of_mixing_old_2,surface_stretching_energy_old_2,distortion_energy_old_2,updated_obj_old_2] =...
                coarse_grain_energy_mem_patch(obj_old_2,type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous,quantal_change,each_mole_content);
            updated_patch_pair_old = [updated_obj_old_1,updated_obj_old_2];
            curvature_energy_old_total = curvature_energy_old_1 + curvature_energy_old_2;
            entropy_of_mixing_old_total = entropy_of_mixing_old_1 + entropy_of_mixing_old_2;
            surface_stretching_energy_old_total = surface_stretching_energy_old_1 + surface_stretching_energy_old_2;
            distortion_energy_old_total = distortion_energy_old_1 + distortion_energy_old_2;
            
            fprintf('>>> Function: metropolis_lipid_exchange_displacement; Calculating Hamiltonians of [%d,%d] new membrane patch pair \n',patch_pair_new(1).ID,patch_pair_new(2).ID);
            % calculate energies of updated pair
            obj_new_1 = patch_pair_new(1);
            obj_new_2 = patch_pair_new(2);
            [curvature_energy_new_1,entropy_of_mixing_new_1,surface_stretching_energy_new_1,distortion_energy_new_1,updated_obj_new_1] =...
                coarse_grain_energy_mem_patch(obj_new_1,type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous,quantal_change,each_mole_content);
            [curvature_energy_new_2,entropy_of_mixing_new_2,surface_stretching_energy_new_2,distortion_energy_new_2,updated_obj_new_2] =...
                coarse_grain_energy_mem_patch(obj_new_2,type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous,quantal_change,each_mole_content);
            updated_patch_pair_new = [updated_obj_new_1,updated_obj_new_2];
            curvature_energy_new_total = curvature_energy_new_1 + curvature_energy_new_2;
            entropy_of_mixing_new_total = entropy_of_mixing_new_1 + entropy_of_mixing_new_2;
            surface_stretching_energy_new_total = surface_stretching_energy_new_1 + surface_stretching_energy_new_2;
            distortion_energy_new_total = distortion_energy_new_1 + distortion_energy_new_2;
            %             [del_entropy] = calculate_change_in_entropy(updated_patch_pair_new,updated_patch_pair_old,each_mole_content,temperature);
            del_entropy = entropy_of_mixing_new_total - entropy_of_mixing_old_total;
            %             [del_entropy] = calculate_change_in_entropy(patch_pair_new,patch_pair_old,each_mole_content);
            % metropolis

            % Calculating internal energy change from curvature and surfACE
            % tension
            [surface_tension_old1,area_face_old1,zero_surface_tension_area_old1] = calculate_changed_surface_tension(updated_obj_old_1,surface_modulus);
            [surface_tension_old2,area_face_old2,zero_surface_tension_area_old2] = calculate_changed_surface_tension(updated_obj_old_2,surface_modulus);
            
            [surface_tension_new1,area_face_new1,zero_surface_tension_area_new1] = calculate_changed_surface_tension(updated_obj_new_1,surface_modulus);
            [surface_tension_new2,area_face_new2,zero_surface_tension_area_new2] = calculate_changed_surface_tension(updated_obj_new_2,surface_modulus);
            
            del_internal_surface_energy1 = surface_tension_new1*( area_face_new1 - area_face_old1 );
            del_internal_surface_energy2 = surface_tension_new2*( area_face_new2 - area_face_old2 );
            
            del_curvature_energy = (curvature_energy_new_1 + curvature_energy_new_2)-(curvature_energy_old_1 + curvature_energy_old_2);
            del_surface_tension_energy = del_internal_surface_energy1 + del_internal_surface_energy2;
            
            del_internal_energy = del_curvature_energy + del_surface_tension_energy;
            
            % doing the entropy and chemical potentials calculations
            num_particle_changed1 = sum(updated_obj_new_1.lipid_ratio_up) - sum(updated_obj_old_1.lipid_ratio_up);
            num_particle_changed2 = sum(updated_obj_new_2.lipid_ratio_up) - sum(updated_obj_old_2.lipid_ratio_up);
            
            [~,chemical_potential_ob1_mole_term] = calculate_chemical_potential(updated_obj_new_1,temperature);
            [~,chemical_potential_ob2_mole_term] = calculate_chemical_potential(updated_obj_new_2,temperature);
            
            del_entropy_chemical_pot_ob1 = - (updated_obj_new_1.chemical_potential_initial + chemical_potential_ob1_mole_term)*num_particle_changed1/temperature;
            del_entropy_chemical_pot_ob2 = - (updated_obj_new_2.chemical_potential_initial + chemical_potential_ob2_mole_term)*num_particle_changed2/temperature;
            
            total_entropy_change_due_to_chemical_pot = del_entropy_chemical_pot_ob1 + del_entropy_chemical_pot_ob2;
            
            
            lipid_comp_ob1_old = updated_obj_old_1.lipid_ratio_up;
            lipid_comp_ob2_old = updated_obj_old_2.lipid_ratio_up;
            sum_type_lipids = lipid_comp_ob1_old + lipid_comp_ob2_old;
            lipid_mat = [lipid_comp_ob1_old;lipid_comp_ob2_old];
            entropy_config_old = sum(sum_type_lipids.*log(sum_type_lipids),2) - sum( lipid_mat(:).*log(lipid_mat(:)) );
            
            lipid_comp_ob1_new = updated_obj_new_1.lipid_ratio_up;
            lipid_comp_ob2_new = updated_obj_new_2.lipid_ratio_up;
            sum_type_lipids = lipid_comp_ob1_new + lipid_comp_ob2_new;
            lipid_mat = [lipid_comp_ob1_new;lipid_comp_ob2_new];
            entropy_config_new = sum(sum_type_lipids.*log(sum_type_lipids),2) - sum( lipid_mat(:).*log(lipid_mat(:)) );
            
            entropy_config_change = entropy_config_new - entropy_config_old;
%             entropy_config_change = 0
            
            free_energy_total = del_internal_energy - temperature*(entropy_config_change + 1*total_entropy_change_due_to_chemical_pot );
             
            decision = ...
                min(1, exp( -(free_energy_total)/temperature ));
                
%             decision = ...
%                 min(1, exp( -(curvature_energy_new_total-curvature_energy_old_total)/temperature )*...
%                 exp( (del_entropy)/temperature )*...
%                 exp(-(surface_stretching_energy_new_total-surface_stretching_energy_old_total)/temperature )*...
%                 exp(-(distortion_energy_new_total-distortion_energy_old_total)/temperature ));
            if obj_new_1.ID==1 || obj_new_2.ID==1
                disp('TRACK\n')
            end
            
            fprintf('>>> Function: metropolis_lipid_exchange_displacement; decision = %d; curvature_factor = %d; mixing_entropy = %d; stretching_factor = %d; distortion_energy = %d \n',...
                decision, exp( -(curvature_energy_new_total-curvature_energy_old_total)/temperature ),exp( (entropy_of_mixing_new_total - entropy_of_mixing_old_total)/temperature ),...
                exp(-(surface_stretching_energy_new_total-surface_stretching_energy_old_total)/temperature ),...
                exp(-(distortion_energy_new_total-distortion_energy_old_total)/temperature ));
            
        end
        
        function [ decision ] = metropolis_protein_infinite_reservoir_move(obj_list)
            
        end
        %% Monte Carlo moves
        function [obj_list] = vertex_displacement_MC_local_move(obj_list,max_displacement,type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous,quantal_change,each_mole_content)
            fprintf('>>> Function: vertex_displacement_MC_local_move; Performing local random moves... \n');
            % Randperm the index of choosing the patch order
            fprintf('>>> Function: vertex_displacement_MC_local_move; Randomizing membrane patch order ... \n');
            all_ID_extract_func = @(index) obj_list(index).ID;
            all_index = 1:size(obj_list,1);
            all_IDs = arrayfun(all_ID_extract_func,all_index);
            randomized_index = randperm(size(obj_list,1));
            displacement_IDs = all_IDs(randomized_index);
            [displacement_IDs_objs,ind_in_list_displacement] = findobjID(obj_list,displacement_IDs);
            all_points = extract_coors_points_ID_specific(obj_list);
            [lipid_comp_upper] = extract_lipids_composition_ID_specific(obj_list,1);
            
            rate_of_acceptance_array = zeros(size(obj_list,1),1);
            % Inside Loop
            for patch_no = 1:size(displacement_IDs_objs,1)
                fprintf('>>> Function: vertex_displacement_MC_local_move; Starting randomizing patch = %d\n',patch_no);
                current_patch = displacement_IDs_objs(patch_no);
                current_patch_neighbours = obj_list(current_patch.neighbours);
                %   update neighbourhood coors from running pos array
                if vecnorm(current_patch.Pos-all_points(current_patch.ID,:)) > 0
                    disp('BaaaM');
                end
                current_patch.Pos = all_points(current_patch.ID,:);
                current_patch.neighbours_coor = all_points(current_patch.neighbours,:);
                current_patch.neighbours_lipid_compositions_up = lipid_comp_upper(current_patch.neighbours,:);
                %   Change position of patch randomly along the orthogonal
                [current_patch] = derive_geometrical_quatities(current_patch);
                normal_vect_current = current_patch.Normal_vect;
                updated_patch = current_patch;
                % jumping distribution choosen Normal
                kick = randn(1,1);
                updated_patch.Pos = updated_patch.Pos + kick*max_displacement*normal_vect_current;
                %   check metropolis
                [decision] = metropolis_local_vertex_displacement(updated_patch,current_patch,current_patch_neighbours,type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous,quantal_change,each_mole_content);
                %   if accepted: update the running pos array
                
                % Calling uniform random number
                u = rand(1);
                if u <= decision
                    fprintf('>>> Function: vertex_displacement_MC_local_move; Move Accepeted decision = %d, u = %d\n\n',decision,u);
                    all_points(updated_patch.ID,:) = updated_patch.Pos;
                    rate_of_acceptance_array(updated_patch.ID) = 1;
                else
                    fprintf('>>> Function: vertex_displacement_MC_local_move; Move Not Accepeted decision = %d, u = %d\n\n',decision,u);
                end
            end
            
            % Update the pos of all patches and neighbours coor from running pos array
            fprintf('>>> Function: vertex_displacement_MC_local_move; Loading the updated pos and neighbours coor to list... \n');
            parfor patch_no = 1:size(obj_list,1)
                current_patch = obj_list(patch_no);
                current_patch.Pos = all_points(current_patch.ID,:);
                current_patch.neighbours_coor = all_points(current_patch.neighbours,:);
                current_patch.neighbours_lipid_compositions_up = lipid_comp_upper(current_patch.neighbours,:);
                obj_list(patch_no) = current_patch;
                %                 fprintf('>>> Function: vertex_displacement_MC_local_move; Loading the updated pos and neighbours coor to patch ID = %d \n',current_patch.ID);
            end
            % update geometrical properties of the whole patch list
            [obj_list] = derive_geometrical_quatities_all(obj_list);
            
            fprintf('>>> Function: vertex_displacement_MC_local_move; Completed; Rate_of_change = %d \n',sum(rate_of_acceptance_array)/length(rate_of_acceptance_array)*100);
        end
        
        function [obj_list] = lipid_exchange_MC_local_move(obj_list,quantal_change,diffusion_rates_ratio,type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous,each_mole_content)
            fprintf('>>> Function: lipid_exchange_MC_local_move; Performing lipid distribution ... \n');
            % Selecting random edges for lipid distribution
            fprintf('>>> Function: lipid_exchange_MC_local_move; Selecting random edges for lipid distribution... \n');
            [surface_mesh] = generate_mesh_from_objlist(obj_list);
            all_triangles = surface_mesh.ConnectivityList;
            all_triangles_sorted = sort(all_triangles,2);
            all_edges_ordered = [all_triangles_sorted(:,1),all_triangles_sorted(:,2);...
                all_triangles_sorted(:,2),all_triangles_sorted(:,3);...
                all_triangles_sorted(:,1),all_triangles_sorted(:,3)];
            num_edges = size(all_edges_ordered,1);
            randomized_edge_order = randperm(round(num_edges/3));
            first_col_patch_idx = all_edges_ordered(randomized_edge_order,1);
            [first_col_patch_objs,~] = findobjID(obj_list,first_col_patch_idx);
            second_col_patch_idx = all_edges_ordered(randomized_edge_order,2);
            [second_col_patch_objs,~] = findobjID(obj_list,second_col_patch_idx);
            obj_pair_list = [first_col_patch_objs,second_col_patch_objs];
            fprintf('>>> Function: lipid_exchange_MC_local_move; %d random edges were selected to distibute lipids \n', round(num_edges/2));
            
            [lipid_comp_upper_current] = extract_lipids_composition_ID_specific(obj_list,1);
            rate_of_acceptance_array = zeros(size(obj_pair_list,1),1);
            % Loop:
            for patch_pair = 1:size(obj_pair_list,1)
                fprintf('>>> Function: lipid_exchange_MC_local_move; Started Diffusion loop no: = %d\n',patch_pair);
                %   update lipid composition in pairs from lipid_comp_list and
                %   also neighbours lipid comp
                current_patch_pair = obj_pair_list(patch_pair,:);
                choose_order = randperm(2);
                patch_1 = current_patch_pair(choose_order(1));
                patch_2 = current_patch_pair(choose_order(2));
                
                patch_1.lipid_ratio_up = lipid_comp_upper_current(patch_1.ID,:);
                patch_1.neighbours_lipid_compositions_up = lipid_comp_upper_current(patch_1.neighbours,:);
                
                patch_2.lipid_ratio_up = lipid_comp_upper_current(patch_2.ID,:);
                patch_2.neighbours_lipid_compositions_up = lipid_comp_upper_current(patch_2.neighbours,:);
                
                %   Distribute lipids among all registered pairs
%                 [new_lipid_composition_patch1_up,new_lipid_composition_patch2_up,exchanged_noter] = ...
%                     distribute_lipids(patch_1,patch_2,quantal_change,2); % mode 1 is choosen to simulate jump distribution
                
                [new_lipid_composition_patch1_up,new_lipid_composition_patch2_up,exchanged_noter] = ...
                distribute_lipids_new2(patch_1,patch_2,quantal_change,diffusion_rates_ratio);
                
                updated_patch_1 = patch_1;
                updated_patch_2 = patch_2;
                
                updated_patch_1.lipid_ratio_up = new_lipid_composition_patch1_up;
                updated_patch_1.neighbours_lipid_compositions_up(updated_patch_1.neighbours==updated_patch_2.ID,:) = new_lipid_composition_patch2_up;
                
                updated_patch_2.lipid_ratio_up = new_lipid_composition_patch2_up;
                updated_patch_2.neighbours_lipid_compositions_up(updated_patch_2.neighbours==updated_patch_1.ID,:) = new_lipid_composition_patch1_up;
                
                
                updated_patch_pair = [updated_patch_1,updated_patch_2];
                current_patch_pair_ordered = [patch_1,patch_2];
                %  Check metropolis, and accordingly update lipid_list
                [decision] = metropolis_lipid_exchange_displacement(updated_patch_pair,current_patch_pair_ordered,type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous,quantal_change,each_mole_content);
                
                u = rand(1);
                if u <= decision*exchanged_noter
                    lipid_comp_upper_current(patch_1.ID,:) = new_lipid_composition_patch1_up;
                    lipid_comp_upper_current(patch_2.ID,:) = new_lipid_composition_patch2_up;
                    fprintf('>>> Function: lipid_exchange_MC_local_move; Move Accepeted decision = %d, exchanged_noter = %d, u = %d\n\n',decision,exchanged_noter,u);
                    rate_of_acceptance_array(patch_pair) = 1;
                else
                    lipid_comp_upper_current(patch_1.ID,:) = patch_1.lipid_ratio_up;
                    lipid_comp_upper_current(patch_2.ID,:) = patch_2.lipid_ratio_up;
                    fprintf('>>> Function: lipid_exchange_MC_local_move; Move Not Accepeted decision = %d, exchanged_noter = %d, u = %d\n\n',decision,exchanged_noter,u);
                end
                
                
            end
            
            % Update changed lipid comp to each patches from lipid_comp_list
            fprintf('>>> Function: lipid_exchange_MC_local_move; Loading the updated lipid list to all patches... \n');
            parfor patch_no = 1:size(obj_list,1)
                current_patch = obj_list(patch_no);
                current_patch.lipid_ratio_up = lipid_comp_upper_current(current_patch.ID,:);
                %                 current_patch.Pos = all_points(current_patch.ID,:);
                %                 current_patch.neighbours_coor = all_points(current_patch.neighbours,:);
                current_patch.neighbours_lipid_compositions_up = lipid_comp_upper_current(current_patch.neighbours,:);
                obj_list(patch_no) = current_patch;
                %                 fprintf('>>> Function: lipid_exchange_MC_local_move; Loading the updated lipid composistions to patch ID = %d \n',current_patch.ID);
            end
            
            % update geometrical properties of the whole patch list
            [obj_list] = derive_geometrical_quatities_all(obj_list);
            fprintf('>>> Function: lipid_exchange_MC_local_move; Completed; Rate_of_change = %d \n',sum(rate_of_acceptance_array)/length(rate_of_acceptance_array)*100);
            
        end
        
        function [] = protein_binding_move(obj_list,protein_obj_all)
            %%look into curvature at each membrane point from obj_list
            [obj_list] = derive_geometrical_quatities_all(obj_list);            
            
            %%protein types 
            Num_types_protein = protein_obj_all.types;
            
            %%Know the histry of bound proteins at mem site.
            all_mem_points_protein_info = @(ind) obj_list(ind).protein_bound_type_ID;
            
            
        end
        %% Accesorry functions
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
                    fprintf('>>> Function: update_facePair; ERROR--COMMON EDGE Not Found\n')
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
        
        function [star_mesh,face_areas,projected_face_areas] = extract_obj_star(obj)
            obj_ID = obj.ID;
            obj_neighbours = obj.neighbours;
            all_IDs = [obj_ID,obj_neighbours'];
            
            obj_pos = obj.Pos;
            obj_neighbours_coor = obj.neighbours_coor;
            all_coors = [obj_pos;obj_neighbours_coor];
            
            obj_faces = obj.faces;
            obj_faces_local = mask_mat_to_sequence(obj_faces,all_IDs);
            star_mesh = triangulation(obj_faces_local,all_coors);
            Normal_vect_current = obj.Normal_vect;
            face_normals = faceNormal(star_mesh);
            face_areas = NaN(size(obj_faces_local,1),1);
            for tri_id = 1:size(obj_faces_local,1)
                current_tri = obj_faces_local(tri_id,:);
                current_points = all_coors(current_tri,:);
                face_areas(tri_id,1) = area_triangle(current_points);
            end
            
            [projected_face_areas] = calculate_projected_area_of_faces(Normal_vect_current,face_normals,face_areas);
        end
        
        function [star_mesh,face_areas] = extract_obj_only_star(obj)
            obj_ID = obj.ID;
            obj_neighbours = obj.neighbours;
            all_IDs = [obj_ID,obj_neighbours'];
            
            obj_pos = obj.Pos;
            obj_neighbours_coor = obj.neighbours_coor;
            all_coors = [obj_pos;obj_neighbours_coor];
            
            obj_faces = obj.faces;
            obj_faces_local = mask_mat_to_sequence(obj_faces,all_IDs);
            star_mesh = triangulation(obj_faces_local,all_coors);
            %             Normal_vect_current = obj.Normal_vect;
            %             face_normals = faceNormal(star_mesh);
            face_areas = NaN(size(obj_faces_local,1),1);
            for tri_id = 1:size(obj_faces_local,1)
                current_tri = obj_faces_local(tri_id,:);
                current_points = all_coors(current_tri,:);
                face_areas(tri_id,1) = area_triangle(current_points);
            end
            
            %             [projected_face_areas] = calculate_projected_area_of_faces(Normal_vect_current,face_normals,face_areas);
        end
        
        function [distances] = calculate_distance_from_edge(obj)
            [star_mesh] = extract_obj_star(obj);
            num_triangles = size(star_mesh.ConnectivityList,1);
            all_points = star_mesh.Points;
            all_triangles = star_mesh.ConnectivityList;
            distances = NaN(num_triangles,1);
            %             fprintf('>>> Function: calculate_distance_from_edge; Calculating the distance of center from the edges \n')
            for tri_id = 1:num_triangles
                current_points = all_points( all_triangles(tri_id,:),: );
                area_current = area_triangle(current_points);
                
                tri_indices = all_triangles(tri_id,:);
                current_edge_points_mask = not(ismember( tri_indices,1 ));
                edge_points_id =  tri_indices( current_edge_points_mask );
                
                edge_points = all_points(edge_points_id,:);
                edge_length = vecnorm(edge_points(1,:)-edge_points(2,:));
                distances(tri_id) = 2*area_current/edge_length;
            end
            
        end
        
        function [surface_mesh] = generate_mesh_from_objlist(obj_list)
            fprintf('>>> Function: generate_mesh_from_objlist; Generating Surface mesh from obj_list\n');
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
            fprintf('>>> Function: generate_mesh_from_objlist; Membrane Mesh generated \n');
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
        
        function [neighbours_num_all] = extract_num_neighbours_ID_specific(obj_list)
            neighbours_num_all = NaN(size(obj_list,1),1);
            ID_extract = @(id)obj_list(id).ID;
            ID_extract_all = arrayfun(ID_extract,1:size(obj_list,1));
            for obj_id = 1:size(obj_list,1)
                ID_current = ID_extract_all(obj_id);
                neighbours_num_all(ID_current,:) = length(obj_list(obj_id).neighbours);
            end
            
        end
        
        function [lipid_comp] = extract_lipids_composition_ID_specific(obj_list,leaflet)
            lipid_comp = NaN(size(obj_list,1),3);
            ID_extract = @(id)obj_list(id).ID;
            ID_extract_all = arrayfun(ID_extract,1:size(obj_list,1));
            for obj_id = 1:size(obj_list,1)
                ID_current = ID_extract_all(obj_id);
                if leaflet == 1
                    lipid_comp(ID_current,:) = obj_list(obj_id).lipid_ratio_up;
                else
                    lipid_comp(ID_current,:) = obj_list(obj_id).lipid_ratio_down;
                end
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
        
        function [IDs] = extract_IDs(obj_list)
            extract_ID_func = @(id) obj_list(id).ID;
            IDs = arrayfun(extract_ID_func,1:size(obj_list,1));
        end
        
        function [normal_vects_all] = extract_Normal_vects(obj_list)
            extract_Normal_func1 = @(id) obj_list(id).Normal_vect(1);
            extract_Normal_func2 = @(id) obj_list(id).Normal_vect(2);
            extract_Normal_func3 = @(id) obj_list(id).Normal_vect(3);
            normal_vects_all = [(arrayfun(extract_Normal_func1,1:size(obj_list,1)))',...
                (arrayfun(extract_Normal_func2,1:size(obj_list,1)))',...
                (arrayfun(extract_Normal_func3,1:size(obj_list,1)))'];
        end
        
        function [obj,ind_in_list_arranged,corresponding_IDs_arranged] = findobjID(obj_list,ID)
            all_ID_extract = @(ind) obj_list(ind).ID;
            all_ID = arrayfun(all_ID_extract,[1:size(obj_list,1)]);
            all_index = (1:size(obj_list,1));
            if length(ID)<=1
                if length(all_ID(all_ID==ID))>1
                    fprintf('>>> Function: findobjID; ERROR multiple objects found \n');
                    obj = NaN;
                    ind_in_list = NaN;
                elseif length(all_ID(all_ID==ID))<1
                    fprintf('>>> Function: findobjID; ERROR no objects found \n');
                    obj = NaN;
                    ind_in_list = NaN;
                else
                    obj = obj_list(all_ID==ID);
                    ind_in_list = all_index(all_ID==ID);
                end
            else
                fprintf('>>> Function: findobjID; objects found \n');
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
        
        function [new_lipid_composition_patch1_up,new_lipid_composition_patch2_up,exchanged_noter] = ...
                distribute_lipids(patch_1,patch_2,quantal_change,type_of_move)
            
            fprintf('>>> Function: distribute_lipids; Exchanging lipids between patch ID = %d to patch_ID = %d \n',patch_1.ID,patch_2.ID);
            
            patch1_lipid_composition_up = patch_1.lipid_ratio_up;
            lipid_indexed_patch1 = [ ones( 1,patch1_lipid_composition_up(1) ) 2*ones( 1,patch1_lipid_composition_up(2) ) 3*ones( 1,patch1_lipid_composition_up(3) ) ];
            randomized_index_patch_1 = randperm(length(lipid_indexed_patch1));
            randomized_lipid_indexed_patch1 = lipid_indexed_patch1(randomized_index_patch_1);
            choosen_lipids_patch1 = randomized_lipid_indexed_patch1(1:quantal_change);
            if type_of_move == 1
                to_move_patch1_to_2 = ...
                    [ length(choosen_lipids_patch1(choosen_lipids_patch1==1)),...
                    length(choosen_lipids_patch1(choosen_lipids_patch1==2)),...
                    length(choosen_lipids_patch1(choosen_lipids_patch1==3)); ];
            else
                to_move_patch1_to_2 = randperm(quantal_change,2);to_move_patch1_to_2(3) = quantal_change - sum(to_move_patch1_to_2);
            end
            
            patch1_num_neighbours = length(patch_1.neighbours);
            patch1_lipid_pool_sharing = patch1_lipid_composition_up/patch1_num_neighbours;
            
            
            patch2_lipid_composition_up = patch_2.lipid_ratio_up;
            lipid_indexed_patch2 = [ ones( 1,patch2_lipid_composition_up(1) ) 2*ones( 1,patch2_lipid_composition_up(2) ) 3*ones( 1,patch2_lipid_composition_up(3) ) ];
            randomized_index_patch_2 = randperm(length(lipid_indexed_patch2));
            randomized_lipid_indexed_patch2 = lipid_indexed_patch2(randomized_index_patch_2);
            choosen_lipids_patch2 = randomized_lipid_indexed_patch2(1:quantal_change);
            if type_of_move == 1
                to_move_patch2_to_1 = ...
                    [ length(choosen_lipids_patch2(choosen_lipids_patch2==1)),...
                    length(choosen_lipids_patch2(choosen_lipids_patch2==2)),...
                    length(choosen_lipids_patch2(choosen_lipids_patch2==3)); ];
            else
                
                to_move_patch2_to_1 = randperm(quantal_change,2);to_move_patch2_to_1(3) = quantal_change - sum(to_move_patch2_to_1);
                
                
                
            end
            
            patch2_num_neighbours = length(patch_2.neighbours);
            patch2_lipid_pool_sharing = patch2_lipid_composition_up/patch2_num_neighbours;
            
            %             to_move_patch1_to_2 = randperm(quantal_change,2);to_move_patch1_to_2(3) = quantal_change - sum(to_move_patch1_to_2);
            %             to_move_patch2_to_1 = randperm(quantal_change,2);to_move_patch2_to_1(3) = quantal_change - sum(to_move_patch2_to_1);
            
            net_interchange_from_1_to_2 = -to_move_patch1_to_2 + to_move_patch2_to_1;
            
            
            fprintf('>>> Function: distribute_lipids; lipid_exchange_numbers = [%d %d %d]\n',...
                net_interchange_from_1_to_2(1), net_interchange_from_1_to_2(2), net_interchange_from_1_to_2(3) );
            
            fprintf('>>> Function: distribute_lipids; Exchanged \n ')
            %                 new_lipid_composition_patch1_up = patch1_lipid_composition_up - net_interchange_from_1_to_2;
            %                 new_lipid_composition_patch2_up = patch2_lipid_composition_up + net_interchange_from_1_to_2;
            %                 exchanged_noter = 1;
            if min( patch1_lipid_pool_sharing - net_interchange_from_1_to_2  )<0 || min( patch2_lipid_pool_sharing + net_interchange_from_1_to_2  )<0
                fprintf('>>> Function: distribute_lipids; NEGATIVE Number of lipids detected \n ')
                exchanged_noter = 0;
                new_lipid_composition_patch1_up = patch1_lipid_composition_up;
                new_lipid_composition_patch2_up = patch2_lipid_composition_up;
            else
                fprintf('>>> Function: distribute_lipids; Exchanged \n ')
                new_lipid_composition_patch1_up = patch1_lipid_composition_up - net_interchange_from_1_to_2;
                new_lipid_composition_patch2_up = patch2_lipid_composition_up + net_interchange_from_1_to_2;
                exchanged_noter = 1;
            end
        end
        
        
        function [new_lipid_composition_patch1_up,new_lipid_composition_patch2_up,exchanged_noter] = ...
                distribute_lipids_new(patch_1,patch_2,quantal_change)
                fprintf('>>> Function: distribute_lipids_new; Distributing lipids based on Jumping distribution\n ');
                % Select lipid patch to start with
                patch_to_select = randperm(2,1);
                
                lipid_packet = NaN(1,3);
                lipid_type_to_select = randperm(3);
                quantal_change_selected_lipid = ceil(abs(quantal_change*randn(1)));
                lipid_packet(lipid_type_to_select(1)) = quantal_change_selected_lipid;
                lipid_packet(lipid_type_to_select(2)) = -randperm(quantal_change_selected_lipid,1);
                lipid_packet(lipid_type_to_select(3)) = -(quantal_change_selected_lipid-abs(lipid_packet(lipid_type_to_select(2))));
                new_lipid_composition_patch1_up = patch_1.lipid_ratio_up;
                new_lipid_composition_patch2_up = patch_2.lipid_ratio_up; 
                
                if patch_to_select == 1                                        
                    new_lipid_composition_patch1_up = new_lipid_composition_patch1_up - lipid_packet;
                    new_lipid_composition_patch2_up = new_lipid_composition_patch2_up + lipid_packet;
                else                                       
                    new_lipid_composition_patch2_up = new_lipid_composition_patch2_up - lipid_packet;
                    new_lipid_composition_patch1_up = new_lipid_composition_patch1_up + lipid_packet;                    
                end
                
                if min(new_lipid_composition_patch1_up)<0  || min(new_lipid_composition_patch2_up)<0
                    new_lipid_composition_patch1_up = patch_1.lipid_ratio_up;
                    new_lipid_composition_patch2_up = patch_2.lipid_ratio_up; 
                    exchanged_noter = 0;
                else
                    exchanged_noter = 1;
                end
        end  
        
        function [new_lipid_composition_patch1_up,new_lipid_composition_patch2_up,exchanged_noter] = ...
                distribute_lipids_new2(patch_1,patch_2,quantal_change,diffusion_rates_ratio)
            
            fprintf('>>> Function: distribute_lipids_new2; Exchanging lipids between patch ID = %d to patch_ID = %d \n',patch_1.ID,patch_2.ID);
            
            patch1_lipid_composition_up = patch_1.lipid_ratio_up;
            patch2_lipid_composition_up = patch_2.lipid_ratio_up;
            
            change_lipid_comps = ceil(quantal_change*diffusion_rates_ratio.*randn(1,3));
            new_lipid_composition_patch1_up = patch1_lipid_composition_up  + change_lipid_comps;
            new_lipid_composition_patch2_up = patch2_lipid_composition_up  - change_lipid_comps;
            
            if min( [ new_lipid_composition_patch1_up new_lipid_composition_patch2_up ] )<=0
                exchanged_noter = 0;
                new_lipid_composition_patch1_up = patch1_lipid_composition_up;
                new_lipid_composition_patch2_up = patch2_lipid_composition_up;
            else                
                exchanged_noter = 1;
            end
            
        end
        
        function [triangle_color_mole_fraction,triangle_color_curvature] = color_code_lipid_density_profile_patch_wise_and_curvature(obj_list)
            %             rgb_color_mat = NaN(size(obj_list,1),3);
            fprintf('>>> Function: color_code_lipid_density_profile_patch_wise; RGB Color coding patches based on mole fraction\n');
            [lipid_comp] = extract_lipids_composition_ID_specific(obj_list,1);
            total_num_lipid_array = sum(lipid_comp,2);
            fprintf('>>> Function: color_code_lipid_density_profile_patch_wise; Calculating the mole fraction\n');
            lipid_comp_mole_fraction = [ lipid_comp(:,1)./total_num_lipid_array lipid_comp(:,2)./total_num_lipid_array lipid_comp(:,3)./total_num_lipid_array  ];
            fprintf('>>> Function: color_code_lipid_density_profile_patch_wise; Generating surface mesh for mesh coloring\n');
            [surface_mesh] = generate_mesh_from_objlist(obj_list);
            triangles = surface_mesh.ConnectivityList;
            fprintf('>>> Function: color_code_lipid_density_profile_patch_wise; Extracting number of neighbours for weighting \n');
            [neighbours_num_all] = extract_num_neighbours_ID_specific(obj_list);
            triangle_color_mole_fraction = NaN(size(triangles,1),3);
            triangle_color_curvature = NaN(size(triangles,1),1);
            fprintf('>>> Function: color_code_lipid_density_profile_patch_wise; Generating mesh color \n');
            
            [mean_curvature,~] = extract_curvature_list(obj_list);
            for tri_ind = 1:size(triangles,1)
                current_triangle = triangles(tri_ind,:);
                lipid_comp_of_vertices = lipid_comp(current_triangle,:);
                weights = neighbours_num_all(current_triangle);
                weighted_lipid_comp = [ lipid_comp_of_vertices(:,1).*weights lipid_comp_of_vertices(:,2).*weights lipid_comp_of_vertices(:,3).*weights ];
                total_lipids_in_current_triangle = sum(weighted_lipid_comp(:));
                lipid_comp_type_wise = sum(weighted_lipid_comp,1);
                lipid_comp_type_wise_mole_fraction = lipid_comp_type_wise/total_lipids_in_current_triangle;
                triangle_color_mole_fraction(tri_ind,:) = lipid_comp_type_wise_mole_fraction;
                triangle_color_curvature(tri_ind,:) = mean(mean_curvature(current_triangle));
            end
            fprintf('>>> Function: color_code_lipid_density_profile_patch_wise; Generating mesh color Done \n\n');
            
            
        end
        
        function [triangle_color_mole_fraction,triangle_color_curvature] = color_code_lipid_density_profile_patch_wise_and_curvature_new(obj_list)
            %             rgb_color_mat = NaN(size(obj_list,1),3);
            fprintf('>>> Function: color_code_lipid_density_profile_patch_wise; RGB Color coding patches based on mole fraction\n');
            [lipid_comp] = extract_lipids_composition_ID_specific(obj_list,1);
            total_num_lipid_array = sum(lipid_comp,2);
            fprintf('>>> Function: color_code_lipid_density_profile_patch_wise; Calculating the mole fraction\n');
            lipid_comp_mole_fraction = [ lipid_comp(:,1)./total_num_lipid_array lipid_comp(:,2)./total_num_lipid_array lipid_comp(:,3)./total_num_lipid_array  ];
            fprintf('>>> Function: color_code_lipid_density_profile_patch_wise; Generating surface mesh for mesh coloring\n');
            [surface_mesh] = generate_mesh_from_objlist(obj_list);
            triangles = surface_mesh.ConnectivityList;
            Points_all = surface_mesh.Points;
            fprintf('>>> Function: color_code_lipid_density_profile_patch_wise; Extracting number of neighbours for weighting \n');
            [neighbours_num_all] = extract_num_neighbours_ID_specific(obj_list);
            triangle_color_mole_fraction = NaN(size(triangles,1),3);
            triangle_color_curvature = NaN(size(triangles,1),1);
            fprintf('>>> Function: color_code_lipid_density_profile_patch_wise; Generating mesh color \n');
            all_face_areas = @(ind) obj_list(ind).Av_vertex;
            all_areas = arrayfun(all_face_areas,[1:size(obj_list)]);
            
            
            [mean_curvature,~] = extract_curvature_list(obj_list);
            for tri_ind = 1:size(triangles,1)
                current_triangle = triangles(tri_ind,:);
                current_triangle_points = Points_all(current_triangle,:);
                area_current_triangle  = area_triangle(current_triangle_points);
                
                lipid_comp_of_vertices = lipid_comp(current_triangle,:);
                
                vertex_1id = current_triangle(1);
                vertex_1id_area_represented = all_areas(vertex_1id);
                lipid_comp_of_vertices1 = ...
                    1/3 * area_current_triangle/vertex_1id_area_represented * lipid_comp(vertex_1id,:);
                
                vertex_2id = current_triangle(2);
                vertex_2id_area_represented = all_areas(vertex_2id);
                lipid_comp_of_vertices2 = ...
                    1/3 * area_current_triangle/vertex_2id_area_represented * lipid_comp(vertex_2id,:);
                
                vertex_3id = current_triangle(3);
                vertex_3id_area_represented = all_areas(vertex_3id);
                lipid_comp_of_vertices3 = ...
                    1/3 * area_current_triangle/vertex_3id_area_represented * lipid_comp(vertex_3id,:);
                
                triangle_lipid_comp = [ lipid_comp_of_vertices1;...
                                        lipid_comp_of_vertices2;...
                                        lipid_comp_of_vertices3];
                sum_triangle_lipid = sum(triangle_lipid_comp,1);
                total_lipids = sum(sum_triangle_lipid);
                sum_triangle_lipid_mole_fraction = sum_triangle_lipid/total_lipids;
                
                
                weights = neighbours_num_all(current_triangle);
                weighted_lipid_comp = [ lipid_comp_of_vertices(:,1).*weights lipid_comp_of_vertices(:,2).*weights lipid_comp_of_vertices(:,3).*weights ];
                total_lipids_in_current_triangle = sum(weighted_lipid_comp(:));
                lipid_comp_type_wise = sum(weighted_lipid_comp,1);
                lipid_comp_type_wise_mole_fraction = lipid_comp_type_wise/total_lipids_in_current_triangle;
                triangle_color_mole_fraction(tri_ind,:) = sum_triangle_lipid_mole_fraction;
                triangle_color_curvature(tri_ind,:) = mean(mean_curvature(current_triangle));
            end
            fprintf('>>> Function: color_code_lipid_density_profile_patch_wise; Generating mesh color Done \n\n');
            
            
        end
        
        function [color_patch] = color_code_patch_mole_fraction_for_scatter_plot(obj_list)
            fprintf('>>> Function: color_code_patch_mole_fraction_for_scatter_plot; RGB Color coding patches based on mole fraction\n');
            [lipid_comp] = extract_lipids_composition_ID_specific(obj_list,1);
            total_num_lipid_array = sum(lipid_comp,2);
            fprintf('>>> Function: color_code_patch_mole_fraction_for_scatter_plot; Calculating the mole fraction\n');
            color_patch = [ lipid_comp(:,1)./total_num_lipid_array lipid_comp(:,2)./total_num_lipid_array lipid_comp(:,3)./total_num_lipid_array  ];
            
        end
        
        function [mean_curvature,principle_curvatures] = extract_curvature_list(obj_list)
            principle_curvatures = NaN(size(obj_list,1),2);
            extract_c1_func = @(id) obj_list(id).c1;
            principle_curvatures(:,1) = arrayfun(extract_c1_func,1:size(obj_list,1));
            extract_c2_func = @(id) obj_list(id).c2;
            principle_curvatures(:,2) = arrayfun(extract_c2_func,1:size(obj_list,1));
            mean_curvature = mean(principle_curvatures,2);
        end
        
        function [obj_list] = distribute_lipids_randomly_global(obj_list)
            fprintf('>>> Function: distribute_lipids_randomly_global; Distributing lipids in each patch randomly \n');
            for obj_ind = 1:size(obj_list,1)
                current_obj = obj_list(obj_ind);
                avg_lipid_proportionality = randperm(100,3);
                [current_obj] = randomly_distribute_lipids(current_obj,avg_lipid_proportionality);
                obj_list(obj_ind) = current_obj;
            end
            %             [obj_list] = randomly_distribute_lipids(obj_list,avg_lipid_proportionality)
            
        end
        
        function [patch_pos_list,curvature_energy_list,entropy_of_mixing_list,surface_stretching_energy_list,distortion_energy_list,H_mean_list,H_spontaneous_list] = ...
                track_patch(obj_list,type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous,quantal_change,each_mole_content)
            patch_pos_list = NaN(size(obj_list,1),3);
            curvature_energy_list = NaN(size(obj_list,1),1);
            entropy_of_mixing_list = NaN(size(obj_list,1),1);
            surface_stretching_energy_list = NaN(size(obj_list,1),1);
            distortion_energy_list = NaN(size(obj_list,1),1);
            H_mean_list = NaN(size(obj_list,1),1);
            H_spontaneous_list = NaN(size(obj_list,1),1);
            for ob_no = 1:size(obj_list,1)
                obj = obj_list(ob_no,1);
                patch_pos_list(ob_no,:) = obj.Pos;
                [curvature_energy_list(ob_no),entropy_of_mixing_list(ob_no),surface_stretching_energy_list(ob_no),distortion_energy_list(ob_no)] =...
                    coarse_grain_energy_mem_patch(obj,type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous,quantal_change,each_mole_content);
                H_mean_list(ob_no,:) = (obj.c1 + obj.c2)/2  ;
                H_spontaneous_list(ob_no,:) =...
                    sum((obj.lipid_ratio_up/sum(obj.lipid_ratio_up)).*type_properties(2,:));
            end
        end
        
        function [c_preferred,c_current,preffered_mole_fractions,actual_mole_fractions,c1_current,c2_current] = ...
                calculate_equillibrium_characteristics(obj,type_properties,temperature,plot_profile,each_mole_content)
            c1_current = obj.c1;
            c2_current = obj.c2;
            c_current = (obj.c1 + obj.c2)/2;
            c_p = type_properties(2,1);
            c_0 = type_properties(2,2);
            c_n = type_properties(2,3);
            k = mean(type_properties(1,:));
            
            [~,face_areas,~] = extract_obj_star(obj);
            total_area = sum(face_areas)/3;
            N_particles = obj.n_particles_up;
            area_of_each_lipid = obj.avg_area_per_lipid;
            %             hamiltonian = @(x,y) 0.5.*k*N_particles*area_of_each_lipid*((2*c_current - (c_p.*x  + c_n.*y + c_0*(1 - x -y))).^2);
            hamiltonian = @(x,y) 0.5.*k*N_particles*area_of_each_lipid*( (c1_current - (c_p.*x  + c_n.*y + c_0*(1 - x -y))).^2 + (c2_current - (c_p.*x  + c_n.*y + c_0*(1 - x -y))).^2 );
            entropy = @(x,y) -N_particles*8.314*temperature*(x.*log(x) + y.*log(y) + (1-x-y).*log(1-x-y))/each_mole_content;
            
            [X_grid,Y_grid]  = meshgrid((0:0.005:1),(0:0.005:1));
            NaN_grid = NaN(size(X_grid));
            NaN_grid(1-(X_grid+Y_grid)>=0)=1;
            energy_well = (hamiltonian(X_grid,Y_grid));
            probability = exp(-hamiltonian(X_grid,Y_grid)/temperature);
            entropy_well = real(entropy(X_grid,Y_grid));
            
            equillibrium_landscape = exp( (entropy_well)/temperature).*probability.*NaN_grid;
            %             equillibrium_landscape = exp( (entropy_well-hamiltonian(X_grid,Y_grid))/temperature);
            [~,ind] = max(equillibrium_landscape(:));
            
            x_fraction = X_grid(ind);
            y_fraction = Y_grid(ind);
            
            c_preferred = x_fraction*c_p + y_fraction*c_n + (1-x_fraction-y_fraction)*c_0;
            actual_mole_fractions = obj.lipid_ratio_up./sum(obj.lipid_ratio_up);
            preffered_mole_fractions = [x_fraction 1-x_fraction-y_fraction y_fraction ];
            if plot_profile == 1
                surf(X_grid,Y_grid,equillibrium_landscape);daspect([1,1,1]);
                % surf(X_grid,Y_grid,entropy_well);hold off
                xlabel('positive curvature');ylabel('negative curvature'); view(2)
            end
        end
        
        function [interesting_ids] = findobj_at_pos(obj_list,pos,rad)
            [IDs] = extract_IDs(obj_list);
            [coors] = extract_coors_points_ID_specific(obj_list);
            reference_coor_mat = (pos'*ones(1,size(coors,1)))';
            pos_vects = coors - reference_coor_mat;
            eucleid_dis = vecnorm(pos_vects,2,2);
            interesting_ids = IDs(eucleid_dis<=rad);
            interesting_eucleid_dis = eucleid_dis(eucleid_dis<=rad);
            [~,sorted_id] = sort(interesting_eucleid_dis);
            interesting_ids = interesting_ids(sorted_id);
        end
        
        function [c1,c2,r,Normal_vect,Av,updated_obj,principal_directions_current] = determine_principal_curvatures(obj)
            % Using curvature estimation method by Taubin
            [star_mesh,face_areas] = extract_obj_only_star(obj);
            faces_current = star_mesh.ConnectivityList;
            face_normals = faceNormal(star_mesh);
            weight_face_areas = face_areas/sum(face_areas);
            Av = sum(face_areas)/3;
            area_based_weight_mat = weight_face_areas*ones(1,3);
            weighted_face_normals = face_normals.*area_based_weight_mat;
            Not_normalized_normal = sum(weighted_face_normals,1);
            Normal_vect = Not_normalized_normal/vecnorm(Not_normalized_normal);
            center_coor = star_mesh.Points(1,:);
            neighbour_coors = star_mesh.Points(2:end,:);
            
            % Calculating the Tij
            I = diag([1,1,1]);
            projection_mat = (I - Normal_vect'*Normal_vect );
            Tij_all = NaN(3,size(neighbour_coors,1));
            kij_all = NaN(size(neighbour_coors,1),1);
            wij_all = NaN(size(neighbour_coors,1),1);
            Mvi_mat = NaN(3,3,size(neighbour_coors,1));
            for j = 1:size(neighbour_coors,1)
                current_neigh = neighbour_coors(j,:);
                real_neighbour_ind = j + 1;
                edge_vect = current_neigh - center_coor;
                
                Tij = projection_mat* edge_vect';
                Tij = Tij/vecnorm(Tij);
                Tij_all(:,j) = Tij;
                
                kij = 2*(Normal_vect*edge_vect')/(vecnorm(edge_vect)^2);
                kij_all(j) = kij;
                
                adjacent_face_mask = double(ismember( faces_current,[1,real_neighbour_ind]));
                adjacent_face_mask_1D = sum(adjacent_face_mask,2)==2;
                adjacent_face_areas = face_areas(adjacent_face_mask_1D);
                wij = sum(adjacent_face_areas);
                wij_all(j) =  wij;
                
                Mvi = wij*kij*( Tij*Tij' );
                Mvi_mat(:,:,j) = Mvi;
            end
            wij_all_normalizing_fact = sum(wij_all);
            Mvi_mat = Mvi_mat./wij_all_normalizing_fact;
            
            Mvi_final = sum(Mvi_mat,3);
            
            % Calculating Householder matrix
            E1 = [1,0,0]';
            condition_check_measure_1 = vecnorm( (E1(:) - Normal_vect(:)) );
            condition_check_measure_2 = vecnorm( (E1(:) + Normal_vect(:)) );
            if condition_check_measure_1>condition_check_measure_2
                Wvi =  E1 - Normal_vect';
                Wvi = Wvi/vecnorm(Wvi);
            else
                Wvi =  E1 + Normal_vect';
                Wvi = Wvi/vecnorm(Wvi);
            end
            
            Qvi = I - 2*(Wvi*Wvi'); % the famous Housholder matrix to allign to Normal vector
            Mvi_transformed = Qvi'*(Mvi_final*Qvi);
            Mvi_minor = -Mvi_transformed(2:end,2:end);
            mp11 = Mvi_minor(1,1);
            mp22 = Mvi_minor(2,2);
            c1 = (3*mp11 - mp22);
            c2 = (3*mp22 - mp11);
            
            principal_curvatures = [c1 c2];
            theta_eigen_direction = ...
                0.5*atan(2*Mvi_minor(1,2)/( Mvi_minor(1,1)-Mvi_minor(2,2) )); %https://arxiv.org/pdf/1306.6291.pdf
            eigen_direction1 = ...
                cos(theta_eigen_direction)*Qvi(:,2) - sin(theta_eigen_direction)*Qvi(:,3);
            eigen_direction2 = ...
                sin(theta_eigen_direction)*Qvi(:,2) + cos(theta_eigen_direction)*Qvi(:,3);
            principal_directions_current = [eigen_direction1,eigen_direction2];
            
            [principal_curvatures,ind] = sort(principal_curvatures);
            principal_directions_current = principal_directions_current(:,ind);
            c1 = principal_curvatures(1);
            c2 = principal_curvatures(2);
            
            r = 1/mean([c1,c2]);
            updated_obj = obj;
            updated_obj.Normal_vect = Normal_vect;
            updated_obj.c1 = c1;
            updated_obj.c2 = c2;
            updated_obj.principal_directions = principal_directions_current;
            updated_obj.Av_vertex = Av;
        end  %%***************
        
        function [principal_curvatures_neighbours,spontaneous_curvature_neigh,N_particles_neigh,obj_neighbours_updated] = determine_principle_curvatures_at_neighbours_dis_at_center(updated_obj,obj_neighbours,type_curvatures)
            % Algorithm: update the new center position in each neighbour
            % Calculate the geometrical properties of each neighbour            
            center_ID = updated_obj.ID;
            center_pos = updated_obj.Pos;
            obj_neighbours_updated = obj_neighbours;
            fprintf('>>> Function:determine_principle_curvatures_at_neighbours_displacement_at_center; updating neighbours\n');
            principal_curvatures_neighbours =  [ NaN(size(obj_neighbours_updated,1),1) NaN(size(obj_neighbours_updated,1),1) ];
            spontaneous_curvature_neigh = NaN(size(obj_neighbours_updated,1),1);
            N_particles_neigh = NaN(size(obj_neighbours_updated,1),1);
            for neighbour_id = 1:size(obj_neighbours,1)
                current_neighbour = obj_neighbours(neighbour_id);
                center_ID_mask_in_neighbour = current_neighbour.neighbours==center_ID;
                current_neighbour.neighbours_coor(center_ID_mask_in_neighbour,:) = center_pos;
                [c1_current,c2_current,r,Normal_vect,Av,updated_neigh] = determine_principal_curvatures(current_neighbour);
                principal_curvatures_neighbours(neighbour_id,:) = [c1_current,c2_current];
                N_particles_neigh(neighbour_id) = sum(current_neighbour.lipid_ratio_up);
                mole_fraction = current_neighbour.lipid_ratio_up/sum(current_neighbour.lipid_ratio_up);
                spontaneous_curvature_neigh(neighbour_id) = sum(type_curvatures.*mole_fraction);
                
                obj_neighbours_updated(neighbour_id) = updated_neigh;
                
            end           
            
        end
        
        function [] = check_curvature_characteristics(obj_list,figure_handle)
            [obj_list] = derive_geometrical_quatities_all(obj_list);
            [surface_mesh] = generate_mesh_from_objlist(obj_list);
        end
        
        function [del_entropy] = calculate_change_in_entropy(updated_obj_pair,old_obj_pair,each_mole_content,temperature)
            % calculating change in entropy inspired from ideal gas mixing
            box1_new = updated_obj_pair(1);
            box2_new = updated_obj_pair(2);
            
            area_box1_new = box1_new.Av_vertex;
            area_box2_new = box2_new.Av_vertex;
            
            N_lipids_total_box1_new = box1_new.n_particles_up;
            N_lipids_total_box2_new = box2_new.n_particles_up;
            
            
            lipids_comp_box1_new = box1_new.lipid_ratio_up;
            mole_fraction_box1_new = lipids_comp_box1_new./sum(lipids_comp_box1_new);
            lipids_comp_box2_new = box2_new.lipid_ratio_up;
            mole_fraction_box2_new = lipids_comp_box2_new./sum(lipids_comp_box2_new);
            
            partial_pressure_new_box1 = mole_fraction_box1_new*N_lipids_total_box1_new/area_box1_new;
            partial_pressure_new_box2 = mole_fraction_box2_new*N_lipids_total_box2_new/area_box2_new;
            
            
            box1_old = old_obj_pair(1);
            box2_old = old_obj_pair(2);
            
            area_box1_old = box1_old.Av_vertex;
            area_box2_old = box2_old.Av_vertex;
            
            N_lipids_total_box1_old = box1_old.n_particles_up;
            N_lipids_total_box2_old = box2_old.n_particles_up;
            
            lipids_comp_box1_old = box1_old.lipid_ratio_up;
            mole_fraction_box1_old = lipids_comp_box1_old./sum(lipids_comp_box1_old);
            lipids_comp_box2_old = box2_old.lipid_ratio_up;
            mole_fraction_box2_old = lipids_comp_box2_old./sum(lipids_comp_box2_old);
            
            partial_pressure_old_box1 = mole_fraction_box1_old*N_lipids_total_box1_old/area_box1_old;
            partial_pressure_old_box2 = mole_fraction_box2_old*N_lipids_total_box2_old/area_box2_old;
            
            % Calculating entropy for three kinds of lipids
            %             del_entropy = @( no_molecules_f,no_molecules_i, p_f, p_i) -1/each_mole_content.*8.314.*log( (p_f^no_molecules_f)/(p_i^no_molecules_i));
            del_entropy = @( no_molecules_f,no_molecules_i, p_f, p_i) -no_molecules_f/each_mole_content.*8.314.*log(p_f)+no_molecules_i/each_mole_content.*8.314.*log(p_i);
            
            
            no_molecules_p_box1_f = lipids_comp_box1_new(1);
            no_molecules_p_box1_i = lipids_comp_box1_old(1);
            del_entropy_p_box1 = del_entropy( no_molecules_p_box1_f,no_molecules_p_box1_i,partial_pressure_new_box1(1),partial_pressure_old_box1(1) );
            
            no_molecules_p_box2_f = lipids_comp_box2_new(1);
            no_molecules_p_box2_i = lipids_comp_box2_old(1);
            del_entropy_p_box2 = del_entropy( no_molecules_p_box2_f,no_molecules_p_box2_i,partial_pressure_new_box2(1),partial_pressure_old_box2(1) );
            
            del_entropy_p = del_entropy_p_box1 + del_entropy_p_box2;
            
            no_molecules_o_box1_f = lipids_comp_box1_new(2);
            no_molecules_o_box1_i = lipids_comp_box1_old(2);
            del_entropy_o_box1 = del_entropy( no_molecules_o_box1_f,no_molecules_o_box1_i,partial_pressure_new_box1(2),partial_pressure_old_box1(2) );
            
            no_molecules_o_box2_f = lipids_comp_box2_new(2);
            no_molecules_o_box2_i = lipids_comp_box2_old(2);
            del_entropy_o_box2 = del_entropy( no_molecules_o_box2_f,no_molecules_o_box2_i,partial_pressure_new_box2(2),partial_pressure_old_box2(2) );
            
            del_entropy_o = del_entropy_o_box1 + del_entropy_o_box2;
            
            no_molecules_n_box1_f = lipids_comp_box1_new(3);
            no_molecules_n_box1_i = lipids_comp_box1_old(3);
            del_entropy_n_box1 = del_entropy( no_molecules_n_box1_f,no_molecules_n_box1_i,partial_pressure_new_box1(3),partial_pressure_old_box1(3) );
            
            no_molecules_n_box2_f = lipids_comp_box2_new(3);
            no_molecules_n_box2_i = lipids_comp_box2_new(3);
            del_entropy_n_box2 = del_entropy( no_molecules_n_box2_f,no_molecules_n_box2_i,partial_pressure_new_box2(3),partial_pressure_old_box2(3) );
            
            del_entropy_n = del_entropy_n_box1 + del_entropy_n_box2;
            
            %% rewriting again
            box1_new_entropy_all_lipids = -mole_fraction_box1_new.*log(partial_pressure_new_box1)*8.314;
            box2_new_entropy_all_lipids = -mole_fraction_box2_new.*log(partial_pressure_new_box2)*8.314;
            box1_old_entropy_all_lipids = -mole_fraction_box1_old.*log(partial_pressure_old_box1)*8.314;
            box2_old_entropy_all_lipids = -mole_fraction_box2_old.*log(partial_pressure_old_box2)*8.314;
            
            initial_entropy = box1_old_entropy_all_lipids + box2_old_entropy_all_lipids;
            final_entropy   = box1_new_entropy_all_lipids + box2_new_entropy_all_lipids;
            
            %             entropy_calc = @(area_fraction,no_lipids,total_number_of_lipids) -area_fraction.*( no_lipids.*8.314.*log(no_lipids/total_number_of_lipids)/each_mole_content );
            del_entropy = -(final_entropy - initial_entropy)*temperature;
            %             del_entropy = del_entropy_p + del_entropy_o + del_entropy_n;
            if del_entropy>0
                disp('gotcha')
            end
        end
        
        function [del_entropy] = calculate_change_in_entropy_configurationBased(updated_obj,old_obj,each_mole_content)
            
            % based on change in entropy due to change in partial pressure
            box1_new = updated_obj;
            box1_old = old_obj;
            
            area_box1_new = box1_new.Av_vertex;
            area_box1_old = box1_old.Av_vertex;
            
            N_lipids_total_box1_new = box1_new.n_particles_up;
            N_lipids_total_box1_old = box1_old.n_particles_up;
            
            lipids_comp_box1_new = box1_new.lipid_ratio_up;
            mole_fraction_box1_new = lipids_comp_box1_new./sum(lipids_comp_box1_new);
            
            lipids_comp_box1_old = box1_old.lipid_ratio_up;
            mole_fraction_box1_old = lipids_comp_box1_old./sum(lipids_comp_box1_old);
            
            partial_pressure_new_box1 = mole_fraction_box1_new*N_lipids_total_box1_new/area_box1_new;
            partial_pressure_old_box1 = mole_fraction_box1_old*N_lipids_total_box1_old/area_box1_old;
            
            del_entropy = @( no_molecules, p_f, p_i) -no_molecules/each_mole_content.*8.314.*log(p_f/p_i);
            no_molecules = N_lipids_total_box1_new; % also = N_lipids_total_box1_old
            
            del_entropy = del_entropy(no_molecules,sum(partial_pressure_new_box1),sum(partial_pressure_old_box1));
            
        end
        
        function [obj_pair] = simulate_diffusion_in_pair(obj_pair,type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous,quantal_change,diffusion_rates_ratio,each_mole_content)
            for t = 1:500
                patch_1 = obj_pair(1);
                patch_2 = obj_pair(2);
                updated_patch_1 = patch_1;
                updated_patch_2 = patch_2;
%                 [new_lipid_composition_patch1_up,new_lipid_composition_patch2_up,exchanged_noter] = ...
%                     distribute_lipids(patch_1,patch_2,quantal_change,2);
                
                [new_lipid_composition_patch1_up,new_lipid_composition_patch2_up,exchanged_noter] = ...
                distribute_lipids_new2(patch_1,patch_2,quantal_change,diffusion_rates_ratio);
                updated_patch_1.lipid_ratio_up = new_lipid_composition_patch1_up;
                updated_patch_1.neighbours_lipid_compositions_up(updated_patch_1.neighbours==updated_patch_2.ID,:) = new_lipid_composition_patch2_up;
                
                updated_patch_2.lipid_ratio_up = new_lipid_composition_patch2_up;
                updated_patch_2.neighbours_lipid_compositions_up(updated_patch_2.neighbours==updated_patch_1.ID,:) = new_lipid_composition_patch1_up;
                
                updated_patch_pair = [updated_patch_1,updated_patch_2];
                %                 current_patch_pair_ordered = [patch_1,patch_2];
                [decision] = ...
                    metropolis_lipid_exchange_displacement(updated_patch_pair,obj_pair,type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous,quantal_change,each_mole_content);
                u = rand(1);
                if u <= decision*exchanged_noter
                    obj_pair = updated_patch_pair;
                end
                fprintf('object1 lipid comp\t');
                sum(obj_pair(1).lipid_ratio_up)
                fprintf('\n');
                
                fprintf('object2 lipid comp\t');
                sum(obj_pair(2).lipid_ratio_up)
                fprintf('\n');
                
                C1_p = obj_pair(1).lipid_ratio_up(1);%/sum(obj_pair(1).lipid_ratio_up(:));
                C2_p = obj_pair(2).lipid_ratio_up(1);%/sum(obj_pair(2).lipid_ratio_up(:));
                C1_o = obj_pair(1).lipid_ratio_up(2);%/sum(obj_pair(1).lipid_ratio_up(:));
                C2_o = obj_pair(2).lipid_ratio_up(2);%/sum(obj_pair(2).lipid_ratio_up(:));
                C1_n = obj_pair(1).lipid_ratio_up(3);%/sum(obj_pair(1).lipid_ratio_up(:));
                C2_n = obj_pair(2).lipid_ratio_up(3);%/sum(obj_pair(2).lipid_ratio_up(:));
                
                %                 scatter3( 0,0,0,500,C1_p,'filled');hold on;scatter3( .1,0,0,500,C1_o,'filled');hold on;scatter3( .2,0,0,500,C1_n,'filled');
                %                 hold on; scatter3( 0,.1,0,500,C2_p,'filled');hold on; scatter3( 0.1,0.1,0,500,C2_o,'filled');hold on; scatter3( 0.2,.1,0,500,C2_n,'filled');
                %                 hold off;
                caxis([0,.5]); daspect([1,1,1]);xlim([-0.2,0.3]);
                spontaneous_curv_1 = sum((obj_pair(1).lipid_ratio_up./sum(obj_pair(1).lipid_ratio_up) ).*type_properties(2,:));
                spontaneous_curv_2 = sum( (obj_pair(2).lipid_ratio_up./sum(obj_pair(2).lipid_ratio_up) ).*type_properties(2,:));
                bar([C1_p,C1_o,C1_n],'r','FaceAlpha',.5);hold on; bar([C2_p,C2_o,C2_n],'b','FaceAlpha',.5);ylim([0,2000]);hold off
                %                 view(2)
%                 current_frame = getframe(gcf);
%                 imwrite(current_frame.cdata,'diffusion_obj_pair40000_reveresed.tif','WriteMode','append');
                pause(.1);
                
                
            end
        end
        
        function [H_mean,H_spontaneous,principal_curv_mat,mole_fraction] = curvature_histogram(obj_list,type_properties,Binlimits,NumBins,plot_or_not)
            fprintf('>>>Function: curvature_histogram; Calculating principal curvatures at each point \n');
            [obj_list] = derive_geometrical_quatities_all(obj_list);
%             c1_noter = NaN(size(obj_list,1),1);
%             c2_noter = NaN(size(obj_list,1),1);
            extract_c1_func = @(id) obj_list(id).c1;
            extract_c2_func = @(id) obj_list(id).c2;
            
            extract_cp_func = @(id) obj_list(id).lipid_ratio_up(1);
            extract_co_func = @(id) obj_list(id).lipid_ratio_up(2);
            extract_cn_func = @(id) obj_list(id).lipid_ratio_up(3);
            
            c1_noter = arrayfun(extract_c1_func,1:size(obj_list,1));
            c2_noter = arrayfun(extract_c2_func,1:size(obj_list,1));
            
            cp_noter = arrayfun(extract_cp_func,1:size(obj_list,1));
            co_noter = arrayfun(extract_co_func,1:size(obj_list,1));
            cn_noter = arrayfun(extract_cn_func,1:size(obj_list,1));
            lipid_ratio_up_all = [ cp_noter' co_noter' cn_noter' ];
            mole_fraction = [ cp_noter'./sum(lipid_ratio_up_all,2) co_noter'./sum(lipid_ratio_up_all,2) cn_noter'./sum(lipid_ratio_up_all,2) ];
            
            principal_curv_mat = [c1_noter',c2_noter'];
            H_mean = mean(principal_curv_mat,2);
            
            lipid_spontaneous_curvatures = type_properties(2,:);
            lipid_spontaneous_curvatures_mat = (lipid_spontaneous_curvatures'*ones(1,size(obj_list,1)))';
            
            H_spontaneous = sum(mole_fraction.*lipid_spontaneous_curvatures_mat,2);
%             hist(axes(fig_handle),[H_mean,H_spontaneous],20);
            
            BinWidth = abs(diff(Binlimits))/NumBins;
            BinEdges = Binlimits(1)-.5*BinWidth:BinWidth:Binlimits(2)+.5*BinWidth;
            Bin_centers = Binlimits(1):BinWidth:Binlimits(2);
            h_counts_hmean = histcounts(H_mean,BinEdges);
%             h_counts_hmean = hist_Hmean_handle.Values;
           h_counts_hspontaneous = histcounts(H_spontaneous,BinEdges);
%             close(hist_HSpontaneous_handle);
%             h_counts_hspontaneous = hist_HSpontaneous_handle.Values;
            
%             axes(fig_handle)
            if plot_or_not ==1
                bar(Bin_centers,h_counts_hmean,'r','FaceAlpha',.5);hold on
                bar(Bin_centers,h_counts_hspontaneous,'g','FaceAlpha',.5);hold off
                xlabel('curvature');ylabel('counts');title('RED =  Hmean(Geometry); GREEN = HSpont(lipids)')
            end
        end
        
        function [curvature_cost_function,diffusion_cost_function] = check_equillibrium_statistics(obj_list,type_properties,plot_yes_or_not)
            [H_mean,H_spontaneous,principal_curv_mat,mole_fraction] = curvature_histogram(obj_list,type_properties,[-1,1],20,0);
            if plot_yes_or_not == 1
            plot(H_mean,H_spontaneous,'*');daspect([1,1,1]); xlabel('Hmean'); ylabel('Hspontaneous');
            end
            
            curvature_cost_function = sqrt(mean((H_mean - 2*H_spontaneous).^2))/mean(H_spontaneous)*100;
            diffusion_cost_function = (std(mole_fraction(:,1)) + std(mole_fraction(:,2)) + std(mole_fraction(:,3)))/mean(mole_fraction(:))*100;
            
        end
        
        function [del_t_diffuse,del_t_move,quantal_change,kick,hypothetical_radius] = predict_parameters(obj_list,diffusion_rate,type_properties,surface_tension,viscosity,temperature)
             
            extract_area_patch_func = @(id) obj_list(id).Av_vertex;
             Av_vertex_all = arrayfun(extract_area_patch_func,1:size(obj_list,1));
             Av_vertex_all_mean = mean(Av_vertex_all);
             hypothetical_radius = sqrt(Av_vertex_all_mean/pi);
             flux_boundary = 2*pi*hypothetical_radius;
             del_t = 1*Av_vertex_all_mean/4/diffusion_rate;
             del_t_diffuse = del_t;
             
             [H_mean,H_spontaneous,principal_curv_mat,mole_fraction] = curvature_histogram(obj_list,type_properties,[-1,1],20,0);
             [lipid_comp] = extract_lipids_composition_ID_specific(obj_list,1);
             std_mole_fraction = mean(std(lipid_comp,1,1));
             gradient_conc = std_mole_fraction/hypothetical_radius/1;
             quantal_change = ceil(gradient_conc*diffusion_rate*2*pi*hypothetical_radius*del_t);
%              quantal_change = ceil(2*pi*hypothetical_radius*std_mole_fraction*del_t/hypothetical_radius/2);
            
%             kick = temperature*Av_vertex_all_mean/16/pi.^3/mean(type_properties(1,:));
%             kick = hypothetical_radius/5;
            q = 1/hypothetical_radius;
            kappa = mean(type_properties(1,:));
            
            relaxation_freq_q = (kappa*q^3 + surface_tension*q)/4/viscosity;
            del_t_move = 1/relaxation_freq_q;
            kick = sqrt(temperature/( kappa*q^4 + surface_tension*q^2));
        end
        
        function [] = plot_curvature_profile(obj_list,num_fig,type_properties,Binlimits,NumBins)
%             [patch_color,triangle_color_curvature] = color_code_lipid_density_profile_patch_wise_and_curvature(obj_list);
            h1 = figure(num_fig);
            h1_sub1 = subplot(2,1,1);
            
            [surface_mesh] = generate_mesh_from_objlist(obj_list); 
            [patch_color,triangle_color_curvature] = color_code_lipid_density_profile_patch_wise_and_curvature(obj_list);    
            axes(h1_sub1);trisurf(surface_mesh,triangle_color_curvature,'FaceAlpha',.3,'EdgeAlpha',0.2); daspect([1,1,1]);title(h1_sub1,'Curvature');colormap(h1_sub1,'jet');
            view([0,0]);hold off;
            h1_sub2 = subplot(2,1,2);
            [H_mean,H_spontaneous,principal_curv_mat,mole_fraction] = curvature_histogram(obj_list,type_properties,Binlimits,NumBins,1);
        end
        
        function [] = patch_geometry_plot_specific(ID,obj_list)
            obj = findobjID(obj_list,ID);
            [star_mesh,face_areas,projected_face_areas] = extract_obj_star(obj);
            surface_mesh = generate_mesh_from_objlist(obj_list);
            trisurf(surface_mesh,'FaceAlpha',0.2,'FaceColor','b');
            hold on; trisurf(star_mesh,'FaceAlpha',.2,'FaceColor','g');
            hold on; scatter3( obj.Pos(1),obj.Pos(2),obj.Pos(3),50,'r','Filled')
            hold on; quiver3( obj.Pos(1),obj.Pos(2),obj.Pos(3),obj.Normal_vect(1),obj.Normal_vect(2),obj.Normal_vect(3),4,'linewidth',2 );
            daspect([1,1,1]);hold off
        end
        
        function [obj,chemical_potential_current] = calculate_chemical_potential(obj,temperature)
            [obj] = derive_geometrical_quatities(obj);
            [star_mesh,face_areas,projected_face_areas] = extract_obj_star(obj);
            area_each_lipid = obj.avg_area_per_lipid;
            num_lipids = sum(obj.lipid_ratio_up);
%             A_tot_projected = sum(projected_face_areas)/3;
            A_tot_projected = sum(face_areas)/3;
            chemical_potential_current = temperature*log(num_lipids*area_each_lipid/A_tot_projected);
            obj.chemical_potential_initial = chemical_potential_current;            
        end
        
        function [surface_tension_current,area_face,zero_surface_tension_area] = calculate_changed_surface_tension(obj,surface_modulus)
            [star_mesh,face_areas,projected_face_areas] = extract_obj_star(obj);
            number_of_lipids_up = sum(obj.lipid_ratio_up);
            area_of_each_lipid = obj.avg_area_per_lipid;
            zero_surface_tension_area = number_of_lipids_up*area_of_each_lipid;
            area_face = sum(face_areas)/3;
            surface_tension_current = surface_modulus*number_of_lipids_up*( area_face/zero_surface_tension_area - 1 )^1;
            
        end
        
        function [obj_list] = populate_patches_with_lipid_type(obj_list,pos_on_membrane,radius,desired_mole_fraction)
            [interesting_ids] = findobj_at_pos(obj_list,pos_on_membrane,radius);
            avg_lipid_area = obj_list(1).avg_area_per_lipid;
            for patch_no = 1:length(interesting_ids)
                obj = obj_list(interesting_ids(patch_no));
                [~,face_areas,projected_face_areas] = extract_obj_star(obj);
                area_patch_current = sum(projected_face_areas)/3;
                total_num_of_lipids = ceil(area_patch_current/avg_lipid_area);
               
                lipid_type_positive_sampling = rand(total_num_of_lipids,1);
                num_lipid_type_positive = length( lipid_type_positive_sampling(lipid_type_positive_sampling<=desired_mole_fraction(1)) );                
                
                lipid_type_negative_sampling = rand(total_num_of_lipids,1);
                num_lipid_type_negative = length( lipid_type_negative_sampling(lipid_type_negative_sampling<=desired_mole_fraction(3)) );
                
                num_lipid_type_zero = total_num_of_lipids - num_lipid_type_positive - num_lipid_type_negative; 
                obj.lipid_ratio_up = [num_lipid_type_positive num_lipid_type_zero num_lipid_type_negative]; 
                
                obj_list(interesting_ids(patch_no)) = obj;
            end
        end
        
        function [obj_list] = distribute_lipids_smoothed_over_multiple_nodes(obj_list,avg_distance, desired_mole_fraction, radius_of_noise)
            % Algorithm
            % Select faces randomly
            [surface_mesh] = generate_mesh_from_objlist(obj_list);
            all_faces = surface_mesh.ConnectivityList;
            all_faces_ids = 1:size(all_faces,1);
            all_points = surface_mesh.Points;
            centroid_x = @(face_ind) mean(all_points(all_faces(face_ind,:),1));
            centroid_y = @(face_ind) mean(all_points(all_faces(face_ind,:),2));
            centroid_z = @(face_ind) mean(all_points(all_faces(face_ind,:),3));
            centroid_all_faces = [ (arrayfun(centroid_x,(1:size(all_faces,1))))',...
                (arrayfun(centroid_y,(1:size(all_faces,1))))',...
                (arrayfun(centroid_z,(1:size(all_faces,1))))'];
            [~,totalArea] = area_volume_surface(surface_mesh.Points',surface_mesh.ConnectivityList');
            number_of_points_to_choose = ceil(totalArea/(4*pi*avg_distance^2));
            randomized_index_of_faces = randperm(size(all_faces,1),number_of_points_to_choose);            
            % populate with lipids to all vertices of a face with a radius 
            selected_faces = all_faces(randomized_index_of_faces,:);
            for selected_fc_id = 1:length(randomized_index_of_faces)
                current_face = selected_faces(selected_fc_id,:);
                current_face_points = all_points(current_face,:);
                current_centroid = mean(current_face_points,1);
                distance_from_other_faces = vecnorm( centroid_all_faces - (current_centroid'*ones(1,size(all_faces,1)))',2,2);  
                effected_neighbour_face_ids = all_faces_ids( distance_from_other_faces<=radius_of_noise );
                distances_from_current_point = distance_from_other_faces( distance_from_other_faces<=radius_of_noise );
                [sorted_distance,IDs] = sort(distances_from_current_point);
                effected_neighbour_face_ids_sorted = effected_neighbour_face_ids(IDs);
                
                for effected_fc_id = 1:length(effected_neighbour_face_ids_sorted)
                    currect_face_id = effected_neighbour_face_ids_sorted(effected_fc_id);
                    current_face_of_interest = all_faces(currect_face_id,:);
                    current_face_points_of_interest = all_points(current_face_of_interest,:);
                    area_current_triangle  = area_triangle(current_face_points_of_interest);
                    
                    obj_1 = obj_list(current_face_of_interest(1));
                    obj_1_lipid_comp_up = obj_1.lipid_ratio_up;
                    obj_1_area = obj_1.Av_vertex;
                    obj_1_lipid_comp_up_fraction_old = floor((area_current_triangle/3/obj_1_area)*obj_1_lipid_comp_up);
                    total_lipids = sum(obj_1_lipid_comp_up_fraction_old);
                    
                    n_positive_lipids_sampling = rand(total_lipids,1);
                    n_positive_lipids = length(n_positive_lipids_sampling(n_positive_lipids_sampling<=desired_mole_fraction(1)));
                    
                    n_neutral_lipids_sampling = rand(total_lipids,1);
                    n_neutral_lipids = length(n_neutral_lipids_sampling(n_neutral_lipids_sampling<=desired_mole_fraction(2)));
                    
                    
                    n_negative_lipids_sampling = rand(total_lipids,1);
                    n_negative_lipids = length(n_negative_lipids_sampling(n_negative_lipids_sampling<=desired_mole_fraction(3)));
                    
                    obj_1_lipid_comp_up_fraction_new = [n_positive_lipids n_neutral_lipids n_negative_lipids];
                    obj_1_lipid_comp_up_new = obj_1_lipid_comp_up - obj_1_lipid_comp_up_fraction_old + obj_1_lipid_comp_up_fraction_new;
                    obj_1.lipid_ratio_up = obj_1_lipid_comp_up_new;
                    
                    obj_list(current_face_of_interest(1)) = obj_1;
                    
                    
                    obj_2 = obj_list(current_face_of_interest(2));
                    obj_2_lipid_comp_up = obj_2.lipid_ratio_up;
                    obj_2_area = obj_2.Av_vertex;
                    obj_2_lipid_comp_up_fraction_old = floor((area_current_triangle/3/obj_2_area)*obj_2_lipid_comp_up);
                    total_lipids = sum(obj_2_lipid_comp_up_fraction_old);
                    
                    n_positive_lipids_sampling = rand(total_lipids,1);
                    n_positive_lipids = length(n_positive_lipids_sampling(n_positive_lipids_sampling<=desired_mole_fraction(1)));
                    
                    n_neutral_lipids_sampling = rand(total_lipids,1);
                    n_neutral_lipids = length(n_neutral_lipids_sampling(n_neutral_lipids_sampling<=desired_mole_fraction(2)));
                    
                    
                    n_negative_lipids_sampling = rand(total_lipids,1);
                    n_negative_lipids = length(n_negative_lipids_sampling(n_negative_lipids_sampling<=desired_mole_fraction(3)));
                    
                    obj_2_lipid_comp_up_fraction_new = [n_positive_lipids n_neutral_lipids n_negative_lipids];
                    obj_2_lipid_comp_up_new = obj_2_lipid_comp_up - obj_2_lipid_comp_up_fraction_old + obj_2_lipid_comp_up_fraction_new;
                    obj_2.lipid_ratio_up = obj_2_lipid_comp_up_new;
                    
                    obj_list(current_face_of_interest(2)) = obj_2;
                    
                    
                    
                    obj_3 = obj_list(current_face_of_interest(3));
                    obj_3_lipid_comp_up = obj_3.lipid_ratio_up;
                    obj_3_area = obj_3.Av_vertex;
                    obj_3_lipid_comp_up_fraction_old = floor((area_current_triangle/3/obj_3_area)*obj_3_lipid_comp_up);
                    total_lipids = sum(obj_3_lipid_comp_up_fraction_old);
                    
                    n_positive_lipids_sampling = rand(total_lipids,1);
                    n_positive_lipids = length(n_positive_lipids_sampling(n_positive_lipids_sampling<=desired_mole_fraction(1)));
                    
                    n_neutral_lipids_sampling = rand(total_lipids,1);
                    n_neutral_lipids = length(n_neutral_lipids_sampling(n_neutral_lipids_sampling<=desired_mole_fraction(2)));
                    
                    
                    n_negative_lipids_sampling = rand(total_lipids,1);
                    n_negative_lipids = length(n_negative_lipids_sampling(n_negative_lipids_sampling<=desired_mole_fraction(3)));
                    
                    obj_3_lipid_comp_up_fraction_new = [n_positive_lipids n_neutral_lipids n_negative_lipids];
                    obj_3_lipid_comp_up_new = obj_3_lipid_comp_up - obj_3_lipid_comp_up_fraction_old + obj_3_lipid_comp_up_fraction_new;
                    obj_3.lipid_ratio_up = obj_3_lipid_comp_up_new;
                    
                    obj_list(current_face_of_interest(3)) = obj_3;
                    
                    
                end
                
                
                
                
%                 area_current_triangle  = area_triangle(current_face_points);                
%                 obj_1 = obj_list(current_face(1));
%                 obj_1_lipid_comp_up = obj_1.lipid_ratio_up;
%                 obj_1_area = obj_1.Av_vertex;
%                 obj_1_lipid_comp_up_fraction_old = area_current_triangle/3/obj_1_area;
                
                
            end
            
        end
        
        function [adjacent_areas_to_edges,total_area_patch] = calculate_adjacent_areas_to_neighbours(obj)
            [star_mesh,face_areas] = extract_obj_only_star(obj);
            total_area_patch = sum(face_areas);
            faces_all = star_mesh.ConnectivityList;
            num_connectors = size(star_mesh.Points,1)-1;
            connector_ids = 2:size(star_mesh.Points,1);
            adjacent_areas_to_edges = NaN(num_connectors,1);
            all_points = star_mesh.Points;
            for connector_ind = 1:num_connectors
                current_connector = [1 connector_ids(connector_ind)];
                connector_in_face_indicator_mat = double(ismember(faces_all,current_connector));
                connector_in_face_indicator_logical = sum(connector_in_face_indicator_mat,2);
                connector_in_face_indicator_logical = connector_in_face_indicator_logical==2;
                corresponding_faces = faces_all(connector_in_face_indicator_logical,:);
                face1 = corresponding_faces(1,:);
                face1_points = all_points(face1,:);
                face1_area = area_triangle(face1_points);
                
                face2 = corresponding_faces(2,:);
                face2_points = all_points(face2,:);
                face2_area = area_triangle(face2_points);
                adjacent_areas_to_edges(connector_ind) = face1_area + face2_area;                
            end
        end
        
        function [energy_stored] = calculate_enthalpic_energy_lipids(obj,lipid_curvatures,enthalpy_per_lipid,instantaneous_curvature)
            lipid_num_fractions = obj.lipid_ratio_up;
%             lipid_mole_fractions = lipid_num_fractions./sum(lipid_num_fractions);
            penalty_curvature_array = 0.5*lipid_num_fractions.*(lipid_curvatures - instantaneous_curvature).^2;
            energy_stored = sum( enthalpy_per_lipid.*penalty_curvature_array );
            
            
        end
        
        function [] = protein_bind_or_not(obj,protein_info_obj,surface_mesh,protein_presence_counter)
            protein_size = protein_info_obj.size_radius;
            
            
        end
        
        
    end
    
end

