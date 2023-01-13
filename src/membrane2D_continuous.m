classdef membrane2D_continuous
    
    properties
        % basic
        ID
        Pos
        
        % Geometry
        Normal_vect
        c1 %principle curvature1
        c2 %principle curvature2
        principal_directions
        Av_vertex
        
        % constants
        faces
        neighbours
        avg_area_per_lipid
        
        % topology
        neighbours_coor
        neighbours_normal_old;
        
        % Composition
        n_particles_up
        neighbours_lipid_compositions_up
    end
    
    methods
        %% Initialize
        function obj = membrane2D_continuous(ID,Pos,neighboursID,faces_current,neighbours_coor_current)
            if nargin ~= 0
                obj.ID = ID;
                obj.Pos = Pos;
                obj.neighbours = neighboursID;
                obj.faces = faces_current;
                obj.neighbours_coor = neighbours_coor_current;
            end
        end
        
        function [obj_list] = initilize_obj_list_container_from_mesh(obj,membrane_mesh)
            obj_list(1) = obj;
            all_pos = membrane_mesh.Points;
            all_IDs = (1:size(all_pos,1))';
            connectivity_mat = membrane_mesh.ConnectivityList;
            %             theta_list = 0:2*pi/num_objs:2*pi-2*pi/num_objs;
            for obj_no = 1:size(all_pos,1)
                
                row_mask_for_neighbour = ...
                    [ connectivity_mat(:,1)==obj_no,connectivity_mat(:,2)==obj_no,connectivity_mat(:,3)==obj_no ];
                row_mask_for_neighbour = sum(row_mask_for_neighbour,2);
                row_mask_for_neighbour = row_mask_for_neighbour==1;
                faces_current = connectivity_mat(row_mask_for_neighbour,:);
                neighboursID = unique(faces_current(:));
                neighboursID = neighboursID(neighboursID~=obj_no);
                neighbours_coor_current = all_pos(neighboursID,:);
                
                obj = membrane2D_continuous(all_IDs(obj_no),all_pos(obj_no,:),neighboursID,faces_current,neighbours_coor_current);
                obj_list(obj_no) = obj;
            end
            obj_list = obj_list';
        end
        
        function [updated_list] = populate_with_lipids_at_equillibrium(obj_list,lipid_head_area,spont_curv_lipids)
            [updated_list] = derive_geometrical_quatities_all(obj_list);
            parfor obj_no = 1:size(updated_list,1)
                current_obj = updated_list(obj_no,1);
                area_surface = current_obj.Av_vertex;
                N_lipids = floor(area_surface/lipid_head_area);
                mole_fraction_type_a = ((current_obj.c1 + current_obj.c2)/2 - spont_curv_lipids(2))/(spont_curv_lipids(1)-spont_curv_lipids(2));
                mole_fraction_type_b = 1-mole_fraction_type_a;
                current_obj.n_particles_up = [ floor(mole_fraction_type_a*N_lipids) floor(mole_fraction_type_b*N_lipids) ];
                current_obj.avg_area_per_lipid = lipid_head_area;
                updated_list(obj_no,1) = current_obj;
            end
        end
        
        function [updated_list] = populate_with_lipids_uniform(obj_list,lipid_head_area,spont_curv_lipids)
            [updated_list] = derive_geometrical_quatities_all(obj_list);
            [mean_curvature,~] = extract_curvature_list(updated_list);
            mean_curvature_all = mean(mean_curvature,1);
            
            parfor obj_no = 1:size(updated_list,1)
                current_obj = updated_list(obj_no,1);
                area_surface = current_obj.Av_vertex;
                N_lipids = floor(area_surface/lipid_head_area);
                mole_fraction_type_a = (mean_curvature_all - spont_curv_lipids(2))/(spont_curv_lipids(1)-spont_curv_lipids(2));
                mole_fraction_type_b = 1-mole_fraction_type_a;
                current_obj.n_particles_up = [ floor(mole_fraction_type_a*N_lipids) floor(mole_fraction_type_b*N_lipids) ];
                current_obj.avg_area_per_lipid = lipid_head_area;
                updated_list(obj_no,1) = current_obj;
            end
        end
        
        
        
        %% Topology
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
        
        %% Energy
        function [free_energy_at_vertex,obj_moved] = calculate_free_energy_at_vertex(obj, obj_pos,obj_lipid_comp,all_pos_list,spont_curv_lipids, kappa, surface_modulous,kappa_coefficient,vanderwaals, kB, temperature)
            obj_moved = obj;
            obj_moved.Pos = obj_pos;
            obj_moved.n_particles_up = obj_lipid_comp;
            obj_moved.neighbours_coor = all_pos_list(obj.neighbours,:);
            [c1_current,c2_current,~,~,Av,obj_moved,~] = determine_principal_curvatures(obj_moved);
            nA = obj_moved.n_particles_up(1);
            nB = obj_moved.n_particles_up(2);
            A0 = (nA+nB)*obj_moved.avg_area_per_lipid;
            per_lipid_area = obj_moved.avg_area_per_lipid;
            
            h_spont = @(nA,nB) nA/(nA + nB)*spont_curv_lipids(1) + nB/(nA + nB)*spont_curv_lipids(2);
            free_energy_hamiltonian = @(c1,c2,nA,nB) 0.5*kappa*(nA + nB)*((c1+c2)/2 - 1*h_spont(nA,nB))^2 + ...
                0.5*kappa*kappa_coefficient*nA*(1*h_spont(nA,nB) - 1*spont_curv_lipids(1))^2 + 0.5*kappa*kappa_coefficient*nB*(1*h_spont(nA,nB) - 1*spont_curv_lipids(2))^2 + ...
                 +vanderwaals*nA*nB/(nA+nB) + ...
                 0.5*surface_modulous*((Av - A0)^2)/A0 + ...
                -temperature*kB*( nA*log(Av/(nA*per_lipid_area)) + nB*log(Av/(nB*per_lipid_area)) + nA+nB);
            
            free_energy_at_vertex =  free_energy_hamiltonian(c1_current,c2_current,nA,nB);
            
            
        end
        
        function [free_energy_of_patch,center_obj] = ...
                calculate_free_energy_patch(center_obj,center_obj_pos,obj_list,all_pos_list,spont_curv_lipids, kappa, surface_modulous,kappa_coefficient,vanderwaals, kB, temperature)
            % update the pos_list with the current center vertex
            all_pos_list(center_obj.ID,:) =  center_obj_pos; %% this is neccesary to ensure that the neighbours neighbour vertex include the updated coordinate of center object
            
            % Calculate free energy at center with the new position
            [free_energy_center,center_obj] = ...
                calculate_free_energy_at_vertex(center_obj, center_obj_pos,center_obj.n_particles_up,all_pos_list,spont_curv_lipids, kappa, surface_modulous,kappa_coefficient,vanderwaals, kB, temperature);
            
            % Update the neighbour_objs neighbour_coors property with
            % the new position of the center obj and calculate energy
            neighbours_ID = center_obj.neighbours;
            neighbours_obj = obj_list(neighbours_ID,:);
            free_energy_neigh = NaN(length(neighbours_ID),1);
            Av_neighbours = NaN(length(neighbours_ID),1);
            %                 center_obj_ID = center_obj.ID;
            for neigh_obj = 1:length(neighbours_ID)
                current_neighbour_obj = neighbours_obj(neigh_obj);
                neigh_obj_pos = all_pos_list(current_neighbour_obj.ID,:);
                %                     current_neighbour_obj.neighbours_coor = all_pos_list( current_neighbour_obj.neighbours,: );
                %                     current_neighbour_obj.neighbours_coor( current_neighbour_obj.neighbours==center_obj_ID,: ) = center_obj_pos;
                [free_energy_neigh(neigh_obj),neighbours_obj(neigh_obj)] = ...
                    calculate_free_energy_at_vertex(current_neighbour_obj, neigh_obj_pos,current_neighbour_obj.n_particles_up,all_pos_list,spont_curv_lipids, kappa, surface_modulous,kappa_coefficient,vanderwaals, kB, temperature);
                Av_neighbours(neigh_obj,1)= neighbours_obj(neigh_obj).Av_vertex;
            end
            % do the sum free energy for the patch based on area
            % contribution
            free_energy_of_patch = free_energy_center + sum(free_energy_neigh);
            
        end
        
        function [penalty_energy] = calculate_composition_gradient_penalty(obj1,obj2,penalty_factor)
            lipid_comp1 = obj1.n_particles_up;
            lipid_comp2 = obj2.n_particles_up;
            
            mole_frac1 = lipid_comp1(:,1)./sum(lipid_comp1);
            mole_frac2 = lipid_comp2(:,1)./sum(lipid_comp2);
            
            total_lipid_type_based = lipid_comp1+lipid_comp2;
            penalty_energy = sum(penalty_factor*total_lipid_type_based.*(mole_frac1 - mole_frac2).^2);
            
        end
        
        function [penalty_energy,penalty_at_each_edge] = calculate_composition_gradient_at_all_edges(obj,lipid_comp_current,all_lipid_comp_list,penalty_factor)
            all_lipid_comp_list(obj.ID,:) = lipid_comp_current;
            [adjacent_areas_to_edges,total_area_patch] = calculate_adjacent_areas_to_neighbours(obj);
            area_patch = total_area_patch/3;
            area_for_edge_in_patch = 1/3*adjacent_areas_to_edges/2;
            weights_for_penalty = area_for_edge_in_patch/area_patch;
            
            
            mole_lipid_comp = lipid_comp_current./sum(lipid_comp_current,2);
            num_neighbours = length(obj.neighbours);
            
            neighbours_lipid_comp = all_lipid_comp_list(obj.neighbours,:);
            mole_neighbour_lipid_comp = [neighbours_lipid_comp(:,1)./sum(neighbours_lipid_comp,2) , neighbours_lipid_comp(:,2)./sum(neighbours_lipid_comp,2)  ];
            
            difference_conc_mat = (mole_lipid_comp'*ones( 1,num_neighbours ))'-mole_neighbour_lipid_comp;
            
            % important note behind the concept: since the molefraction of
            % lipid1 is constrained with mole fraction of lipid2, the
            % inhomogenity at boundary is described by the difference in
            % one type of lipid
            
            grad_type1_lipid_mole = difference_conc_mat(:,1).^2;
            penalty_at_each_edge = penalty_factor*sum(lipid_comp_current)*weights_for_penalty.*grad_type1_lipid_mole;
            penalty_energy = sum(penalty_at_each_edge);
            
            
        end
        
        function [penalty_energy,penalty_of_objs] = ...
                calculate_composition_gradient_penalty_ofPair(obj1,lipid_comp_current1,obj2,lipid_comp_current2,all_lipid_comp_list,penalty_factor)
            all_lipid_comp_list(obj1.ID,:) = lipid_comp_current1;
            all_lipid_comp_list(obj2.ID,:) = lipid_comp_current2;
            
            obj1.neighbours_lipid_compositions_up = all_lipid_comp_list(obj1.neighbours,:);
            obj2.neighbours_lipid_compositions_up = all_lipid_comp_list(obj2.neighbours,:);
            
%             [~,grad_vect_along_surf_obj1] = calculate_gradient_least_square_method(obj1);
%             [~,grad_vect_along_surf_obj2] = calculate_gradient_least_square_method(obj2);
            [grad_vect_along_surf_obj1] = calculate_gradient_green_gauss_method(obj1);
            [grad_vect_along_surf_obj2] = calculate_gradient_green_gauss_method(obj2);
            
            grad_mag_obj1 = vecnorm(grad_vect_along_surf_obj1,2,1);
            grad_mag_obj2 = vecnorm(grad_vect_along_surf_obj2,2,1);
            
            N_lipid_obj1 = sum(lipid_comp_current1);
            N_lipid_obj2 = sum(lipid_comp_current2);
            
            penalty_obj1 = penalty_factor*N_lipid_obj1*grad_mag_obj1^2;
            penalty_obj2 = penalty_factor*N_lipid_obj2*grad_mag_obj2^2;
            
            penalty_energy = penalty_obj1 + penalty_obj2;
            penalty_of_objs = [ penalty_obj1 penalty_obj2 ];
        end
        
        %% MC moves
        
        function [obj_list] = orthogonal_move(obj_list,step_size,spont_curv_lipids, kappa, surface_modulous,kappa_coefficient,vanderwaals, kB, temperature)
            % randomly sample points
            total_objs = size(obj_list,1);
            randomized_order = randperm(total_objs);
            all_pos_list = extract_coors_points_ID_specific(obj_list);
            
            % update pos_list based on metropolis
            for obj_seq = 1:length(randomized_order)
                current_ID = randomized_order(obj_seq);
                current_obj = obj_list(current_ID,1);
                current_obj_pos = all_pos_list(current_ID,:);
                
                % Calculate patch energy before move
                [initial_free_energy_of_patch,current_obj] = ...
                    calculate_free_energy_patch(current_obj,current_obj_pos,obj_list,all_pos_list,spont_curv_lipids, kappa, surface_modulous,kappa_coefficient,vanderwaals, kB, temperature);
                normal_direction_current = current_obj.Normal_vect;
                unit_normal_current = normal_direction_current./vecnorm(normal_direction_current,2,2);
                
                % make the move using jumping probability
                kick_displacement = step_size*randn(1,1);
                updated_obj_pos = current_obj_pos + kick_displacement*unit_normal_current;
                [final_free_energy_of_patch,updated_obj] = ...
                    calculate_free_energy_patch(current_obj,updated_obj_pos,obj_list,all_pos_list,spont_curv_lipids, kappa, surface_modulous,kappa_coefficient,vanderwaals, kB, temperature);
                
                % calculate del_free_energy and apply metropolis
                del_free_energy = final_free_energy_of_patch - initial_free_energy_of_patch;
                decision = min( [1,exp( -del_free_energy/(kB*temperature) )] );
                
                uniform_sampler = rand(1,1);
                if uniform_sampler <= decision
                    all_pos_list(current_ID,:) = updated_obj_pos;
                    fprintf('FUNCTION: orthogonal_move; UPDATED since U(%d)<= decision(%d)\n',uniform_sampler,decision);
                else
                    fprintf('FUNCTION: orthogonal_move; NOT UPDATED since U(%d)> decision(%d)\n',uniform_sampler,decision);
                end
            end
            
            % update obj_list Pos and neighbour_coor property
            parfor obj_no = 1:size(obj_list,1)
                current_obj = obj_list(obj_no,1);
                updated_obj = current_obj;
                updated_obj.Pos = all_pos_list(current_obj.ID,:);
                updated_obj.neighbours_coor = all_pos_list(current_obj.neighbours,:);
                [~,~,~,~,~,updated_obj,~] = determine_principal_curvatures(updated_obj);
                obj_list(obj_no,1) = updated_obj;
                fprintf('FUNCTION: orthogonal_move; LOADING NEW PROPERTIES ID(%d)\n',current_obj.ID);
            end
            
        end
        
        function [obj_list] = exchange_move(obj_list,particle_packet_size,spont_curv_lipids, kappa, surface_modulous,penalty_factor,kappa_coefficient,vanderwaals, kB, temperature)
            % choose random objects
            total_objs = size(obj_list,1);
            randomized_order_first = randperm(total_objs);
            % choose neighbour from each obj
            ID_pairs = [randomized_order_first',NaN(length(randomized_order_first),1)];
            for obj_seq = 1:length(randomized_order_first)
                current_obj = obj_list(randomized_order_first(obj_seq),1);
                rand_neigh_ID = current_obj.neighbours( randperm( length(current_obj.neighbours),1 ) );
                ID_pairs(obj_seq,2) = rand_neigh_ID;
            end
            
            [all_lipid_comp] = extract_lipids_composition_ID_specific(obj_list);
            all_pos_list = extract_coors_points_ID_specific(obj_list);
            
            for pair_no = 1:size(ID_pairs,1)
                current_obj_1 = obj_list(ID_pairs(pair_no,1),1);
                current_obj_2 = obj_list(ID_pairs(pair_no,2),1);
                
                current_obj1_particles = all_lipid_comp(current_obj_1.ID,:);
                current_obj2_particles = all_lipid_comp(current_obj_2.ID,:);
                % calculate initial free energy of each pair
                [initial_free_energy_obj1,current_obj_1] = ...
                    calculate_free_energy_at_vertex(current_obj_1, current_obj_1.Pos,current_obj1_particles,all_pos_list,spont_curv_lipids, kappa, surface_modulous,kappa_coefficient,vanderwaals, kB, temperature);
                [initial_free_energy_obj2,current_obj_2] = ...
                    calculate_free_energy_at_vertex(current_obj_2, current_obj_2.Pos,current_obj2_particles,all_pos_list,spont_curv_lipids, kappa, surface_modulous,kappa_coefficient,vanderwaals, kB, temperature);
                
                %                 [initial_penalty_energy] = calculate_composition_gradient_penalty(current_obj_1,current_obj_2,penalty_factor);
                [initial_penalty_energy,initial_penalty_of_objs] = ...
                    calculate_composition_gradient_penalty_ofPair(current_obj_1,current_obj1_particles,current_obj_2,current_obj2_particles,all_lipid_comp,penalty_factor);
                
                
                % calculate final free energy of each pair after exchange
                packet_to_transfer_from1_to_2 = round( particle_packet_size*[ randn(1,1),randn(1,1) ] );
                updated_obj1_particles = current_obj1_particles - packet_to_transfer_from1_to_2;
                updated_obj2_particles = current_obj2_particles + packet_to_transfer_from1_to_2;
                
                particle_num_negative_check = [ updated_obj1_particles updated_obj2_particles ];
                is_there_negative = min(particle_num_negative_check)<0;
                
                [final_free_energy_obj1,current_obj_1] = ...
                    calculate_free_energy_at_vertex(current_obj_1, current_obj_1.Pos,updated_obj1_particles,all_pos_list,spont_curv_lipids, kappa, surface_modulous,kappa_coefficient,vanderwaals, kB, temperature);
                [final_free_energy_obj2,current_obj_2] = ...
                    calculate_free_energy_at_vertex(current_obj_2, current_obj_2.Pos,updated_obj2_particles,all_pos_list,spont_curv_lipids, kappa, surface_modulous,kappa_coefficient,vanderwaals, kB, temperature);
                
                %                 [final_penalty_energy] = calculate_composition_gradient_penalty(current_obj_1,current_obj_2,penalty_factor);
                [final_penalty_energy,final_penalty_of_objs] = ...
                    calculate_composition_gradient_penalty_ofPair(current_obj_1,updated_obj1_particles,current_obj_2,updated_obj2_particles,all_lipid_comp,penalty_factor);
                
                
                % apply metropolis
                del_free_energy = (final_free_energy_obj1 + final_free_energy_obj2 + final_penalty_energy) - (initial_free_energy_obj1 + initial_free_energy_obj2+initial_penalty_energy);
                
                decision = min( [1, exp( -del_free_energy/(kB*temperature) ) ] );
                uniform_sampler = rand(1,1);
                if uniform_sampler <= decision && not(is_there_negative)
                    all_lipid_comp(current_obj_1.ID,:) = updated_obj1_particles;
                    all_lipid_comp(current_obj_2.ID,:) = updated_obj2_particles;
                    fprintf('FUNCTION: exchange_move; UPDATED since U(%d)<= decision(%d)\n',uniform_sampler,decision);
                else
                    fprintf('FUNCTION: exchange_move; NOT UPDATED since U(%d)> decision(%d) or lipid_numbers is negative (is_there_negative %d)\n',uniform_sampler,decision,is_there_negative);
                end
                
                
                
            end
            
            % load lipid_comp_list into obj_list
            parfor obj_no = 1:size(obj_list,1)
                current_obj = obj_list(obj_no,1);
                updated_obj = current_obj;
                updated_obj.n_particles_up = all_lipid_comp(current_obj.ID,:);
                %                 updated_obj.neighbours_coor = all_pos_list(current_obj.neighbours,:);
                [~,~,~,~,~,updated_obj,~] = determine_principal_curvatures(updated_obj);
                obj_list(obj_no,1) = updated_obj;
                fprintf('FUNCTION: exchange_move; LOADING NEW PROPERTIES ID(%d)\n',current_obj.ID);
            end
            
        end
        
        
        
        
        %% Accesory
        function [IDs] = extract_IDs(obj_list)
            extract_ID_func = @(id) obj_list(id).ID;
            IDs = arrayfun(extract_ID_func,1:size(obj_list,1));
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
        
        function [normal_vects_all] = extract_Normal_vects(obj_list)
            extract_Normal_func1 = @(id) obj_list(id).Normal_vect(1);
            extract_Normal_func2 = @(id) obj_list(id).Normal_vect(2);
            extract_Normal_func3 = @(id) obj_list(id).Normal_vect(3);
            normal_vects_all = [(arrayfun(extract_Normal_func1,1:size(obj_list,1)))',...
                (arrayfun(extract_Normal_func2,1:size(obj_list,1)))',...
                (arrayfun(extract_Normal_func3,1:size(obj_list,1)))'];
        end
        
        function [mean_curvature,principle_curvatures] = extract_curvature_list(obj_list)
            principle_curvatures = NaN(size(obj_list,1),2);
            extract_c1_func = @(id) obj_list(id).c1;
            principle_curvatures(:,1) = arrayfun(extract_c1_func,1:size(obj_list,1));
            extract_c2_func = @(id) obj_list(id).c2;
            principle_curvatures(:,2) = arrayfun(extract_c2_func,1:size(obj_list,1));
            mean_curvature = mean(principle_curvatures,2);
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
        
        
        function [neighbours_num_all] = extract_num_neighbours_ID_specific(obj_list)
            neighbours_num_all = NaN(size(obj_list,1),1);
            ID_extract = @(id)obj_list(id).ID;
            ID_extract_all = arrayfun(ID_extract,1:size(obj_list,1));
            for obj_id = 1:size(obj_list,1)
                ID_current = ID_extract_all(obj_id);
                neighbours_num_all(ID_current,:) = length(obj_list(obj_id).neighbours);
            end
            
        end
        
        function [lipid_comp] = extract_lipids_composition_ID_specific(obj_list)
            lipid_comp = NaN(size(obj_list,1),2);
            ID_extract = @(id)obj_list(id).ID;
            ID_extract_all = arrayfun(ID_extract,1:size(obj_list,1));
            for obj_id = 1:size(obj_list,1)
                ID_current = ID_extract_all(obj_id);
                lipid_comp(ID_current,:) = obj_list(obj_id).n_particles_up;
            end
            
        end
        
        function [triangle_color_curvature] = color_code_patch_wise_curvature(obj_list)
            
            fprintf('>>> Function: color_code_curvature_profile_patch_wise; Generating surface mesh for mesh coloring\n');
            [surface_mesh] = generate_mesh_from_objlist(obj_list);
            triangles = surface_mesh.ConnectivityList;
            
            triangle_color_curvature = NaN(size(triangles,1),1);
            
            [mean_curvature,~] = extract_curvature_list(obj_list);
            for tri_ind = 1:size(triangles,1)
                current_triangle = triangles(tri_ind,:);
                triangle_color_curvature(tri_ind,:) = mean(mean_curvature(current_triangle));
            end
            fprintf('>>> Function: color_code_lipid_density_profile_patch_wise; Generating mesh color Done \n\n');
            
            
        end
        
        function [triangle_color_mole_fraction] = color_code_patch_wise_composition(obj_list)
            %             rgb_color_mat = NaN(size(obj_list,1),3);
            fprintf('>>> Function: color_code_lipid_density_profile_patch_wise; RGB Color coding patches based on mole fraction\n');
            [lipid_comp] = extract_lipids_composition_ID_specific(obj_list);
            %             total_num_lipid_array = sum(lipid_comp,2);
            fprintf('>>> Function: color_code_lipid_density_profile_patch_wise; Calculating the mole fraction\n');
            %             lipid_comp_mole_fraction = [ lipid_comp(:,1)./total_num_lipid_array lipid_comp(:,2)./total_num_lipid_array lipid_comp(:,3)./total_num_lipid_array  ];
            fprintf('>>> Function: color_code_lipid_density_profile_patch_wise; Generating surface mesh for mesh coloring\n');
            [surface_mesh] = generate_mesh_from_objlist(obj_list);
            triangles = surface_mesh.ConnectivityList;
            Points_all = surface_mesh.Points;
            fprintf('>>> Function: color_code_lipid_density_profile_patch_wise; Extracting number of neighbours for weighting \n');
            %             [neighbours_num_all] = extract_num_neighbours_ID_specific(obj_list);
            triangle_color_mole_fraction = NaN(size(triangles,1),2);
            %             triangle_color_curvature = NaN(size(triangles,1),1);
            fprintf('>>> Function: color_code_lipid_density_profile_patch_wise; Generating mesh color \n');
            all_face_areas = @(ind) obj_list(ind).Av_vertex;
            all_areas = arrayfun(all_face_areas,[1:size(obj_list)]);
            
            
            %             [mean_curvature,~] = extract_curvature_list(obj_list);
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
                
                
                %                 weights = neighbours_num_all(current_triangle);
                %                 weighted_lipid_comp = [ lipid_comp_of_vertices(:,1).*weights lipid_comp_of_vertices(:,2).*weights lipid_comp_of_vertices(:,3).*weights ];
                %                 total_lipids_in_current_triangle = sum(weighted_lipid_comp(:));
                %                 lipid_comp_type_wise = sum(weighted_lipid_comp,1);
                %                 lipid_comp_type_wise_mole_fraction = lipid_comp_type_wise/total_lipids_in_current_triangle;
                triangle_color_mole_fraction(tri_ind,:) = sum_triangle_lipid_mole_fraction;
                %                 triangle_color_curvature(tri_ind,:) = mean(mean_curvature(current_triangle));
            end
            fprintf('>>> Function: color_code_lipid_density_profile_patch_wise; Generating mesh color Done \n\n');
            
            
        end
        
        function [coors,principle_curvatures,hmean,lipid_comp,mole_lipid_comp,order_parameter]  = extract_patch_properties(obj_list,spont_curv_lipids)
            [coors] = extract_coors_points_ID_specific(obj_list);
            [hmean,principle_curvatures] = extract_curvature_list(obj_list);
            [lipid_comp] = extract_lipids_composition_ID_specific(obj_list);
            mole_lipid_comp = [ lipid_comp(:,1)./sum(lipid_comp,2),lipid_comp(:,2)./sum(lipid_comp,2) ];
            order_parameter =  sum(mole_lipid_comp.*(spont_curv_lipids'*ones(1,size(obj_list,1)))' ,2 );
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
        
        function [grad_vect,grad_vect_along_surf] = calculate_gradient_least_square_method(obj)
            %% reference https://www.youtube.com/watch?v=7ymFkxx2R_k&t=708s
            center_coor = obj.Pos;
            neighbour_coor_current = obj.neighbours_coor;
            
            distance_mat = neighbour_coor_current - ( center_coor'*ones(1,size(neighbour_coor_current,1)) )';
            lipid_comp_neighbour = obj.neighbours_lipid_compositions_up;
            
            mole_lipid_comp_neighbour_type1 = lipid_comp_neighbour(:,1)./sum(lipid_comp_neighbour,2);
            mole_lipid_comp_center_type1 = obj.n_particles_up(1)/sum(obj.n_particles_up(:));
            
            difference_conc = abs(mole_lipid_comp_neighbour_type1 - mole_lipid_comp_center_type1);
            weight_mat = diag( 1./vecnorm(distance_mat,2,2) );
            
            G_mat = distance_mat'*weight_mat'*weight_mat*distance_mat;
            
            grad_vect = inv(G_mat)*distance_mat'*weight_mat'*weight_mat*difference_conc;
            
            grad_vect_mag = vecnorm(grad_vect,2,1);
            direction_on_surface = grad_vect - dot(obj.Normal_vect(:),grad_vect(:))*(obj.Normal_vect)';
            unit_direction_vect = direction_on_surface./vecnorm(direction_on_surface,2,1);
            
            grad_vect_along_surf = grad_vect_mag*unit_direction_vect;
            
            
            
        end
        
        function [grad_vect,grad_vect_along_surf] = calculate_gradient_least_square_method_modified(obj)
            %% reference https://www.youtube.com/watch?v=7ymFkxx2R_k&t=708s
            center_coor = obj.Pos;
            neighbour_coor_current = obj.neighbours_coor;
            
            distance_mat = neighbour_coor_current - ( center_coor'*ones(1,size(neighbour_coor_current,1)) )';
            lipid_comp_neighbour = obj.neighbours_lipid_compositions_up;
            
            mole_lipid_comp_neighbour_type2 = lipid_comp_neighbour(:,2)./sum(lipid_comp_neighbour,2);
            mole_lipid_comp_center_type2 = obj.n_particles_up(2)/sum(obj.n_particles_up(:));
            
            difference_conc = (mole_lipid_comp_neighbour_type2 - mole_lipid_comp_center_type2);
            weight_mat = diag( 1./vecnorm(distance_mat,2,2) );
            
            G_mat = distance_mat'*weight_mat'*weight_mat*distance_mat;
            
            grad_vect = inv(G_mat)*distance_mat'*weight_mat'*weight_mat*difference_conc;
            
            grad_vect_mag = vecnorm(grad_vect,2,1);
            direction_on_surface = grad_vect - dot(obj.Normal_vect(:),grad_vect(:))*(obj.Normal_vect)';
            unit_direction_vect = direction_on_surface./vecnorm(direction_on_surface,2,1);
            
            grad_vect_along_surf = grad_vect_mag*unit_direction_vect;
            
            
            
        end
        
        
        function [grad_vect_along_surf] = calculate_gradient_green_gauss_method(obj)
            %% reference https://www.youtube.com/watch?v=7ymFkxx2R_k&t=708s
            center_coor = obj.Pos;
            neighbour_coor_current = obj.neighbours_coor;
            
            distance_mat = neighbour_coor_current - ( center_coor'*ones(1,size(neighbour_coor_current,1)) )';
            lipid_comp_neighbour = obj.neighbours_lipid_compositions_up;
            
%             mole_lipid_comp_neighbour_type1 = lipid_comp_neighbour(:,1)./sum(lipid_comp_neighbour,2);
%             mole_lipid_comp_center_type1 = obj.n_particles_up(1)/sum(obj.n_particles_up(:));
%             
%             difference_conc = abs(mole_lipid_comp_neighbour_type1 - mole_lipid_comp_center_type1);
%             weight_mat = diag( 1./vecnorm(distance_mat,2,2) );
%             
%             G_mat = distance_mat'*weight_mat'*weight_mat*distance_mat;
%             
%             grad_vect = inv(G_mat)*distance_mat'*weight_mat'*weight_mat*difference_conc;
%             
%             grad_vect_mag = vecnorm(grad_vect,2,1);
%             direction_on_surface = grad_vect - dot(obj.Normal_vect(:),grad_vect(:))*(obj.Normal_vect)';
%             unit_direction_vect = direction_on_surface./vecnorm(direction_on_surface,2,1);
            
%             grad_vect_along_surf = grad_vect_mag*unit_direction_vect;        
            
            
            % Calculate the edge ID pairs
            current_faces = obj.faces;
            edge_pairs = current_faces;
            edge_pairs(current_faces==obj.ID) = 0;
            edge_pairs = sort(edge_pairs,2);
            edge_pairs = edge_pairs(:,2:3);
            
            %Calculate the normal face vectors of the edges
            edgeVertex1ID = edge_pairs(:,1);
            [~,edgeVertex1IDinNeigList] = ismember(edgeVertex1ID,obj.neighbours); 
            edgeVertex1Pos = neighbour_coor_current(edgeVertex1IDinNeigList,:);
            
            edgeVertex2ID = edge_pairs(:,2);
            [~,edgeVertex2IDinNeigList] = ismember(edgeVertex2ID,obj.neighbours); 
            edgeVertex2Pos = neighbour_coor_current(edgeVertex2IDinNeigList,:);
            
            edge_12_vects = edgeVertex2Pos - edgeVertex1Pos;
            unit_edge_12_vects = ...
                [edge_12_vects(:,1)./vecnorm(edge_12_vects,2,2),...
                edge_12_vects(:,2)./vecnorm(edge_12_vects,2,2),...
                edge_12_vects(:,3)./vecnorm(edge_12_vects,2,2)]; %% The set of vectors from the first vertex of each edge to the 2nd vertex
            
            edge_10_vects = (obj.Pos'*ones(1,length(obj.neighbours)))' - edgeVertex1Pos; %% The set of vectors from the first vertex of each edge to the center vertex
            
            point_on_edge_from_center_as_normal = edgeVertex1Pos + (( dot(edge_10_vects,unit_edge_12_vects,2) )*ones(1,3)).*unit_edge_12_vects;
            normal_direction_vector_fromcenter = point_on_edge_from_center_as_normal - (obj.Pos'*ones(1,length(obj.neighbours)))';
            
            
            normal_vect_repeated = ((obj.Normal_vect)'*ones(1,length(obj.neighbours)))';
            normal_direction_vector_fromcenter_surface_plane = normal_direction_vector_fromcenter ...
                                              - (dot(normal_direction_vector_fromcenter,normal_vect_repeated,2)*ones(1,3)).*normal_vect_repeated;
                                          
            magnitude = vecnorm(normal_direction_vector_fromcenter_surface_plane,2,2);                              
            unit_normal_toEdge = [ normal_direction_vector_fromcenter_surface_plane(:,1)./magnitude, ...
                normal_direction_vector_fromcenter_surface_plane(:,2)./magnitude, ...
                normal_direction_vector_fromcenter_surface_plane(:,3)./magnitude];
          
           % Calculate the scaler at the edge centers by linear interpolation
           scaler1Pos = lipid_comp_neighbour(edgeVertex1IDinNeigList,:);
           scaler2Pos = lipid_comp_neighbour(edgeVertex2IDinNeigList,:);
           scaler_at_edge_center = (scaler1Pos + scaler2Pos)/2;
           scaler_at_edge_center_type2 = scaler_at_edge_center(:,2)./sum(scaler_at_edge_center,2);
           
           
           %Calculate the scaler flux
           edge_length = vecnorm(edge_12_vects,2,2);
           flux_at_edge_mag = edge_length.*scaler_at_edge_center_type2;
           flux_vect = (flux_at_edge_mag*[1,1,1]).*unit_normal_toEdge;
           net_flux = sum(flux_vect,1);
           grad_vect_along_surf = net_flux/obj.Av_vertex; 
        end
        
        
        function [scaler_interpolated,interpolated_mesh,surface_mesh] = interpolate_vertex_scalers_to_mesh(obj_list,scaler_field)
            [surface_mesh] = generate_mesh_from_objlist(obj_list);
            circum_centers = circumcenter(surface_mesh);
            triangles = surface_mesh.ConnectivityList;
            all_points = surface_mesh.Points;
            
            % know the number of edges
            triangles_sorted = sort(triangles,2);
            edges_real_mesh = [triangles_sorted(:,1),triangles_sorted(:,2); ...
                triangles_sorted(:,2),triangles_sorted(:,3); ...
                triangles_sorted(:,1),triangles_sorted(:,3)];
            edges_real_mesh_unique = unique(edges_real_mesh,'rows');
            total_num_points_new_mesh = size(edges_real_mesh_unique,1) + size(circum_centers,1) + size(all_points,1);
            
            
            % loop
            interpolated_connectivity_list = NaN( 6*size(triangles,1),3);
            scaler_interpolated = NaN( size(interpolated_connectivity_list,1),size(scaler_field,2));
            %all_points_new_mesh = NaN( total_num_points_new_mesh,3 );
            %all_points_new_mesh(1:size(all_points,1),:) = all_points;
            edge_center_coors = (all_points(edges_real_mesh_unique(:,1) ,:) + all_points(edges_real_mesh_unique(:,2) ,:))/2;
            
            all_points_new_mesh = [ all_points;circum_centers;edge_center_coors ] ;
            
            
            for tri = 1:size(triangles,1)
                current_tri = triangles(tri,:);
                current_circum_center = circum_centers(tri,:);
                
                point_1 = all_points(current_tri(1),:);scaler_field1 = scaler_field(current_tri(1),:);
                point_2 = all_points(current_tri(2),:);scaler_field2 = scaler_field(current_tri(2),:);
                point_3 = all_points(current_tri(3),:);scaler_field3 = scaler_field(current_tri(3),:);
                
                edge12_center = (point_2 + point_1)/2; dis_circum_e12 = vecnorm( (current_circum_center - edge12_center) ,2,2);
                edge23_center = (point_3 + point_2)/2; dis_circum_e23 = vecnorm( (current_circum_center - edge23_center) ,2,2);
                edge31_center = (point_1 + point_3)/2; dis_circum_e31 = vecnorm( (current_circum_center - edge31_center) ,2,2);
                
                scaler_field_value_ed12 = (scaler_field( current_tri(1),: )+scaler_field( current_tri(2),: ))/2;
                scaler_field_value_ed23 = (scaler_field( current_tri(2),: )+scaler_field( current_tri(3),: ))/2;
                scaler_field_value_ed31 = (scaler_field( current_tri(3),: )+scaler_field( current_tri(1),: ))/2;
                
                scaler_field_value_circum = 1/dis_circum_e12*scaler_field_value_ed12 + 1/dis_circum_e23*scaler_field_value_ed23 + 1/dis_circum_e31*scaler_field_value_ed31;
                scaler_field_value_circum = scaler_field_value_circum./sum( [1/dis_circum_e12,1/dis_circum_e23,1/dis_circum_e31 ] );
                
                % make the mesh
                all_points_current = [ point_1;... %1
                    point_2;... %2
                    point_3;... %3
                    current_circum_center;... %4
                    edge12_center;...         %5
                    edge23_center;...         %6
                    edge31_center ];          %7
                [~,point_ids] = ismember(all_points_current,all_points_new_mesh,'rows');
                
                interpolated_connectivity_list(6*(tri-1) + 1,:) = [ point_ids(1) point_ids(5) point_ids(4)  ];
                scaler_interpolated(6*(tri-1) + 1,:) = mean([ scaler_field1;scaler_field_value_ed12;scaler_field_value_circum],1);
                
                
                interpolated_connectivity_list(6*(tri-1) + 2,:) = [ point_ids(5) point_ids(2) point_ids(4)  ];
                scaler_interpolated(6*(tri-1) + 2,:) = mean([ scaler_field_value_ed12;scaler_field2;scaler_field_value_circum],1);
                
                interpolated_connectivity_list(6*(tri-1) + 3,:) = [ point_ids(2) point_ids(6) point_ids(4)  ];
                scaler_interpolated(6*(tri-1) + 3,:) = mean([ scaler_field2;scaler_field_value_ed23;scaler_field_value_circum],1);
                
                interpolated_connectivity_list(6*(tri-1) + 4,:) = [ point_ids(6) point_ids(3) point_ids(4)  ];
                scaler_interpolated(6*(tri-1) + 4,:) = mean([ scaler_field_value_ed23;scaler_field3;scaler_field_value_circum],1);
                
                interpolated_connectivity_list(6*(tri-1) + 5,:) = [ point_ids(3) point_ids(7) point_ids(4)  ];
                scaler_interpolated(6*(tri-1) + 5,:) = mean([ scaler_field3;scaler_field_value_ed31;scaler_field_value_circum],1);
                
                interpolated_connectivity_list(6*(tri-1) + 6,:) = [ point_ids(7) point_ids(1) point_ids(4)  ];
                scaler_interpolated(6*(tri-1) + 6,:) = mean([ scaler_field_value_ed31;scaler_field1;scaler_field_value_circum],1);
                fprintf('FUNCTION: interpolate_vertex_scalers_to_mesh percentage_complete = %d\n',tri/size(triangles,1)*100);
                
            end
            
            interpolated_mesh = triangulation(interpolated_connectivity_list,all_points_new_mesh);
            
            
        end
        
        function [avg_free_energy_at_vertex,mole_lipid_comp,order_parameter,hmean,free_energy_curve,H_space,comp_space]  = ...
                mark_obj_free_energy_plot(obj_list,spont_curv_lipids, kappa, surface_modulous,kappa_coefficient,vanderwaals, kB, temperature)
            %             [coors] = extract_coors_points_ID_specific(obj_list);
            avg_free_energy_at_vertex = NaN(size(obj_list,1),1);
            %             order_parameter = NaN(size(obj_list,1),1);
            [coors,~,hmean,~,mole_lipid_comp,order_parameter]  = extract_patch_properties(obj_list,spont_curv_lipids);
            for obj_no = 1:size(obj_list,1)
                obj = obj_list(obj_no,1);
                [avg_free_energy_at_vertex(obj_no)] = ...
                    calculate_free_energy_at_vertex(obj, obj.Pos,obj.n_particles_up,coors,spont_curv_lipids, kappa, surface_modulous,kappa_coefficient,vanderwaals, kB, temperature);
                avg_free_energy_at_vertex(obj_no) = avg_free_energy_at_vertex(obj_no)/sum(obj.n_particles_up);
            end
            H_space = spont_curv_lipids(1):diff(spont_curv_lipids)/100:spont_curv_lipids(2);
            comp_space = 0:0.01:1;
            [H_mesh,comp_mesh] = meshgrid(H_space,comp_space);
            free_energy_curve = 1*0.5*kappa.*( H_mesh - spont_curv_lipids(1).*comp_mesh - spont_curv_lipids(2).*(1-comp_mesh) ).^2 + ...
                0.5*kappa*comp_mesh.*( spont_curv_lipids(1).*comp_mesh + spont_curv_lipids(2).*(1-comp_mesh) - spont_curv_lipids(1)).^2 + ...
                0.5*kappa*(1-comp_mesh).*( spont_curv_lipids(1).*comp_mesh + spont_curv_lipids(2).*(1-comp_mesh) - spont_curv_lipids(2)).^2 + ...
                +vanderwaals*comp_mesh.*(1-comp_mesh)...
                -temperature.*kB.*( - comp_mesh.*log(comp_mesh) -  (1-comp_mesh).*log(1-comp_mesh) + 1);
            
        end
    end
    
end