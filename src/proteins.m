classdef proteins
    
   properties
       ID
       Pos
       structure
       size_protein
       in_or_out
   end
   
   methods
       
       %% Initialize functions       
       function obj = proteins(ID)
            if nargin ~= 0
                obj.ID = ID;
            end
       end
       
       function [obj_list] = proteins_list(first_obj,num_of_particles)
            obj_list(1,1) =  first_obj;
            for particle = 2:num_of_particles
                obj = proteins(particle);
                obj_list(particle,1) =  obj;
            end
       end
       
       function [updated_obj_list] = load_properties_of_proteins(obj_list,ID_list,Pos_list,structure_list,size_list,in_or_out_list)
           fprintf('Function:load_properties_of_proteins; loading properties of proteins\n');
           updated_obj_list = obj_list;
           for protein_ind = 1:length(ID_list)
              current_obj = obj_list(protein_ind);
              current_obj.ID = ID_list(protein_ind);
              current_obj.Pos = Pos_list(protein_ind,:);
              current_obj.structure = structure_list(protein_ind,:);
              current_obj.size_protein = size_list(protein_ind);
              current_obj.in_or_out = in_or_out_list(protein_ind);
              updated_obj_list(protein_ind)  = current_obj;
           end
           fprintf('Function:load_properties_of_proteins; loading properties COMPLETE of proteins\n');
            
       end
       
       %% Dynamics 
       
       %% Accecories functions
       function [] = protein_in_or_out(obj_list,membrane_mesh,outside_cell_coor,resolution)
           all_points_membrane_mesh = incenter(membrane_mesh);
           all_points_idx = 1:size(all_points_membrane_mesh,1);
           face_normal = faceNormal(membrane_mesh);
%            triangulation_surface = membrane_mesh.ConnectivityList;
%            calc_distance = @()
           for protein_id = 1:size(obj_list,1)
               current_protein = obj_list(protein_id,1);
               current_protein_pos = current_protein.Pos;
               
               direction_vect = outside_cell_coor - current_protein_pos;
               unit_direction_vect = direction_vect./vecnorm(direction_vect,2,2);
               
               x_coor_line = current_protein_pos(1):resolution*unit_direction_vect(1):outside_cell_coor(1);
               y_coor_line = current_protein_pos(2):resolution*unit_direction_vect(2):outside_cell_coor(2);
               z_coor_line = current_protein_pos(3):resolution*unit_direction_vect(3):outside_cell_coor(3);
               
               line_segment = [x_coor_line' y_coor_line' z_coor_line'];
               line_cross_indicator = zeros(size(line_segment,1),1);
               for point_in_seg = 1:size(line_segment,1)
                   current_point = line_segment(point_in_seg,:);
                   reference_mat = (current_point'*ones(1,size(all_points_membrane_mesh,1)))';
                   distance_array = vecnorm(all_points_membrane_mesh-reference_mat,2,2);
                   
                   %% Find points in its neighbourhood
                   neighbourhood_points_mask = distance_array<=2*resolution;
                   
                   %% from triangulation find a point of a random neighbour to know the surface
                   if sum(neighbourhood_points_mask)==0
                       line_cross_indicator(point_in_seg) = 0;
                   else
                       membrane_points_detected = all_points_idx(neighbourhood_points_mask);
                       corresponding_distances = distance_array(neighbourhood_points_mask);
                       [~,min_id] = min(corresponding_distances);
                       nearest_neighbour = membrane_points_detected(min_id(1));
                       if abs(dot(unit_direction_vect,face_normal(nearest_neighbour,:)))<=0.2
                           line_cross_indicator(point_in_seg) = 1;
                       else
                           line_cross_indicator(point_in_seg) = 0;
                       end
                        
                       
                   end
                       
                   
                   %% check orthogonality criteria
                   
                   %% determine if inside or outside
                                      
                   
                   
               end
               
               
               
           end
           
           
           
       end
       
       
        
        
       
       
   end
           
end