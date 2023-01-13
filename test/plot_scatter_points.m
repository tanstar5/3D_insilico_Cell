% digits(5)
membrane_mesh = cell1.membrane.GlobalMesh;
[membrane_mesh_list] = ...
    membrane_particles_list(membrane_particle(1),size( membrane_mesh.Points,1 ));

[membrane_mesh_list] = assign_local_info_to_objects_fromMesh(membrane_mesh_list,membrane_mesh);
trisurf(membrane_mesh,'FaceColor','red','FaceAlpha',.5);daspect([1,1,1]);
% [membrane_mesh_list,list_id_logical] = assign_type_poles(membrane_mesh_list,8,2);
[membrane_mesh_list] = assign_type_randomly_conserved(membrane_mesh_list,[0,2,3]);
temperature = 1000;
kick_displacement = 0.05;
hardball_rad = 0.2;


[membrane_mesh_list] = assign_type_properties(membrane_mesh_list,[0,2,3],[1/1000,1/500,-1/500],[1,1,1]);
[membrane_mesh_list,tot_energy] = calculate_energy_stored(membrane_mesh_list,hardball_rad,100);
% mean_curvature_color = @(id) mean(membrane_mesh_list(id).DarbouxFrame{2,1});
% mean_curvature_at
[mean_curvature] = mean_curvature_calc(membrane_mesh_list);
energy_extract = @(id) membrane_mesh_list(id).energy_stored;
energy_all_mesh = arrayfun(energy_extract,1:size(membrane_mesh_list,1));
% color = zeros(size(list_id_logical));
% color(list_id_logical==1) = 10;
% color(list_id_logical==0) = -10;
hold on; scatter3(membrane_mesh.Points(:,1),membrane_mesh.Points(:,2),membrane_mesh.Points(:,3),...
    10,energy_all_mesh,'filled');
caxis([min(energy_all_mesh) max(energy_all_mesh)]);

% obj = membrane_mesh_list(5366);
% %% plotting neighbouring vertex
% for ne = 1:length(obj.neighbours)
% hold on;
% scatter3(membrane_mesh.Points(obj.neighbours(ne),1),...
%     membrane_mesh.Points(obj.neighbours(ne),2),...
%     membrane_mesh.Points(obj.neighbours(ne),3),160,'filled');
% 
% end
% 
% %% plotting normals
% hold on;
% face_centroid = NaN(size(obj.faces,1),3);
% for fc = 1:size(obj.faces,1)
%     current_face = obj.faces(fc,:);
%     fc_points = membrane_mesh.Points(current_face,:);
%     face_centroid(fc,:) = mean(fc_points,1);
% end
% quiver3(face_centroid(:,1), face_centroid(:,2), face_centroid(:,3), obj.faceNormals(:,1), obj.faceNormals(:,2), obj.faceNormals(:,3));
% 
% 
% hold on;
% scatter3(obj.Pos(1),...
%     obj.Pos(2),...
%     obj.Pos(3),160);
hold off;

