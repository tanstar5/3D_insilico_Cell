function [barbed_end_mesh] = initialize_barbed_end_mesh_from_pointed_end(pointed_end_mesh,length_actin_filaments)
%% Calculating the face direction vectors 
Face_direction_vects = faceNormal(pointed_end_mesh,[1:size(pointed_end_mesh.ConnectivityList,1)]');
pointed_end_mesh_points = pointed_end_mesh.Points;
%% Calculating the neighbourhood idx for each point in the pointed end mesh
faces = pointed_end_mesh.ConnectivityList;
face_idx = (1:size(faces,1))';
points = pointed_end_mesh.Points;
points_idx = (1:size(points,1))';
connected_faces_noter = cell(size(points,1),1);
direction_vects = NaN(size(points,1),3);
barbed_mesh_points = NaN(size(points,1),3);

for point_idth = 1:size(points,1)
    current_points_idx = points_idx(point_idth);
    presence_mask = ismember(faces,current_points_idx);
    presence_mask_logical = presence_mask(:,1)+ presence_mask(:,2)+ presence_mask(:,3);
    connected_faces_noter{point_idth,1} = face_idx(logical(presence_mask_logical));
    direction_vects(point_idth,:) = sum(Face_direction_vects(connected_faces_noter{point_idth,1},:),1);
    direction_vects(point_idth,:) = direction_vects(point_idth,:)/vecnorm(direction_vects(point_idth,:),2,2);
    barbed_mesh_points(point_idth,:) = pointed_end_mesh_points(point_idth,:) + mean( length_actin_filaments(connected_faces_noter{point_idth,1},:),1).*direction_vects(point_idth,:);
end

barbed_end_mesh = triangulation(faces,barbed_mesh_points);
end