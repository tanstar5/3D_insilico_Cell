% diffused_membrane_mesh_list = membrane_mesh_list;
% 
% [diffused_membrane_mesh_list] = diffuse_fluctuation_orthogonal(diffused_membrane_mesh_list,hardball_rad,kick_displacement,100,temperature,pressure,surface_tension);
% [diffused_membrane_mesh_list] = diffuse_vertex_flip(diffused_membrane_mesh_list,temperature,hardball_rad,pressure,surface_tension);
% [diffused_membrane_mesh_list] = diffuse_edge_flip(diffused_membrane_mesh_list,hardball_rad,100,temperature,pressure,surface_tension);
% 
% diffused_membrane_mesh_list(1).Normal_vect = -diffused_membrane_mesh_list(1).Normal_vect;
obj = diffused_membrane_mesh_list(2);

[caz,cel] = view;
obj.Pos = obj.Pos - 0.3*obj.Normal_vect;
obj.Pos;
% figure()

scatter3(obj.Pos(1),obj.Pos(2),obj.Pos(3),'filled');
view(caz,cel);
hold on;scatter3(obj.neighbours_coor(:,1),obj.neighbours_coor(:,2),obj.neighbours_coor(:,3));hold off

% plotting the object with neighbours
all_IDs = [obj.ID,obj.neighbours'];
all_coors = [obj.Pos;obj.neighbours_coor];
current_faces = obj.faces;
faces_local = mask_mat_to_sequence(current_faces,all_IDs);

current_surf = triangulation(faces_local,all_coors);
hold on; trisurf(current_surf,'FaceColor','blue','FaceAlpha',.3);
hold on;quiver3(obj.Pos(1),obj.Pos(2),obj.Pos(3),obj.Normal_vect(1),obj.Normal_vect(2),obj.Normal_vect(3));
principle_vect1 = obj.DarbouxFrame{1,1}(:,1);
principle_vect2 = obj.DarbouxFrame{1,1}(:,2);
hold on;quiver3(obj.Pos(1),obj.Pos(2),obj.Pos(3),principle_vect1(1),principle_vect1(2),principle_vect1(3));
hold on;quiver3(obj.Pos(1),obj.Pos(2),obj.Pos(3),principle_vect2(1),principle_vect2(2),principle_vect2(3));
view(2);
daspect([1,1,1]);view(caz,cel);hold off;
[obj] = calculate_topological_characteristics(obj,hardball_rad,100);
[principle_curvatures] = extract_curvature(obj);
principle_curvatures(1,:)