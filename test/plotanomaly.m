obj = findobjID(diffused_membrane_mesh_list,4571);
hold on; scatter3(obj.Pos(1),obj.Pos(2),obj.Pos(3),40,'filled');
hold on; quiver3(obj.Pos(1),obj.Pos(2),obj.Pos(3),obj.Normal_vect(1),obj.Normal_vect(2),obj.Normal_vect(3));

% %% checking the angles
% edge1_tri1 = edge_to_remove_objs(1).Pos - edge_to_form_objs(1).Pos;
% edge1_tri1_unit = edge1_tri1/vecnorm(edge1_tri1);
% edge2_tri1 = edge_to_remove_objs(2).Pos - edge_to_form_objs(1).Pos;
% edge2_tri1_unit = edge2_tri1/vecnorm(edge2_tri1);
% 
% angle_tri1 = acos(dot(edge1_tri1_unit,edge2_tri1_unit))*180/pi
% angle_bisector_tri1 = (edge1_tri1_unit+edge2_tri1_unit)/1;
% angle_bisector_tri1 = angle_bisector_tri1/vecnorm(angle_bisector_tri1);
% edge_to_form = edge_to_form_objs(2).Pos - edge_to_form_objs(1).Pos;
% edge_to_form_unit = edge_to_form/vecnorm(edge_to_form);
% angle_with_bisector_toformedge = acos(dot(angle_bisector_tri1,edge_to_form_unit))*180/pi
% 
% 
% decision = angle_with_bisector_toformedge<angle_tri1/2