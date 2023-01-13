temperature = 1000;
kick_displacement = 0.05;
hardball_rad = 0.2;

trisurf(membrane_mesh,'FaceColor','red','FaceAlpha',.5);daspect([1,1,1]);
hold on; scatter3(membrane_mesh.Points(:,1),membrane_mesh.Points(:,2),membrane_mesh.Points(:,3),...
    10,'r','filled');
% [diffused_membrane_mesh_list] = diffuse_fluctuation(membrane_mesh_list,hardball_rad,kick_displacement,100,temperature);
diffused_membrane_mesh_list = membrane_mesh_list;
% [coors] = extract_coors_points(diffused_membrane_mesh_list);
% hold on; scatter3(coors(:,1),coors(:,2),coors(:,3),...
%     10,'b','filled');hold off;
[diffused_mesh] = generate_mesh_from_objlist(membrane_mesh_list);
hold on; trisurf(diffused_mesh,'FaceColor','blue','FaceAlpha',.5);hold off;

for df = 1:100
    [diffused_membrane_mesh_list] = diffuse_fluctuation(diffused_membrane_mesh_list,hardball_rad,kick_displacement,100,temperature);
    [diffused_membrane_mesh_list] = diffuse_vertex_flip(diffused_membrane_mesh_list,temperature);
    [diffused_mesh] = generate_mesh_from_objlist(diffused_membrane_mesh_list);
    trisurf(membrane_mesh,'FaceColor','red','FaceAlpha',.5);daspect([1,1,1]);
    hold on;
    trisurf(diffused_mesh,'FaceColor','blue','FaceAlpha',.5);daspect([1,1,1]);hold off;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    energy_extract = @(id) diffused_membrane_mesh_list(id).energy_stored;
%     energy_all_mesh = arrayfun(energy_extract,1:size(diffused_membrane_mesh_list,1));
    [type_all] = extract_type(diffused_membrane_mesh_list);
%     color = zeros(size(list_id_logical));
%     color(list_id_logical==1) = 10;
%     color(list_id_logical==0) = -10;
    hold on; scatter3(diffused_mesh.Points(:,1),diffused_mesh.Points(:,2),diffused_mesh.Points(:,3),...
        10,type_all,'filled');
    caxis([min(type_all) max(type_all)]);hold off;
    view(2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    pause(.1);
end