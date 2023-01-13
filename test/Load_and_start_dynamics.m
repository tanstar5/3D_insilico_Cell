%%% Protocol for simulating Cell

%% Load and Initialize cell body

cell1 = active_cell_object(1);
actin_length = .05;
cell1 = create_initial_actin_mesh_pointed_ends_from_stl(cell1,1,actin_length); %% function active_cell_bject
hold on;
unit_membrane_mesh_length = .15; 

cell1 = create_initial_membrane_mesh_barbed_ends(cell1,1,unit_membrane_mesh_length);
hold off

%% Initialize Membrane mesh
membrane_mesh = cell1.membrane.GlobalMesh;
[membrane_mesh_list] = ...
    membrane_particles_list(membrane_particle(1),size( membrane_mesh.Points,1 ));

[membrane_mesh_list] = assign_local_info_to_objects_fromMesh(membrane_mesh_list,membrane_mesh);
trisurf(membrane_mesh,'FaceColor','red','FaceAlpha',.5);daspect([1,1,1]);
% [membrane_mesh_list,list_id_logical] = assign_type_poles(membrane_mesh_list,8,2);


%% Assign and measure initial membrane properties
temperature = 1;
pressure = 0;
hardball_rad = unit_membrane_mesh_length /6;
kick_displacement = hardball_rad*2;
type_list = [0,1,2];
type_Co = [1/5,0,-1/5];
type_rigidities = [50,10,50];
[membrane_mesh_list] = assign_type_randomly_conserved(membrane_mesh_list,type_list);
% [membrane_mesh_list,list_id_logical] = assign_type_poles(membrane_mesh_list,1,8,2);
[membrane_mesh_list] = assign_type_properties(membrane_mesh_list,type_list,type_Co,type_rigidities);
[membrane_mesh_list,tot_energy] = calculate_energy_stored(membrane_mesh_list,hardball_rad,100);
[mean_curvature] = mean_curvature_calc(membrane_mesh_list);
energy_extract = @(id) membrane_mesh_list(id).energy_stored;
energy_all_mesh = arrayfun(energy_extract,1:size(membrane_mesh_list,1));
hold on; scatter3(membrane_mesh.Points(:,1),membrane_mesh.Points(:,2),membrane_mesh.Points(:,3),...
    10,energy_all_mesh,'filled');
caxis([min(energy_all_mesh) max(energy_all_mesh)]);hold off

%% Start Dynamics
% temperature = 1;
pressure = 5;

[diffused_mesh] = generate_mesh_from_objlist(membrane_mesh_list);
hold on; trisurf(diffused_mesh,'FaceColor','blue','FaceAlpha',.5);hold off;
diffused_membrane_mesh_list = membrane_mesh_list;
surface_tension = 0.01;

for df = 1:10000
    tic;
    figure(1);
    fprintf('Loop t = %d Started\n',df);
%     diffuse_fluctuation_orthogonal(obj_list,hardball_radius,kick_displacement,overlapping_energy,temperature_current,pressure,surface_tension)
    [diffused_membrane_mesh_list] = diffuse_fluctuation_orthogonal(diffused_membrane_mesh_list,hardball_rad,kick_displacement,100,temperature,pressure,surface_tension);
    

%     [diffused_membrane_mesh_list] = diffuse_edge_flip(diffused_membrane_mesh_list,hardball_rad,100,temperature,pressure,surface_tension);
    [diffused_membrane_mesh_list] = diffuse_vertex_flip(diffused_membrane_mesh_list,temperature,hardball_rad,pressure,surface_tension);
    [diffused_membrane_mesh_list] = diffuse_edge_flip(diffused_membrane_mesh_list,hardball_rad,100,temperature,pressure,surface_tension);
    [diffused_membrane_mesh_list] = diffuse_edge_flip(diffused_membrane_mesh_list,hardball_rad,100,temperature,pressure,surface_tension);
    [diffused_membrane_mesh_list] = diffuse_edge_flip(diffused_membrane_mesh_list,hardball_rad,100,temperature,pressure,surface_tension);
%     [diffused_membrane_mesh_list] = diffuse_vertex_flip(diffused_membrane_mesh_list,temperature,pressure);
%     [diffused_membrane_mesh_list] = diffuse_edge_flip(diffused_membrane_mesh_list,hardball_rad,100,temperature,pressure);
%     [diffused_membrane_mesh_list] = diffuse_vertex_flip(diffused_membrane_mesh_list,temperature,pressure);
%     [diffused_membrane_mesh_list] = diffuse_edge_flip(diffused_membrane_mesh_list,hardball_rad,100,temperature,pressure);
%     [diffused_membrane_mesh_list] = allign_faceOrder_based_on_faceNormal(diffused_membrane_mesh_list);
    [diffused_mesh] = generate_mesh_from_objlist(diffused_membrane_mesh_list);
    [face_curvature] = generate_curvature_plotSurface(diffused_membrane_mesh_list,diffused_mesh);
    
    
%     trisurf(membrane_mesh,'FaceColor','red','FaceAlpha',.1);daspect([1,1,1]);
%     hold on;
%     [type_all] = extract_type(diffused_membrane_mesh_list);
    [type_all_Ho] = extract_spontaneousCurvature(diffused_membrane_mesh_list);
    face_curvature_color = (face_curvature-min(face_curvature(:)))/(max(face_curvature(:))-min(face_curvature(:)))*( max(type_all_Ho)-min(type_all_Ho) )+min(type_all_Ho);
    trisurf(diffused_mesh,face_curvature,'FaceAlpha',.6);daspect([1,1,1]);hold off;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    energy_extract = @(id) diffused_membrane_mesh_list(id).energy_stored;
    
    hold on; scatter3(diffused_mesh.Points(:,1),diffused_mesh.Points(:,2),diffused_mesh.Points(:,3),...
        30,type_all_Ho,'filled');
    caxis([min(type_all_Ho) max(type_all_Ho)]);hold off;
    view([-180,4]);xlim([-3.5,3.5]);ylim([-0.5,6.5]);zlim([-3.5,3.5]);%colormap(parula(length(type_list)));
    
    current_frame = getframe(gcf);
    imwrite(current_frame.cdata,'Membrane.tif','WriteMode','append');
    
    [free_energy,curvature_energy] =...
        calculate_energy_ofType(diffused_membrane_mesh_list,type_list,pressure,surface_tension,temperature,hardball_rad,1,3);
    current_frame = getframe(gcf);
    imwrite(current_frame.cdata,'EnergyPlot.tif','WriteMode','append');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    pause(.1);
    toc;
    fprintf('Loop t = %d Completed\n\n\n',df);
end


