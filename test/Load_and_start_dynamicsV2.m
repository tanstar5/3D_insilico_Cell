%% Load and Initialize cell body
cell1 = active_cell_object(1);
actin_length = .05;
actin_mesh_max_length = 5.0;
actin_mesh_min_length = 0.6*actin_mesh_max_length;
cell1 = create_initial_actin_mesh_pointed_ends_from_stl(cell1,1,actin_length,actin_mesh_max_length,actin_mesh_min_length); %% function active_cell_bject
hold on;
unit_membrane_mesh_length = actin_mesh_min_length/5;
cell1 = create_initial_membrane_mesh_barbed_ends(cell1,1,unit_membrane_mesh_length);
hold off

%% Initialize Membrane mesh
membrane_mesh = cell1.membrane.GlobalMesh;
[membrane_mesh_list] = ...
    membrane_particles_list(membrane_particle(1),size( membrane_mesh.Points,1 ));

[membrane_mesh_list] = assign_local_info_to_objects_fromMesh(membrane_mesh_list,membrane_mesh);
hold on;trisurf(membrane_mesh,'FaceColor','red','FaceAlpha',.5);daspect([1,1,1]);

%% Assign and measure initial membrane properties
temperature = 65;
pressure = 60;
surface_tension = 00;
hardball_rad = unit_membrane_mesh_length*.30;
% kick_displacement = hardball_rad*.2;
kick_displacement = (unit_membrane_mesh_length-2*hardball_rad)/3;
normal_displacement = 2*kick_displacement;
volume_constraint_coefficient = 65*4;
V_initial = 4/3*pi*15^3;
k_a = 100;
% hardball_rad = kick_displacement*.3;

% type_list =       [0,     1,  2,     3,    4     5,  6];
% type_Co =         [-1/10, 0,  0,     1/10, 1/5, 1/3, 1/2];
% type_rigidities = [200,     100,  100,     200,    200,   200    500];

type_list =       [0,     1     2      3];
type_Co =         [-1/5,  0     0      1/5];
type_rigidities = [150,   50    50     150 ];

[membrane_mesh_list] = assign_type_randomly_conserved(membrane_mesh_list,type_list);
% [membrane_mesh_list,list_id_logical] = assign_type_poles(membrane_mesh_list,1,8,2);
[membrane_mesh_list] = assign_type_properties(membrane_mesh_list,type_list,type_Co,type_rigidities);
% [membrane_mesh_list,tot_energy] = calculate_energy_stored(membrane_mesh_list,hardball_rad,100);
[mean_curvature] = mean_curvature_calc(membrane_mesh_list);
% energy_extract = @(id) membrane_mesh_list(id).energy_stored;
% energy_all_mesh = arrayfun(energy_extract,1:size(membrane_mesh_list,1));
% hold on; scatter3(membrane_mesh.Points(:,1),membrane_mesh.Points(:,2),membrane_mesh.Points(:,3),...
%     10,energy_all_mesh,'filled');
% caxis([min(energy_all_mesh) max(energy_all_mesh)]);hold off
[diffused_mesh] = generate_mesh_from_objlist(membrane_mesh_list);

hold on; trisurf(diffused_mesh,'FaceColor','blue','FaceAlpha',.5);hold off;
diffused_membrane_mesh_list = membrane_mesh_list;
[diffused_membrane_mesh_list] = update_neighbours_normals(diffused_membrane_mesh_list);
[diffused_membrane_mesh_list] = allign_faceOrder_based_on_faceNormal(diffused_membrane_mesh_list);

%% Select randomly all type of objects and plot for energy comparison
[type_all] = extract_type(diffused_membrane_mesh_list);
all_IDs = extract_IDs(diffused_membrane_mesh_list);
no_of_types = length(type_list);
[type_all_Ho] = extract_spontaneousCurvature(diffused_membrane_mesh_list);
hold on; scatter3(diffused_mesh.Points(:,1),diffused_mesh.Points(:,2),diffused_mesh.Points(:,3),...
    15,type_all_Ho,'filled');
for type_ind = 1:no_of_types
    IDs_current_type = all_IDs(type_all==type_list(type_ind));
    random_id_selected = IDs_current_type(randi(length(IDs_current_type)));
    obj_to_plot = findobjID(diffused_membrane_mesh_list,random_id_selected);
    obj_to_plot.ID
    compare_energies(obj_to_plot,kick_displacement*10,pressure,surface_tension,temperature,hardball_rad,4);
    hold on
end
hold off;

%% Running iterations for diffusion
for df = 0:10000
    h1 = figure(1);
    fprintf('Loop t = %d Started\n',df);
    
    %% post processing
    [diffused_mesh] = generate_mesh_from_objlist(diffused_membrane_mesh_list);
    [face_curvature] = generate_curvature_plotSurface(diffused_membrane_mesh_list,diffused_mesh);
    [type_all_Ho] = extract_spontaneousCurvature(diffused_membrane_mesh_list);
    face_curvature_color = (face_curvature-min(face_curvature(:)))/(max(face_curvature(:))-min(face_curvature(:)))*( max(type_all_Ho)-min(type_all_Ho) )+min(type_all_Ho);
    
    %% Plotting
    %     if mod(df,5) == 0
    trisurf(membrane_mesh,'FaceColor','red','FaceAlpha',.3,'EdgeAlpha',0.2);daspect([1,1,1]);hold on;
    %     end
    trisurf(diffused_mesh,face_curvature,'FaceAlpha',.6,'EdgeAlpha',0.2);daspect([1,1,1]);hold off;
    if mod(df,1) == 0
        hold on; scatter3(diffused_mesh.Points(:,1),diffused_mesh.Points(:,2),diffused_mesh.Points(:,3),...
            20,type_all_Ho,'filled');
    end
    caxis([min(type_all_Ho) max(type_all_Ho)]);hold off;
    view([90,0]);xlim([0,18]);ylim([-15,15]);zlim([-15,15]);
    
    %% Saving into file
    set(h1, 'Position', [10 10 900 900]);
    current_frame = getframe(h1);
    imwrite(current_frame.cdata,'Membrane.tif','WriteMode','append');
    
    [free_energy,curvature_energy] =...
        calculate_energy_ofType(diffused_membrane_mesh_list,type_list,pressure,surface_tension,temperature,hardball_rad,1,3);
    set(gcf, 'Position', [10 10 500 500]);
    current_frame = getframe(gcf);
    imwrite(current_frame.cdata,'EnergyPlot.tif','WriteMode','append');
    if df == 0
       [V_initial,area_total] = area_volume_surface(diffused_mesh.Points',diffused_mesh.ConnectivityList'); 
       area_reduced = area_total/size(diffused_membrane_mesh_list,1);
    end
    
    %% checking mesh condition
    if mod(df,3) == 0
        h4 = figure(4);
        close(h4)
        all_IDs = extract_IDs(diffused_membrane_mesh_list);
        for type_ind = 1:no_of_types
            IDs_current_type = all_IDs(type_all==type_list(type_ind));
            random_id_selected = IDs_current_type(randi(length(IDs_current_type)));
            obj_to_plot = findobjID(diffused_membrane_mesh_list,random_id_selected);
            obj_to_plot.ID
            compare_energies(obj_to_plot,kick_displacement,pressure,surface_tension,temperature,hardball_rad,4);
            hold on
        end
        [diffused_membrane_mesh_list] = ...
                enhance_structures_based_on_curvature(diffused_membrane_mesh_list,normal_displacement,surface_tension,hardball_rad,100,pressure,temperature,volume_constraint_coefficient,V_initial);
    end
    
    pause(.1);
    %% Diffusing
    [diffused_membrane_mesh_list,acceptance_rate_orthogonal] = diffuse_fluctuation_orthogonal_instantaneous(diffused_membrane_mesh_list,hardball_rad,kick_displacement,100,temperature,pressure,surface_tension,area_reduced,k_a);
    [diffused_membrane_mesh_list,acceptance_rate_lateralDiffusion] = diffuse_vertex_flip(diffused_membrane_mesh_list,temperature,hardball_rad,pressure,surface_tension);
    [diffused_membrane_mesh_list,acceptance_rate_fluidity] = diffuse_edge_flip(diffused_membrane_mesh_list,hardball_rad,100,temperature,pressure,surface_tension);
    fprintf('Loop t = %d Completed\n\n\n',df);
end







