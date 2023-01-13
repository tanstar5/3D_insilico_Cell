input_file_name = 'circle_rad_20_and_5multi.stl';
destinationFolder = 'C:\Users\tanum\OneDrive\Documents\PhD work\workfolder_surface\Self_Organization_project\active_cell_version42_two_component_gradient_termUpdate\outputThesis';

%% Physical parameters
kB = 1;
lipid_head_area = 0.5; % nm^2
kappa = 27*kB*300;
surface_modulous = kappa/4*12; %% relation between kappa and surface modulous is kappa = surface_molous*t^2 where t is bilayer thickness (wikipedia: search Lipid bilayer mechanics)
% or in popular literature kappa =
% surface_modulous*t^2/12
kappa_coefficient = 1;
vanderwaals = 2*kB*300;

temperature = 300;
penalty_factor = 8*kB*300;

radius = 300;
spont_curv_lipids = [-0.2+1/radius,1/radius+0.2];
spont_curv_lipids = [-1/60,1/60];

time_steps = 1;
step_size = 0.3;
particle_packet_size = 30;
patch_min_edge = 1.4;
patch_max_edge = 1.4;


obj = ...
    parameters(kB,lipid_head_area,kappa,surface_modulous,kappa_coefficient,vanderwaals,penalty_factor, ...
    temperature,spont_curv_lipids,radius,time_steps,step_size,particle_packet_size,patch_min_edge,patch_max_edge);

%% Create parameters for all experimental conditions
parameters_obj = obj;
% parameters_obj1 = obj;
% parameters_obj2 = obj;
% parameters_obj3 = obj;
% parameters_obj4 = obj;
% parameters_obj5 = obj;
% parameters_obj6 = obj;
% [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj,input_file_name,destinationFolder,5,1,1,1);
% % [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj1,input_file_name,destinationFolder,3,1,1);
% [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj,input_file_name,destinationFolder,5,1,1,0);
% delete(gcp('nocreate')); %delete the current pool

% [patch_min_edge_scaled,patch_max_edge_scaled,obj_list] = run_simulationv2(parameters_obj,input_file_name,destinationFolder,5,1,0,0);
[membrane_mesh] = import_surface_mesh_from_stl(input_file_name,patch_min_edge,patch_max_edge);
% trisurf(membrane_mesh),daspect([1,1,1]);
membrane_mesh_points = membrane_mesh.Points;
membrane_mesh_points_scaled = membrane_mesh_points/20*radius;
membrane_mesh = triangulation(membrane_mesh.ConnectivityList,membrane_mesh_points_scaled);
patch_min_edge_scaled  = patch_min_edge/20*radius;
patch_max_edge_scaled  = patch_max_edge/20*radius;

obj_mesh = membrane2D_continuous(1,[],[],[],[]);
[obj_list] = initilize_obj_list_container_from_mesh(obj_mesh,membrane_mesh);
[obj_list] = populate_with_lipids_at_equillibrium(obj_list,lipid_head_area,spont_curv_lipids);

[IDs] = extract_IDs(obj_list);
[lipid_comp] = extract_lipids_composition_ID_specific(obj_list);
[coors,principle_curvatures,hmean,lipid_comp,mole_lipid_comp,order_parameter]  = extract_patch_properties(obj_list,spont_curv_lipids);
all_scaler_mat = [ order_parameter,mole_lipid_comp(:,2),hmean];
[scaler_interpolated,interpolated_mesh,original_mesh] = interpolate_vertex_scalers_to_mesh(obj_list,all_scaler_mat);




pos_objs = zeros(size(obj_list,1),3);
grad_vect_along_surf_all_lq = zeros(size(obj_list,1),3);
grad_vect_along_surf_all_gauss = zeros(size(obj_list,1),3);
for i = 1:size(obj_list,1)
    obj = obj_list(i,1);
    obj.neighbours_lipid_compositions_up = lipid_comp(obj.neighbours,:);
    pos_objs(i,:) = obj.Pos;
    [~,grad_vect_along_surf_lq] = calculate_gradient_least_square_method_modified(obj);
    [grad_vect_along_surf_gauss] = calculate_gradient_green_gauss_method(obj);
    grad_vect_along_surf_all_lq(i,:) = grad_vect_along_surf_lq;
    grad_vect_along_surf_all_gauss(i,:) = grad_vect_along_surf_gauss;
end
trisurf(original_mesh,'FaceColor','r','FaceAlpha',0,'EdgeAlpha',0.1); daspect([1,1,1]);
hold on;trisurf(interpolated_mesh,scaler_interpolated(:,2),'FaceAlpha',.9,'EdgeAlpha',0.0);caxis([0.2,0.8]);colorbar;colormap(jet)
% hold on;q = quiver3(pos_objs(:,1),pos_objs(:,2),pos_objs(:,3),grad_vect_along_surf_all_lq(:,1),grad_vect_along_surf_all_lq(:,2),grad_vect_along_surf_all_lq(:,3),4,'r');q.Marker = '.';q.MaxHeadSize = 5; q.LineWidth = .5 ; daspect([1,1,1]);hold off
hold on;q = quiver3(pos_objs(:,1),pos_objs(:,2),pos_objs(:,3),grad_vect_along_surf_all_gauss(:,1),grad_vect_along_surf_all_gauss(:,2),grad_vect_along_surf_all_gauss(:,3),4,'b');q.Marker = '.';q.MaxHeadSize = 5; q.LineWidth = .5 ; daspect([1,1,1]);hold off
view(90,90);
set(gca,'FontSize',16)