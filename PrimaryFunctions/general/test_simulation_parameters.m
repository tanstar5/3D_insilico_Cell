function [] = test_simulation_parameters(parameters_obj,input_file_name,destinationFolder,figure_no)

% constants
kB = parameters_obj.kB;
lipid_head_area = parameters_obj.lipid_head_area;

% physical parameters
kappa = parameters_obj.kappa;
surface_modulous = parameters_obj.surface_modulous;
penalty_factor = parameters_obj.penalty_factor;
temperature = parameters_obj.temperature;
spont_curv_lipids = parameters_obj.spont_curv_lipids;
radius = parameters_obj.radius;
kappa_coefficient = parameters_obj.kappa_coefficient;
vanderwaals = parameters_obj.vanderwaals;

% simulation parameters
time_steps = parameters_obj.time_steps;
step_size = parameters_obj.step_size;
particle_packet_size = parameters_obj.particle_packet_size;
patch_min_edge = parameters_obj.patch_min_edge;
patch_max_edge = parameters_obj.patch_max_edge;
fig_limits = radius + 0.20*radius;

output_file_name = sprintf( 'test_exp_no_%d_temperature_%d_radius_%d.fig',figure_no,temperature,radius);

if ~exist(destinationFolder, 'dir')
    mkdir(destinationFolder);
end
fullDestinationFileName = fullfile(destinationFolder, output_file_name);


[membrane_mesh] = import_surface_mesh_from_stl(input_file_name,patch_min_edge,patch_max_edge);
% trisurf(membrane_mesh),daspect([1,1,1]);
%% Due to inconsistency in stl file we scale it keepin the number of points same based on radius
membrane_mesh_points = membrane_mesh.Points;
membrane_mesh_points_scaled = membrane_mesh_points/20*radius;
membrane_mesh = triangulation(membrane_mesh.ConnectivityList,membrane_mesh_points_scaled);
patch_min_edge_scaled  = patch_min_edge/20*radius;
patch_max_edge_scaled  = patch_max_edge/20*radius;

obj = membrane2D_continuous(1,[],[],[],[]);
[obj_list] = initilize_obj_list_container_from_mesh(obj,membrane_mesh);
[obj_list] = populate_with_lipids_at_equillibrium(obj_list,lipid_head_area,spont_curv_lipids);

% check geometric calculations
h1 = figure(figure_no);
[obj_list] = derive_geometrical_quatities_all(obj_list);
% [mean_curvature,principle_curvatures] = extract_curvature_list(obj_list);

[coors,principle_curvatures,hmean,lipid_comp,mole_lipid_comp,order_parameter]  = extract_patch_properties(obj_list,spont_curv_lipids);
all_scaler_mat = [ order_parameter,mole_lipid_comp(:,2),hmean];
[scaler_interpolated,interpolated_mesh,original_mesh] = interpolate_vertex_scalers_to_mesh(obj_list,all_scaler_mat);
%% plot the results
h1_sub1 = subplot(4,3,[7,8,10,11]);
axes(h1_sub1);trisurf(original_mesh,'FaceColor','r','FaceAlpha',0,'EdgeAlpha',0.1); daspect([1,1,1]);title(h1_sub1,'Curvature');colormap(h1_sub1,'parula');hold off;
hold on;trisurf(interpolated_mesh,scaler_interpolated(:,3),'FaceAlpha',1,'EdgeAlpha',0.1);caxis(spont_curv_lipids);colorbar;colormap(parula);hold off
ylim([-fig_limits,fig_limits]);xlim([-fig_limits,fig_limits]);zlim([-fig_limits,fig_limits]);
view(0,0);ylim([-fig_limits,0]);xlim([-fig_limits,fig_limits]);zlim([-fig_limits,fig_limits]);hold off;

h1_sub2 = subplot(4,3,3);
axes(h1_sub2);trisurf(original_mesh,'FaceColor','r','FaceAlpha',0,'EdgeAlpha',0.1); daspect([1,1,1]);title(h1_sub2,'order parameter');colormap(h1_sub2,'parula');hold off;
hold on;trisurf(interpolated_mesh,scaler_interpolated(:,1),'FaceAlpha',1,'EdgeAlpha',0.1);caxis(spont_curv_lipids);colorbar;colormap(parula);hold off
ylim([-fig_limits,fig_limits]);xlim([-fig_limits,fig_limits]);zlim([-fig_limits,fig_limits]);hold off;

h1_sub3 = subplot(4,3,[1,2,4,5]);
axes(h1_sub3);trisurf(original_mesh,'FaceColor','r','FaceAlpha',0,'EdgeAlpha',0.1); daspect([1,1,1]);title(h1_sub3,'composition');colormap(h1_sub3,'parula');hold off;
hold on;trisurf(interpolated_mesh,scaler_interpolated(:,2),'FaceAlpha',1,'EdgeAlpha',0.1);caxis([0.2,0.8]);colorbar;colormap(parula);hold off
view(0,0);ylim([-fig_limits,0]);xlim([-fig_limits,fig_limits]);zlim([-fig_limits,fig_limits]);caxis([0.2,0.8]);colorbar;hold off;
pause(0.1);

[avg_free_energy_at_vertex,mole_lipid_comp,order_parameter,hmean,free_energy_curve,H_space,comp_space]  = ...
    mark_obj_free_energy_plot(obj_list,spont_curv_lipids, kappa, surface_modulous,kappa_coefficient,vanderwaals, kB, temperature);




h1_sub4 = subplot(4,3,[6,9,12]);
surf(H_space,comp_space,free_energy_curve,'FaceColor','b','FaceAlpha',0.5);
hold on; scatter3(hmean,mole_lipid_comp(:,1),avg_free_energy_at_vertex,20,'g');hold off ;view(0,0);
xlabel('curvature');ylabel('mole_fraction');zlabel('avg free energy');

suptitle(sprintf('experiment%d frame no:%d temp: %d radius: %d\n spontcurv [%d %d]',figure_no,0,temperature,radius,spont_curv_lipids(1),spont_curv_lipids(2)))
%% Saving into file
savefig(h1,fullDestinationFileName);

end