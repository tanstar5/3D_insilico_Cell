 function [patch_min_edge_scaled,patch_max_edge_scaled,obj_list] = run_simulationv2(parameters_obj,input_file_name,destinationFolder,frame_steps,figure_no,uniformly_distribute_lipids,orthogonal_move_switch)
% constants
kB = parameters_obj.kB;
lipid_head_area = parameters_obj.lipid_head_area;

% physical parameters
kappa = parameters_obj.kappa;
surface_modulous = parameters_obj.surface_modulous;
penalty_factor = parameters_obj.penalty_factor;
kappa_coefficient = parameters_obj.kappa_coefficient;
vanderwaals = parameters_obj.vanderwaals;

temperature = parameters_obj.temperature;
spont_curv_lipids = parameters_obj.spont_curv_lipids;
radius = parameters_obj.radius;

% simulation parameters
time_steps = parameters_obj.time_steps;
step_size = parameters_obj.step_size;
particle_packet_size = parameters_obj.particle_packet_size;
patch_min_edge = parameters_obj.patch_min_edge;
patch_max_edge = parameters_obj.patch_max_edge;
fig_limits = radius + 0.20*radius;

output_file_name = sprintf( 'output_exp_no_%d_temperature_%d_radius_%d.tif',figure_no,temperature,radius);
output_fig_file_name = sprintf( 'output_exp_no_%d_temperature_%d_radius_%d.fig',figure_no,temperature,radius);

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
if uniformly_distribute_lipids == 0
    [obj_list] = populate_with_lipids_at_equillibrium(obj_list,lipid_head_area,spont_curv_lipids);
else
    [obj_list] = populate_with_lipids_uniform(obj_list,lipid_head_area,spont_curv_lipids);
end

% check geometric calculations
h1 = figure(figure_no);
[obj_list] = derive_geometrical_quatities_all(obj_list);
[normal_vects_all] = extract_Normal_vects(obj_list);
% [mean_curvature,principle_curvatures] = extract_curvature_list(obj_list);
[triangle_color_curvature] = color_code_patch_wise_curvature(obj_list);
[surface_mesh] = generate_mesh_from_objlist(obj_list);
h1_sub1 = subplot(2,2,1);
axes(h1_sub1);trisurf(surface_mesh,triangle_color_curvature,'FaceAlpha',.3,'EdgeAlpha',0.2); daspect([1,1,1]);title(h1_sub1,'Curvature');colormap(h1_sub1,'jet');hold off;
h1_sub2 = subplot(2,2,2);
histogram(h1_sub2,triangle_color_curvature)

[triangle_color_mole_fraction] = color_code_patch_wise_composition(obj_list);
h1_sub3 = subplot(2,2,3);
axes(h1_sub3);trisurf(surface_mesh,triangle_color_mole_fraction(:,2),'FaceAlpha',.3,'EdgeAlpha',0.2); daspect([1,1,1]);title(h1_sub3,'Composition');colormap(h1_sub3,'jet');hold off;

for t = 1:time_steps
    if orthogonal_move_switch==1
        [obj_list] = orthogonal_move(obj_list,step_size,spont_curv_lipids, kappa, surface_modulous,kappa_coefficient,vanderwaals, kB, temperature);
    end
    [obj_list] = exchange_move(obj_list,particle_packet_size,spont_curv_lipids, kappa, surface_modulous,penalty_factor,kappa_coefficient,vanderwaals, kB, temperature);
    [coors,principle_curvatures,hmean,lipid_comp,mole_lipid_comp,order_parameter]  = extract_patch_properties(obj_list,spont_curv_lipids);
    
    %% interpolating the scalers in an interpolated mesh
    if mod(t-1,frame_steps)==0
        all_scaler_mat = [ order_parameter,mole_lipid_comp(:,2),hmean];
        [scaler_interpolated,interpolated_mesh,original_mesh] = interpolate_vertex_scalers_to_mesh(obj_list,all_scaler_mat);
        %% plot the results
%         h1_sub1 = subplot(4,3,[7,8,10,11]);
        h1_sub1 = subplot(2,2,3);
        axes(h1_sub1);trisurf(original_mesh,'FaceColor','r','FaceAlpha',0,'EdgeAlpha',0.1); daspect([1,1,1]);title(h1_sub1,'Curvature');colormap(h1_sub1,'jet');hold off;
        hold on;trisurf(interpolated_mesh,scaler_interpolated(:,3),'FaceAlpha',1,'EdgeAlpha',0.1);caxis(spont_curv_lipids);colorbar;colormap(jet);hold off
        ylim([-fig_limits,fig_limits]);xlim([-fig_limits,fig_limits]);zlim([-fig_limits,fig_limits]);
        view(0,0);ylim([-fig_limits,0]);xlim([-fig_limits,fig_limits]);zlim([-fig_limits,fig_limits]);hold off;
        
%         h1_sub2 = subplot(4,3,3);
        h1_sub2 = subplot(2,2,2);
        axes(h1_sub2);trisurf(original_mesh,'FaceColor','r','FaceAlpha',0,'EdgeAlpha',0.1); daspect([1,1,1]);title(h1_sub2,'order parameter');colormap(h1_sub2,'jet');hold off;
        hold on;trisurf(interpolated_mesh,scaler_interpolated(:,1),'FaceAlpha',1,'EdgeAlpha',0.1);caxis(spont_curv_lipids);colorbar;colormap(jet);hold off
        ylim([-fig_limits,fig_limits]);xlim([-fig_limits,fig_limits]);zlim([-fig_limits,fig_limits]);hold off;
        
%         h1_sub3 = subplot(4,3,[1,2,4,5]);
        h1_sub3 = subplot(2,2,1);
        axes(h1_sub3);trisurf(original_mesh,'FaceColor','r','FaceAlpha',0,'EdgeAlpha',0.1); daspect([1,1,1]);title(h1_sub3,'composition');colormap(h1_sub3,'jet');hold off;
        hold on;trisurf(interpolated_mesh,scaler_interpolated(:,2),'FaceAlpha',1,'EdgeAlpha',0.1);caxis([0.2,0.8]);colorbar;colormap(jet);hold off
        view(0,0);ylim([-fig_limits,0]);xlim([-fig_limits,fig_limits]);zlim([-fig_limits,fig_limits]);caxis([0.2,0.8]);colorbar;hold off;
        pause(0.1);
        
        [avg_free_energy_at_vertex,mole_lipid_comp,order_parameter,hmean,free_energy_curve,H_space,comp_space]  = ...
                mark_obj_free_energy_plot(obj_list,spont_curv_lipids, kappa, surface_modulous,kappa_coefficient,vanderwaals, kB, temperature);
        
%         h1_sub4 = subplot(4,3,[6,9,12]);
        h1_sub4 = subplot(2,2,4);
        surf(H_space,comp_space,free_energy_curve,'FaceColor','b','FaceAlpha',0.5);
        hold on; scatter3(hmean,mole_lipid_comp(:,1),avg_free_energy_at_vertex,20,'g');hold off ;view(0,0);
        xlabel('curvature');ylabel('mole_fraction');zlabel('avg free energy');
        
        suptitle(sprintf('experiment%d frame no:%d temp: %d radius: %d\n spontcurv [%d %d]',figure_no,t,temperature,radius,spont_curv_lipids(1),spont_curv_lipids(2)))
        fprintf('\n\nLOOP NO %d COMPLETED\n\n',t);
        
        %% Saving into file
        set(h1, 'Position', [100   55  1200  1200]);
        current_frame = getframe(h1);
        imwrite(current_frame.cdata,fullDestinationFileName,'WriteMode','append');
        
    end
    
end
fullDestinationFileName = fullfile(destinationFolder, output_fig_file_name);
savefig(h1,fullDestinationFileName);


end