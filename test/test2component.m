%% File to save output operations
output_file_name = 'penalty_abs_change_least_square_0p05_rad_300_realistic.tif';
destinationFolder = 'C:\Users\tanum\OneDrive\Documents\PhD work\workfolder_surface\Self_Organization_project\active_cell_version41_two_component_gradient_term\outputs';

%%
if ~exist(destinationFolder, 'dir')
  mkdir(destinationFolder);
end
fullDestinationFileName = fullfile(destinationFolder, output_file_name);
%%
radius = 300;
time_steps = 5000;
step_size = 0.3;
particle_packet_size = 30;


lipid_head_area = 0.5; % in nm^2
kB = 1;
temperature = 100;
kappa = 60*kB*300;
surface_modulous = kappa/4*12; %% relation between kappa and surface modulous is kappa = surface_molous*t^2 where t is bilayer thickness (wikipedia: search Lipid bilayer mechanics)
                            % or in popular literature kappa =
                            % surface_modulous*t^2/12
penalty_factor = 0.05*kB*temperature;                            
spont_curv_lipids = [-0.2+1/radius,1/radius+0.2];

%% Initilize membrane from mesh
filename_mesh = 'circle_rad_20.stl';
fig_limits = radius + 10;
patch_min_edge = 1.4;  %3;
patch_max_edge = 1.4;  %3;
[membrane_mesh] = import_surface_mesh_from_stl(filename_mesh,patch_min_edge,patch_max_edge);trisurf(membrane_mesh),daspect([1,1,1]);

%% Due to inconsistency in stl file we scale it keepin the number of points same based on radius
membrane_mesh_points = membrane_mesh.Points;
membrane_mesh_points_scaled = membrane_mesh_points/20*radius;
membrane_mesh = triangulation(membrane_mesh.ConnectivityList,membrane_mesh_points_scaled);
patch_min_edge_scaled  = patch_min_edge/20*radius
patch_max_edge_scaled  = patch_max_edge/20*radius

obj = membrane2D_continuous(1,[],[],[],[]);
[obj_list] = initilize_obj_list_container_from_mesh(obj,membrane_mesh);
[obj_list] = populate_with_lipids_at_equillibrium(obj_list,lipid_head_area,spont_curv_lipids);

% check geometric calculations
[obj_list] = derive_geometrical_quatities_all(obj_list);
[normal_vects_all] = extract_Normal_vects(obj_list);
[mean_curvature,principle_curvatures] = extract_curvature_list(obj_list);
[triangle_color_curvature] = color_code_patch_wise_curvature(obj_list);
[surface_mesh] = generate_mesh_from_objlist(obj_list);
h1_sub1 = subplot(2,2,1);
axes(h1_sub1);trisurf(surface_mesh,triangle_color_curvature,'FaceAlpha',.3,'EdgeAlpha',0.2); daspect([1,1,1]);title(h1_sub1,'Curvature');colormap(h1_sub1,'jet');hold off;
h1_sub2 = subplot(2,2,2);
histogram(h1_sub2,triangle_color_curvature)

[triangle_color_mole_fraction] = color_code_patch_wise_composition(obj_list);
h1_sub3 = subplot(2,2,3);
axes(h1_sub3);trisurf(surface_mesh,triangle_color_mole_fraction(:,2),'FaceAlpha',.3,'EdgeAlpha',0.2); daspect([1,1,1]);title(h1_sub3,'Composition');colormap(h1_sub3,'jet');hold off;

%% Evolve
h1 = figure(1);
for t = 1:time_steps   
    [obj_list] = orthogonal_move(obj_list,step_size,spont_curv_lipids, kappa, surface_modulous, kB, temperature);
    [obj_list] = exchange_move(obj_list,particle_packet_size,spont_curv_lipids, kappa, surface_modulous,penalty_factor, kB, temperature);
    [coors,principle_curvatures,hmean,lipid_comp,mole_lipid_comp,order_parameter]  = extract_patch_properties(obj_list,spont_curv_lipids);
    
    %% interpolating the scalers in an interpolated mesh 
    all_scaler_mat = [ order_parameter,mole_lipid_comp(:,2),hmean];
    [scaler_interpolated,interpolated_mesh,original_mesh] = interpolate_vertex_scalers_to_mesh(obj_list,all_scaler_mat);
    
    %% plot the results
    h1_sub1 = subplot(2,2,2);
    axes(h1_sub1);trisurf(original_mesh,'FaceColor','r','FaceAlpha',0,'EdgeAlpha',0.1); daspect([1,1,1]);title(h1_sub1,'Curvature');colormap(h1_sub1,'parula');hold off;
    hold on;trisurf(interpolated_mesh,scaler_interpolated(:,3),'FaceAlpha',1,'EdgeAlpha',0.1);caxis(spont_curv_lipids);colorbar;colormap(parula);hold off
    ylim([-fig_limits,fig_limits]);xlim([-fig_limits,fig_limits]);zlim([-fig_limits,fig_limits]);hold off;
    
    h1_sub2 = subplot(2,2,4);
    axes(h1_sub2);trisurf(original_mesh,'FaceColor','r','FaceAlpha',0,'EdgeAlpha',0.1); daspect([1,1,1]);title(h1_sub2,'order parameter');colormap(h1_sub2,'parula');hold off;
    hold on;trisurf(interpolated_mesh,scaler_interpolated(:,1),'FaceAlpha',1,'EdgeAlpha',0.1);caxis(spont_curv_lipids);colorbar;colormap(parula);hold off
    ylim([-fig_limits,fig_limits]);xlim([-fig_limits,fig_limits]);zlim([-fig_limits,fig_limits]);hold off;
    
    h1_sub3 = subplot(2,2,[1,3]);
    axes(h1_sub3);trisurf(original_mesh,'FaceColor','r','FaceAlpha',0,'EdgeAlpha',0.1); daspect([1,1,1]);title(h1_sub3,'composition');colormap(h1_sub3,'parula');hold off;
    hold on;trisurf(interpolated_mesh,scaler_interpolated(:,2),'FaceAlpha',1,'EdgeAlpha',0.1);caxis([0.2,0.8]);colorbar;colormap(parula);hold off
    view(0,0);ylim([-fig_limits,0]);xlim([-fig_limits,fig_limits]);zlim([-fig_limits,fig_limits]);caxis([0.2,0.8]);colorbar;hold off;
    pause(0.1);   
    
% %     [triangle_color_curvature] = color_code_patch_wise_curvature(obj_list);
%     [surface_mesh] = generate_mesh_from_objlist(obj_list);
%     h1_sub1 = subplot(2,2,2);
%     axes(h1_sub1);trisurf(surface_mesh,'FaceColor','r','FaceAlpha',.3,'EdgeAlpha',0.2); daspect([1,1,1]);title(h1_sub1,'Curvature');colormap(h1_sub1,'jet');hold off;
%     hold on; scatter3( coors(:,1),coors(:,2),coors(:,3),5,hmean,'filled' );caxis([0.02-0.1,0.02+0.1]);colorbar;
%     ylim([-fig_limits,fig_limits]);xlim([-fig_limits,fig_limits]);zlim([-fig_limits,fig_limits]);hold off;
%     
% %     [triangle_color_mole_fraction] = color_code_patch_wise_composition(obj_list);
%     h1_sub2 = subplot(2,2,4);
%     axes(h1_sub2);trisurf(surface_mesh,'FaceColor','r','FaceAlpha',.3,'EdgeAlpha',0.2); daspect([1,1,1]);title(h1_sub2,'Composition');colormap(h1_sub2,'jet');hold off;
%     hold on; scatter3( coors(:,1),coors(:,2),coors(:,3),5,mole_lipid_comp(:,2),'filled' );caxis([0.2,0.8]);colorbar;
%     ylim([-fig_limits,fig_limits]);xlim([-fig_limits,fig_limits]);zlim([-fig_limits,fig_limits]);hold off;
%     
%     
%     h1_sub3 = subplot(2,2,[1,3]);
% %     h1_sub3 = subplot(1,1,1);
% %     order_parameter = triangle_color_mole_fraction(:,1)*spont_curv_lipids(1) + triangle_color_mole_fraction(:,2)*spont_curv_lipids(2); 
%     axes(h1_sub3);trisurf(surface_mesh,'FaceColor','r','FaceAlpha',.3,'EdgeAlpha',0.2); daspect([1,1,1]);title(h1_sub3,'Order Parameter');colormap(h1_sub3,'jet');
%     hold on; scatter3( coors(:,1),coors(:,2),coors(:,3),8,order_parameter,'filled' );hold off;
%     view(0,0);ylim([-fig_limits,0]);xlim([-fig_limits,fig_limits]);zlim([-fig_limits,fig_limits]);caxis([0.02-0.2,0.02+0.2]);colorbar;colormap(hot);hold off;
%     pause(0.1);
    fprintf('\n\nLOOP NO %d COMPLETED\n\n',t);
    
    %% Saving into file    
    set(h1, 'Position', [10 10 900 450]);    
    current_frame = getframe(h1);
    imwrite(current_frame.cdata,fullDestinationFileName,'WriteMode','append');
    
end


%% try
% updated_obj_list=orthogonal_move(obj_list,.1,[-0.2,0.2], 60*300, 10, 1, 300);
% [triangle_color_curvature] = color_code_patch_wise_curvature(updated_obj_list);
% [surface_mesh] = generate_mesh_from_objlist(updated_obj_list);
% h1_sub1 = subplot(1,1,1);
% axes(h1_sub1);trisurf(surface_mesh,triangle_color_curvature,'FaceAlpha',.3,'EdgeAlpha',0.2); daspect([1,1,1]);title(h1_sub1,'Curvature');colormap(h1_sub1,'jet');hold off;
% 
% updated_obj_list=exchange_move(obj_list,10,[-0.2,0.2], 200, 10, 1, .1);
% [triangle_color_mole_fraction] = color_code_patch_wise_composition(updated_obj_list);
% h1_sub3 = subplot(1,1,1);
% axes(h1_sub3);trisurf(surface_mesh,triangle_color_mole_fraction(:,2),'FaceAlpha',.3,'EdgeAlpha',0.2); daspect([1,1,1]);title(h1_sub3,'Composition');colormap(h1_sub3,'jet');hold off;
