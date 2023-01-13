%% Initialize mesh from stl
patch_min_edge = 0.55;
patch_max_edge = 0.55;
[membrane_mesh] = import_surface_mesh_from_stl('cell_tether.stl',patch_min_edge,patch_max_edge);trisurf(membrane_mesh),daspect([1,1,1]);

[obj_list] = membrane_patch_list(membrane_patch(1),size(membrane_mesh.Points,1));
[obj_list] = load_spatial_properties_from_mesh(obj_list,membrane_mesh);
[obj_list] = derive_geometrical_quatities_all(obj_list);

[obj_list] = determine_num_particles_per_patch_basedArea(obj_list,50);
[obj_list] = randomly_distribute_lipids(obj_list,[.33,.33,.33]);
% [obj_list] = distribute_lipids_randomly_global(obj_list);
% type_properties = [20000,20000,20000;... %% curvature based accumulation
%                    1/2,0,-1/10;...
%                    1,1,1];

type_properties = [20*650,20*650,20*650;...
    1/8,0,-1/8;...
    1,1,1];
temperature = 6;
gaussian_modulus = 1;
surface_modulus = 20;
distortion_modulous = 5;
quantal_change = 10;



no_of_loops = 100000;
curvature_energy_list = NaN(no_of_loops,1);
entropy_of_mixing_list = NaN(no_of_loops,1);
surface_stretching_energy_list = NaN(no_of_loops,1);
distortion_energy_list = NaN(no_of_loops,1);
H_mean_list = NaN(no_of_loops,1);
H_spontaneous_list = NaN(no_of_loops,1);

% [obj_list] = vertex_displacement_MC_local_move(obj_list,0*0.3*(obj_list(1).Av_vertex).^.5,type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous,quantal_change);
for loop_no = 1:no_of_loops
    
    fprintf('\n LOOP NO = %d STARTED \n',loop_no);
%     if mod(loop_no,1) == 0
        [obj_list] = vertex_displacement_MC_local_move(obj_list,1*0.3*(obj_list(1).Av_vertex).^.5,type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous,quantal_change);
%     end
    if mod(loop_no,1) == 0
        [obj_list] = lipid_exchange_MC_local_move(obj_list,quantal_change, type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous);
    end
    [surface_mesh] = generate_mesh_from_objlist(obj_list);    
    
    [patch_color,triangle_color_curvature] = color_code_lipid_density_profile_patch_wise_and_curvature(obj_list);
    h1 = figure(1);
    if mod(loop_no,1) == 0
%         trisurf(surface_mesh,triangle_color_curvature,'FaceAlpha',.2); daspect([1,1,1]);hold off 
%     else
        trisurf(surface_mesh,patch_color(:,3)./patch_color(:,1),'FaceAlpha',.2); daspect([1,1,1]);hold off
%         trisurf(surface_mesh,patch_color(:,1)./patch_color(:,3),'FaceAlpha',.2); daspect([1,1,1]);hold off
    end
    
%     if mod(loop_no,1) == 0
%         [color_patch] = color_code_patch_mole_fraction_for_scatter_plot(obj_list);
%         hold on; scatter3(surface_mesh.Points(:,1),surface_mesh.Points(:,2),surface_mesh.Points(:,3),40,color_patch(:,1)./color_patch(:,3),'filled');caxis([0 1.5]) ;hold off
%     end
    view([0,-90]);
    %% Saving into file
    
    set(h1, 'Position', [10 10 900 900]);
    xlim([-10,15]);ylim([-12,12]);zlim([-15,0]);
    current_frame = getframe(h1);
    imwrite(current_frame.cdata,'Membrane.tif','WriteMode','append');
    
    % tracking a patch
    
    [patch_pos_list,curvature_energy_list(loop_no),entropy_of_mixing_list(loop_no),surface_stretching_energy_list(loop_no),distortion_energy_list(loop_no),H_mean_list(loop_no),H_spontaneous_list(loop_no)] = ...
                track_patch(obj_list(1),type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous,quantal_change);
            
    h2 = figure(2);
    h2_sub1 = subplot(2,2,1);
    time_axis = 1:loop_no;
    plot(h2_sub1,time_axis,curvature_energy_list(time_axis),'-*b',time_axis,entropy_of_mixing_list(time_axis),'-*r');
    h2_sub2 = subplot(2,2,2);
    plot(h2_sub2,1:3,obj_list(1).lipid_ratio_up);
    h3_sub2 = subplot(2,2,[3,4]);
    plot(H_mean_list,H_spontaneous_list,'-*');
    pause(0.1)
    
end