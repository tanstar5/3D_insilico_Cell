%% Inputs
filename_mesh = 'circle.stl';
type_properties = [1*40,1*40,1*40;...
                   1/20,1/50,-1/20;...
                   1,   1,    1];
temperature = 1;
gaussian_modulus = 0;
surface_modulus = 40;
distortion_modulous = 0;
no_particles = 5000;
diffusion_rate = 100;
viscosity = 10;
% kick = 0.012;
% quantal_change = 20;
no_of_loops = 2800;
each_mole_content = 0.01;

%% Initialize mesh from stl
patch_min_edge = 3.5;
patch_max_edge = 4;
[membrane_mesh] = import_surface_mesh_from_stl(filename_mesh,patch_min_edge,patch_max_edge);trisurf(membrane_mesh),daspect([1,1,1]);

%% Make the patch list from mesh
[obj_list] = membrane_patch_list(membrane_patch(1),size(membrane_mesh.Points,1));
[obj_list,total_surface_area] = load_spatial_properties_from_mesh(obj_list,membrane_mesh);
[obj_list] = derive_geometrical_quatities_all(obj_list);

[obj_list] = determine_num_particles_per_patch_basedAreaCorrected(obj_list,no_particles,total_surface_area,temperature);
[obj_list] = randomly_distribute_lipids(obj_list,[0.33,0.33,0.33]);
% [obj_list] = populate_patches_with_lipid_type(obj_list,[0 ,0,50],10,[0.70,0.15,0.15])
% [obj_list] = populate_patches_with_lipid_type(obj_list,[46.55,9.807,-1.3],5,[0.15,0.15,0.70])
% [obj_list] = populate_patches_with_lipid_type(obj_list,[-2,55,0],5,[0.15,0.15,0.70])
[obj_list] = distribute_lipids_smoothed_over_multiple_nodes(obj_list,5, [0.65,0.30,0.05], 6);
[obj_list] = distribute_lipids_smoothed_over_multiple_nodes(obj_list,4, [0.05,0.30,0.65], 5);
%% check thermodynamics
trisurf(membrane_mesh,'FaceAlpha',.3),daspect([1,1,1]);


%% Simulation parameters
quantal_change = 40;
diffusion_rates_ratio = [1,1,1];
kick = patch_min_edge/15;
percentage_diffuse_call = 50;
percentage_move_call = 50;


indexing_100_calls = [ ones(1,ceil(percentage_diffuse_call)),2*ones(1,ceil(percentage_move_call)) ];

for loop_no = 1:no_of_loops
    pick = randperm(length(indexing_100_calls),1);    
    
    %% plot surface mesh
    h1 = figure(1);
    h1_sub1 = subplot(1,2,1);
    h1_sub2 = subplot(1,2,2);
    [surface_mesh] = generate_mesh_from_objlist(obj_list); 
    [patch_color,triangle_color_curvature] = color_code_lipid_density_profile_patch_wise_and_curvature_new(obj_list);    
    axes(h1_sub1);trisurf(surface_mesh,triangle_color_curvature,'FaceAlpha',.3,'EdgeAlpha',0.2); daspect([1,1,1]);title(h1_sub1,'Curvature');colormap(h1_sub1,'jet');hold off;
%     axes(h1_sub1);trisurf(surface_mesh,patch_color(:,1)+patch_color(:,2),'FaceAlpha',.3,'EdgeAlpha',0.2); daspect([1,1,1]);title(h1_sub1,'molefraction typePositive');caxis([0.4,0.6]);colormap(h1_sub1,'jet');hold off;
    xlim(h1_sub1,[-0,70]);ylim(h1_sub1,[-70,70]);zlim(h1_sub1,[-70,70]);
    view(h1_sub1,[90,0]);    
    
    axes(h1_sub2);trisurf(surface_mesh,patch_color(:,1),'FaceAlpha',.3,'EdgeAlpha',0.2); daspect([1,1,1]);title(h1_sub2,'molefraction HighlyPositive');caxis([0.20,0.8]);colormap(h1_sub2,'jet');hold off;
    xlim(h1_sub2,[-0,70]);ylim(h1_sub2,[-70,70]);zlim(h1_sub2,[-70,70]);
%     xlim(h1_sub2,x_lim);ylim(h1_sub2,y_lim);zlim(h1_sub2,z_lim);
    view(h1_sub2,[90,0]);
    %% Saving into file    
    set(h1, 'Position', [10 10 900 450]);    
    current_frame = getframe(h1);
    imwrite(current_frame.cdata,'Membrane.tif','WriteMode','append');   
    pause(0.01);
    call_move  = indexing_100_calls(pick);
    
%     if call_move == 1
        [obj_list] = lipid_exchange_MC_local_move(obj_list,quantal_change,diffusion_rates_ratio, type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous,each_mole_content);
%     else
        [obj_list] = vertex_displacement_MC_local_move(obj_list,kick,type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous,quantal_change,each_mole_content);
%     end
    
end