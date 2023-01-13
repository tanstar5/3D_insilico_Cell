%% Inputs
filename_mesh = 'baumgart.stl';
type_properties = [1*44444,1*44444,1*44444;...
                   1/10,1/50,-1/4;...
                   1,   1,    1];
temperature = 1;
gaussian_modulus = 20;
surface_modulus = 4;
distortion_modulous = 50;
no_particles = 1000;
diffusion_rate = 100;
viscosity = 10;
% kick = 0.012;
% quantal_change = 20;
no_of_loops = 100000;
each_mole_content = 600;

%% Initialize mesh from stl
patch_min_edge = 2;
patch_max_edge = 2;
[membrane_mesh] = import_surface_mesh_from_stl(filename_mesh,patch_min_edge,patch_max_edge);trisurf(membrane_mesh),daspect([1,1,1]);

%% Make the patch list from mesh
[obj_list] = membrane_patch_list(membrane_patch(1),size(membrane_mesh.Points,1));
[obj_list,total_surface_area] = load_spatial_properties_from_mesh(obj_list,membrane_mesh);
[obj_list] = derive_geometrical_quatities_all(obj_list);

[obj_list] = determine_num_particles_per_patch_basedAreaCorrected(obj_list,no_particles,total_surface_area);
[obj_list] = randomly_distribute_lipids(obj_list,[.33,.33,.33]);
% [obj_list] = distribute_lipids_randomly_global(obj_list);


%% check thermodynamics
trisurf(membrane_mesh,'FaceAlpha',.3),daspect([1,1,1]);
% [interesting_ids] = findobj_at_pos(obj_list,[-4.394 0.8875 3.979],.5);
[interesting_ids] = findobj_at_pos(obj_list,[42 16 19],3);

%% Predict simulation run parameters from patch length
[del_t_diffuse,del_t_move,quantal_change,kick,hypothetical_radius] = predict_parameters(obj_list,diffusion_rate,type_properties,surface_modulus,viscosity,temperature);
freq_t_diffuse_call = 1/del_t_diffuse;
freq_t_move_call = 1/del_t_move;
percentage_diffuse_call = freq_t_diffuse_call/(freq_t_diffuse_call+freq_t_move_call)*100;
percentage_move_call = freq_t_move_call/(freq_t_diffuse_call+freq_t_move_call)*100;


indexing_100_calls = [ ones(1,ceil(percentage_diffuse_call)),2*ones(1,ceil(percentage_move_call)) ];

quantal_change = 10
diffusion_rates_ratio = [1,1,1];

for loop_no = 1:no_of_loops
    pick = randperm(length(indexing_100_calls),1);
%     call_move  = indexing_100_calls(pick);
    call_move = 1
    if call_move == 1
        [obj_list] = lipid_exchange_MC_local_move(obj_list,quantal_change,diffusion_rates_ratio, type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous,each_mole_content);
    else
        [obj_list] = vertex_displacement_MC_local_move(obj_list,kick,type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous,quantal_change,each_mole_content);
    end
    
    %% plot surface mesh
    h1 = figure(1);
    h1_sub1 = subplot(1,2,1);
    h1_sub2 = subplot(1,2,2);
    [surface_mesh] = generate_mesh_from_objlist(obj_list); 
    [patch_color,triangle_color_curvature] = color_code_lipid_density_profile_patch_wise_and_curvature_new(obj_list);    
    axes(h1_sub1);trisurf(surface_mesh,triangle_color_curvature,'FaceAlpha',.3,'EdgeAlpha',0.2); daspect([1,1,1]);title(h1_sub1,'Curvature');colormap(h1_sub1,'jet');hold off;
    hold on; quiver3( obj_list(interesting_ids(1)).Pos(1),obj_list(interesting_ids(1)).Pos(2),obj_list(interesting_ids(1)).Pos(3), obj_list(interesting_ids(1)).Normal_vect(1), obj_list(interesting_ids(1)).Normal_vect(2), obj_list(interesting_ids(1)).Normal_vect(3) );hold off;
%     xlim(h1_sub1,x_lim);ylim(h1_sub1,y_lim);zlim(h1_sub1,z_lim);
    xlim(h1_sub1,[-0,70]);ylim(h1_sub1,[-30,30]);zlim(h1_sub1,[-30,30]);
    view(h1_sub1,[90,0]);    
    
    axes(h1_sub2);trisurf(surface_mesh,patch_color(:,3)./(patch_color(:,1)),'FaceAlpha',.3,'EdgeAlpha',0); daspect([1,1,1]);title(h1_sub2,'molefraction');colormap(h1_sub2,'jet');hold off;
    xlim(h1_sub2,[-0,70]);ylim(h1_sub2,[-30,30]);zlim(h1_sub2,[-30,30]);
%     xlim(h1_sub2,x_lim);ylim(h1_sub2,y_lim);zlim(h1_sub2,z_lim);
    view(h1_sub2,[90,0]);
    %% Saving into file    
    set(h1, 'Position', [10 10 900 900]);    
    current_frame = getframe(h1);
    imwrite(current_frame.cdata,'Membrane.tif','WriteMode','append');   
    
    
end