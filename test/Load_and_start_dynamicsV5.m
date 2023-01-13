%% Inputs
filename_mesh = 'cell_multi_19.stl';
type_properties = [1*40000,1*40000,1*40000;...
    1/.5,1/10,-1/6;...
    1,1,1];
temperature = 5;
gaussian_modulus = 20;
surface_modulus = 40;
distortion_modulous = 50;
no_particles = 1000;
kick = 0.012;
quantal_change = 20;
no_of_loops = 100000;
each_mole_content = 600;



%% Initialize mesh from stl
patch_min_edge = 0.6;
patch_max_edge = 0.6;
[membrane_mesh] = import_surface_mesh_from_stl(filename_mesh,patch_min_edge,patch_max_edge);trisurf(membrane_mesh),daspect([1,1,1]);

%% Make the patch list from mesh
[obj_list] = membrane_patch_list(membrane_patch(1),size(membrane_mesh.Points,1));
[obj_list] = load_spatial_properties_from_mesh(obj_list,membrane_mesh);
[obj_list] = derive_geometrical_quatities_all(obj_list);

[obj_list] = determine_num_particles_per_patch_basedArea(obj_list,no_particles);
[obj_list] = randomly_distribute_lipids(obj_list,[.33,.33,.33]);
% [obj_list] = distribute_lipids_randomly_global(obj_list);

%% check thermodynamics
trisurf(membrane_mesh,'FaceAlpha',.3),daspect([1,1,1]);
[az,el] = view;
[x_lim] = xlim;
[y_lim] = ylim;
[z_lim] = zlim;
% [interesting_ids] = findobj_at_pos(obj_list,[-4.394 0.8875 3.979],.5);
[interesting_ids] = findobj_at_pos(obj_list,[-0.7015 10.59 3.292],.5);

[c_preferred,c_current,preffered_mole_fractions,actual_mole_fractions] = ...
                calculate_equillibrium_characteristics(obj_list(interesting_ids(1)),type_properties,temperature,1,each_mole_content)

%% Loop
curvature_energy_list = NaN(no_of_loops,1);
entropy_of_mixing_list = NaN(no_of_loops,1);
surface_stretching_energy_list = NaN(no_of_loops,1);
distortion_energy_list = NaN(no_of_loops,1);
total_energy = NaN(no_of_loops,1);
H_mean_list = NaN(no_of_loops,1);
H_spontaneous_list = NaN(no_of_loops,1);
c1_current = NaN(no_of_loops,1);
c2_current = NaN(no_of_loops,1);
c_preferred = NaN(no_of_loops,1);
c_current = NaN(no_of_loops,1);
preffered_mole_fractions = NaN(no_of_loops,3);
actual_mole_fractions = NaN(no_of_loops,3);

%% prediction

% [obj_list] = vertex_displacement_MC_local_move(obj_list,0*0.3*(obj_list(1).Av_vertex).^.5,type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous,quantal_change,each_mole_content);
% [obj_list] = vertex_displacement_MC_local_move(obj_list,0*kick,type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous,quantal_change,each_mole_content);
for loop_no = 1:no_of_loops
    
    fprintf('\n LOOP NO = %d STARTED \n',loop_no);
    %% MC Moves
    if mod(loop_no,1)==0
         
    end
    [obj_list] = lipid_exchange_MC_local_move(obj_list,quantal_change, type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous,each_mole_content);
    
    %% plot surface mesh
    h1 = figure(1);
    h1_sub1 = subplot(1,2,1);
    h1_sub2 = subplot(1,2,2);
    [surface_mesh] = generate_mesh_from_objlist(obj_list); 
    [patch_color,triangle_color_curvature] = color_code_lipid_density_profile_patch_wise_and_curvature(obj_list);    
    axes(h1_sub1);trisurf(surface_mesh,triangle_color_curvature,'FaceAlpha',.3,'EdgeAlpha',0.2); daspect([1,1,1]);title(h1_sub1,'Curvature');colormap(h1_sub1,'jet');hold off;
    hold on; quiver3( obj_list(interesting_ids(1)).Pos(1),obj_list(interesting_ids(1)).Pos(2),obj_list(interesting_ids(1)).Pos(3), obj_list(interesting_ids(1)).Normal_vect(1), obj_list(interesting_ids(1)).Normal_vect(2), obj_list(interesting_ids(1)).Normal_vect(3) );hold off;
    xlim(h1_sub1,x_lim);ylim(h1_sub1,y_lim);zlim(h1_sub1,z_lim);
    view(h1_sub1,[az,el]);    
    
    axes(h1_sub2);trisurf(surface_mesh,patch_color(:,3)./sum(patch_color(:,:),2),'FaceAlpha',.3,'EdgeAlpha',0); daspect([1,1,1]);title(h1_sub2,'molefraction');caxis([0,.4]);colormap(h1_sub2,'jet');hold off;
%     xlim(h1_sub2,[-10,15]);ylim(h1_sub2,[-12,12]);zlim(h1_sub2,[-15,0]);
    xlim(h1_sub2,x_lim);ylim(h1_sub2,y_lim);zlim(h1_sub2,z_lim);
    view(h1_sub2,[az,el]);
    %% Saving into file    
    set(h1, 'Position', [10 10 900 900]);    
    current_frame = getframe(h1);
    imwrite(current_frame.cdata,'Membrane.tif','WriteMode','append');
    
    % tracking a patch
%     [c_preferred,c_expected,preffered_mole_fractions,actual_mole_fractions] = ...
%         calculate_equillibrium_characteristics(obj,type_properties,temperature,plot_profile);
    
    [patch_pos_list,curvature_energy_list(loop_no),entropy_of_mixing_list(loop_no),surface_stretching_energy_list(loop_no),distortion_energy_list(loop_no),H_mean_list(loop_no),H_spontaneous_list(loop_no)] = ...
                track_patch(obj_list(interesting_ids(1)),type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous,quantal_change,each_mole_content);
    total_energy(loop_no) =  curvature_energy_list(loop_no)+ surface_stretching_energy_list(loop_no)+ distortion_energy_list(loop_no);
    
    h2 = figure(2);
    h2_sub1 = subplot(2,2,1);
    time_axis = 1:loop_no;
    plot(h2_sub1,time_axis,curvature_energy_list(time_axis),'-*b',time_axis,entropy_of_mixing_list(time_axis),'-*m',...
        time_axis,surface_stretching_energy_list(time_axis),'-*g',time_axis,distortion_energy_list(time_axis),'-*c',...
        time_axis,total_energy(time_axis),'-r');
    h2_sub2 = subplot(2,2,2);
    plot(h2_sub2,1:3,obj_list(interesting_ids(1)).lipid_ratio_up);
    h2_sub4 = subplot(2,2,4);
    plot(H_mean_list,H_spontaneous_list,'-*');
    h2_sub3 = subplot(2,2,3);
    if loop_no >=2
        plot(h2_sub3,2:loop_no,diff(curvature_energy_list(1:loop_no)),'-*b',2:loop_no,diff(entropy_of_mixing_list(1:loop_no)),'-*m',...
        2:loop_no,diff(surface_stretching_energy_list(1:loop_no)),'-*g',2:loop_no,diff(distortion_energy_list(1:loop_no)),'-*c',...
        2:loop_no,diff(total_energy(1:loop_no)),'-r');
    end
    
    
    %% Making RGB channel for three lipid types
    h3 = figure(3);
    h3_sub = subplot(1,1,1);
    trisurf(surface_mesh,patch_color(:,1)./sum(patch_color,2),'FaceAlpha',.33,'EdgeAlpha',0); daspect([1,1,1]);view(h3_sub,[-90,0]);colormap(h3_sub,'gray');
    set(h3, 'Position', [10 10 450 450]);
    current_frame = getframe(h3);
    imwrite(current_frame.cdata,'Membrane_RED.tif','WriteMode','append');
    
    trisurf(surface_mesh,patch_color(:,2)./sum(patch_color,2),'FaceAlpha',.33,'EdgeAlpha',0); daspect([1,1,1]);view(h3_sub,[-90,0]);colormap(h3_sub,'gray');
    set(h3, 'Position', [10 10 450 450]);
    current_frame = getframe(h3);
    imwrite(current_frame.cdata,'Membrane_GREEN.tif','WriteMode','append');
    
    trisurf(surface_mesh,patch_color(:,3)./sum(patch_color,2),'FaceAlpha',.33,'EdgeAlpha',0); daspect([1,1,1]);view(h3_sub,[-90,0]);colormap(h3_sub,'gray');
    set(h3, 'Position', [10 10 450 450]);
    current_frame = getframe(h3);
    imwrite(current_frame.cdata,'Membrane_BLUE.tif','WriteMode','append');
    
    h4 = figure(4);
    [c_preferred(loop_no),c_current(loop_no),preffered_mole_fractions(loop_no,:),actual_mole_fractions(loop_no,:),c1_current(loop_no),c2_current(loop_no)] = calculate_equillibrium_characteristics(obj_list(interesting_ids(1)),type_properties,temperature,1,each_mole_content);
    c_preferred(loop_no) = sum(obj_list(interesting_ids(1)).lipid_ratio_up.*type_properties(2,:))/sum(obj_list(interesting_ids(1)).lipid_ratio_up);
    h5 = figure(5);
    h5_sub1 = subplot(2,1,1);
    plot(h5_sub1,time_axis,c_preferred(time_axis),'-r',time_axis,c_current(time_axis),'-*r', time_axis,c1_current(time_axis),'-*g',time_axis,c2_current(time_axis),'-*b');
    h5_sub2 = subplot(2,1,2);
    plot(h5_sub2,time_axis,preffered_mole_fractions(time_axis,1),'-+r',time_axis,preffered_mole_fractions(time_axis,2),'-or',time_axis,preffered_mole_fractions(time_axis,3),'-xr',...
        time_axis,actual_mole_fractions(time_axis,1),'-+b',time_axis,actual_mole_fractions(time_axis,2),'-ob',time_axis,actual_mole_fractions(time_axis,3),'-xb');
    
    pause(0.1)
end