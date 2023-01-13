no_of_loops = 1000
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
    
    axes(h1_sub2);trisurf(surface_mesh,patch_color(:,1)+patch_color(:,2),'FaceAlpha',.3,'EdgeAlpha',0.2); daspect([1,1,1]);title(h1_sub2,'molefraction HighlyPositive');caxis([0.20,0.8]);colormap(h1_sub2,'jet');hold off;
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