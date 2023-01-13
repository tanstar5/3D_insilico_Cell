no_of_loops = 100000;
for loop_no = 1:no_of_loops
    
    fprintf('\n LOOP NO = %d STARTED \n',loop_no);
    
    [obj_list] = vertex_displacement_MC_local_move(obj_list,0*0.2*(obj_list(1).Av_vertex).^.5,type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous,quantal_change);
    [obj_list] = lipid_exchange_MC_local_move(obj_list,quantal_change, type_properties,temperature,gaussian_modulus,surface_modulus,distortion_modulous);
    
    [surface_mesh] = generate_mesh_from_objlist(obj_list);
%     [patch_color] = color_code_lipid_density_profile_patch_wise(obj_list);
    
    [~,triangle_color_curvature] = color_code_lipid_density_profile_patch_wise_and_curvature(obj_list);   
    trisurf(surface_mesh,triangle_color_curvature,'FaceAlpha',.2); daspect([1,1,1])
    
    [color_patch] = color_code_patch_mole_fraction_for_scatter_plot(obj_list);
    hold on; scatter3(surface_mesh.Points(:,1),surface_mesh.Points(:,2),surface_mesh.Points(:,3),20,color_patch(:,1)./color_patch(:,3),'filled'); hold off
    
    pause(0.1)
end
