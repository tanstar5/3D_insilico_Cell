h1_sub3 = subplot(2,2,[1,3]);
    
    [coors,principle_curvatures,hmean,lipid_comp,mole_lipid_comp,order_parameter]  = extract_patch_properties(obj_list,spont_curv_lipids);
%     order_parameter = triangle_color_mole_fraction(:,1)*spont_curv_lipids(1) + triangle_color_mole_fraction(:,2)*spont_curv_lipids(2); 
    
    fprintf('\n\nLOOP NO %d COMPLETED\n\n',t);
    
    [lipid_comp] = extract_lipids_composition_ID_specific(obj_list);
    obj = obj_list(350);
    obj.neighbours_lipid_compositions_up = lipid_comp(obj.neighbours,:);
    
    [grad_vect,grad_vect_along_surf] = calculate_gradient_least_square_method(obj);
   
    h1_sub3 = subplot(1,1,1);
    hold on;trisurf(interpolated_mesh,scaler_interpolated(:,2),'FaceAlpha',1,'EdgeAlpha',0.1);caxis([0.2,0.8]);colorbar;colormap(parula);daspect([1,1,1]);hold off
    view(0,0);ylim([-300,0]);xlim([-300,300]);zlim([-300,300]);caxis([0.2,0.8]);colorbar;hold off;
    pause(0.1);
    
    hold on;scatter3(obj.Pos(1),obj.Pos(2),obj.Pos(3),100,'r','filled');hold off
    grant_along_plane = grad_vect - dot(obj.Normal_vect(:),grad_vect(:))*(obj.Normal_vect)';
    hold on; quiver3(obj.Pos(1),obj.Pos(2),obj.Pos(3),grad_vect_along_surf(1),grad_vect_along_surf(2),grad_vect_along_surf(3),20,'LineWidth',3);colormap(parula);hold off
%     hold on; quiver3(obj.Pos(1),obj.Pos(2),obj.Pos(3),grant_along_plane(1),grant_along_plane(2),grant_along_plane(3),10,'LineWidth',10);hold off
    
    