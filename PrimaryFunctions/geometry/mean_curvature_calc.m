function [mean_curvature] = mean_curvature_calc(membrane_mesh_list)
mean_curvature = NaN(size(membrane_mesh_list,1),1);
parfor obj_id = 1:size(membrane_mesh_list,1)
    curvatures = sort(membrane_mesh_list(obj_id).DarbouxFrame{2,1})
    mean_curvature(obj_id) = (mean(curvatures(2:3)));
end
end