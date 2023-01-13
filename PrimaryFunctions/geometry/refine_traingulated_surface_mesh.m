function [refined_surface_mesh] = refine_traingulated_surface_mesh(surface_mesh)
elements_mat = surface_mesh.ConnectivityList;
points_all = surface_mesh.Points;

refined_mesh_elements = cell(size(elements_mat,1),1);

for tr = 1:size(elements_mat,1)
    element = elements_mat(:,tr);
    nodes_real = points_all(element,:);
    
    P1 = nodes_real(1,:);
    P2 = nodes_real(2,:);
    P3 = nodes_real(3,:);
    [output] = equation_of_plane(P1,P2,P3);
    nodes_transformed = output.transformed_coor;
    nodes_transformed_2D = nodes_transformed(:,1:2);
    
    model = createpde();
    geometryFromMesh(model,nodes_transformed_2D',element');
    mesh_default = generateMesh(model);
    
    nodes_new_transformed2D = mesh_default.Nodes;
    nodes_new_transformed2D = nodes_new_transformed2D';  
    
    
    
end
model = createpde();
% geometryFromMesh(model,nodes,elements);
end