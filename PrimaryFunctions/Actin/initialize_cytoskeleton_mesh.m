function [cytoskeleton_mesh_list] = initialize_cytoskeleton_mesh(barbed_end_mesh)
    barbed_end_points = barbed_end_mesh.Points;
    barbed_end_triangles = barbed_end_mesh.ConnectivityList;
    cytoskeleton_mesh_list(size(barbed_end_triangles,1),1) = cytoskeleton_mesh(1);
    N_IDs = neighbors(barbed_end_mesh);
    for tr = 1:size(barbed_end_triangles,1)
        cytoskeleton_mesh_obj = cytoskeleton_mesh(tr);
        cytoskeleton_mesh_obj = assign_points(cytoskeleton_mesh_obj, barbed_end_points(barbed_end_triangles(tr,:),:)   );
        cytoskeleton_mesh_obj = assign_triangle(cytoskeleton_mesh_obj,barbed_end_triangles(tr,:));
        cytoskeleton_mesh_obj = assign_neighbours(cytoskeleton_mesh_obj,N_IDs(tr,:));
        cytoskeleton_mesh_obj = assign_neighbour_triangles(cytoskeleton_mesh_obj,barbed_end_triangles(N_IDs(tr,:),:));
        cytoskeleton_mesh_obj = assign_edges(cytoskeleton_mesh_obj,barbed_end_triangles(tr,:));
        cytoskeleton_mesh_list(tr) = cytoskeleton_mesh_obj;
    end
end