all_Points = membrane_mesh.Points;
for tr = 1:size(membrane_mesh.ConnectivityList,1)
    tri_id = membrane_mesh.ConnectivityList(tr,:);
    Points = all_Points(tri_id,:);
    area_triangle(Points)
end


% area_triangle(Points)