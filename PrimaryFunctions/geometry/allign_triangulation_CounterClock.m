function [alligned_triangulation] = allign_triangulation_CounterClock(triangulation_toallign,reference_point)
    Points_triangulation = triangulation_toallign.Points;
    triangulation_mat = triangulation_toallign.ConnectivityList;
    face_normals = faceNormal(triangulation_toallign);
    triangulation_mat_new = NaN(size(triangulation_mat));
    for tri_id = 1:size(triangulation_mat,1)
        tri_current = triangulation_mat(tri_id,:);
        face_normal_current = face_normals(tri_id,:);
        
        tri_points_current = Points_triangulation(tri_current,:);
        centroid_current = mean(tri_points_current,1);
        
        outward_vect = centroid_current - reference_point;
        
        allignment = sign( dot(outward_vect,face_normal_current) );
        
        if allignment<=0
             triangulation_mat_new(tri_id,:) = fliplr( triangulation_mat(tri_id,:) );
             fprintf('    >> Triangulation ANTICLOCKWISED\n');
        else
             triangulation_mat_new(tri_id,:) = triangulation_mat(tri_id,:) ;
             fprintf('    >> Triangulation CORRECT\n');
        end    
    end
    alligned_triangulation = triangulation(triangulation_mat_new,Points_triangulation);
    
end