function [points_transformed] = inverse_transform_to_XYZ_coor(points,BasisVectors,origin)
%% This function converts the Points (Pei Pej Pek) where Pei Pej Pek are column vectors
%% representing the coordinates of the points in the old frame represented by the basis 
%% represented by the basis vectors ei, ej, ek (row vectors) represented in x (1,0,0) y (0,1,0) z (0,0,1)
%% coordinate system
transformation_mat = BasisVectors';
points_transformed = NaN(size(points));
parfor pts = 1:size(points,1)
    transformed_coor = transformation_mat*points(pts,:)';
    points_transformed(pts,:) = transformed_coor' + origin;    
end
end