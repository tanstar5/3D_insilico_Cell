function [points_transformed] = transform_coor_to_plane(Points_to_transform,BasisVectors,origin_coor)
unit_basis_vect1 = BasisVectors(1,:);
unit_basis_vect2 = BasisVectors(2,:);
unit_basis_vect3 = BasisVectors(3,:);

points_transformed = NaN(size(Points_to_transform));
for pt = 1:size(Points_to_transform,1)
%     P = Points_to_transform(pt,:);
    points_transformed(pt,1) = dot(Points_to_transform(pt,:) - origin_coor,unit_basis_vect1);
    points_transformed(pt,2) = dot(Points_to_transform(pt,:) - origin_coor,unit_basis_vect2);
    points_transformed(pt,3) = dot(Points_to_transform(pt,:) - origin_coor,unit_basis_vect3);
end

end