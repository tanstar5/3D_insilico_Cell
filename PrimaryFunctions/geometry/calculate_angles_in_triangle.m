function [angles] = calculate_angles_in_triangle(Points)
edge1_2 = Points(2,:) - Points(1,:);
edge1_2_unit = edge1_2/vecnorm(edge1_2);
edge1_3 = Points(3,:) - Points(1,:);
edge1_3_unit = edge1_3/vecnorm(edge1_3);
angle_1 = acos( dot(edge1_3_unit,edge1_2_unit) );

edge2_1 = Points(1,:) - Points(2,:);
edge2_1_unit = edge2_1/vecnorm(edge2_1);
edge2_3 = Points(3,:) - Points(2,:);
edge2_3_unit = edge2_3/vecnorm(edge2_3);
angle_2 = acos( dot(edge2_3_unit,edge2_1_unit) );

angle_3 = pi - angle_1 - angle_2;
angles = [angle_1,angle_2,angle_3];

end