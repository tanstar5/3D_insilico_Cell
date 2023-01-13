function [scales_triangular_points,incenter_rad] = scale_triangle_based_on_barycentric(Points,inward_length)
A = Points(1,:);
B = Points(2,:);
C = Points(3,:);

a = vecnorm(C-B,2,2);
b = vecnorm(A-C,2,2);
c = vecnorm(B-A,2,2);

perimeter = a + b + c;
semi_perimeter = perimeter/2;
area = (semi_perimeter*( semi_perimeter-a )*( semi_perimeter-b )*...
                ( semi_perimeter-c ))^(.5);
incenter_rad = area*2/perimeter;
constant_factor = inward_length/perimeter/incenter_rad;

A_p_bary = [ 1-(b+c)*constant_factor, b*constant_factor, c*constant_factor ];
B_p_bary = [ a*constant_factor, 1-(a+c)*constant_factor, c*constant_factor ];
C_p_bary = [ a*constant_factor, b*constant_factor, 1-(a+b)*constant_factor ];

A_p = A_p_bary(1)*A + A_p_bary(2)*B + A_p_bary(3)*C;
B_p = B_p_bary(1)*A + B_p_bary(2)*B + B_p_bary(3)*C;
C_p = C_p_bary(1)*A + C_p_bary(2)*B + C_p_bary(3)*C;

scales_triangular_points = [A_p;B_p;C_p];


end