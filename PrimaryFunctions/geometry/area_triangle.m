function [area] = area_triangle(Points)
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
end