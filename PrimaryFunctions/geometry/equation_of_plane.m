function [output] = equation_of_plane(P1,P2,P3)
    %% equation of the form ax + by + cz + d = 0
    %% output = [a b c d];
    %% input 3 points P1 P2 P3 defining the plane in 3D
    Ax = P1(1);
    Bx = P2(1);
    Cx = P3(1);
    
    Ay = P1(2);
    By = P2(2);
    Cy = P3(2);
    
    Az = P1(3);
    Bz = P2(3);
    Cz = P3(3);
    
    a = (By - Ay)*(Cz - Az) - (Cy - Ay)*(Bz - Az);
    b = (Bz - Az)*(Cx - Ax) - (Cz - Az)*(Bx - Ax);
    c = (Bx - Ax)*(Cy - Ay) - (Cx - Ax)*(By - Ay);
    d = -(a*Ax + b*Ay + c*Az);    
    output_plane(1,1:4) = [a b c d];
    output.plane = output_plane;
    %% Finding vector transformation/components of three points on the plane
    %% unit_basis_vect3 is the vector perpendicular to the plane (so zero component)
    origin_coor = (P1 + P2 + P3)/3 ;
    unit_basis_vect1 = (P1 - origin_coor)/vecnorm((P1 - origin_coor),2,2);
    unit_basis_vect3 = cross(P2 - P1,P3 - P2);
    unit_basis_vect3 = unit_basis_vect3/vecnorm(unit_basis_vect3,2,2);
    unit_basis_vect2 = cross(unit_basis_vect1,unit_basis_vect3);
    
    P1_transformed(1,1) = dot(P1-origin_coor,unit_basis_vect1);
    P1_transformed(1,2) = dot(P1-origin_coor,unit_basis_vect2);
    P1_transformed(1,3) = dot(P1-origin_coor,unit_basis_vect3);
    
    P2_transformed(1,1) = dot(P2-origin_coor,unit_basis_vect1);
    P2_transformed(1,2) = dot(P2-origin_coor,unit_basis_vect2);
    P2_transformed(1,3) = dot(P2-origin_coor,unit_basis_vect3);
    
    P3_transformed(1,1) = dot(P3-origin_coor,unit_basis_vect1);
    P3_transformed(1,2) = dot(P3-origin_coor,unit_basis_vect2);
    P3_transformed(1,3) = dot(P3-origin_coor,unit_basis_vect3);
    
    transformed_vects = [P1_transformed;P2_transformed;P3_transformed];
    output.transformed_coor = transformed_vects;
    
    output.BasisVectors = [unit_basis_vect1;unit_basis_vect2;unit_basis_vect3];
    
    
    
end