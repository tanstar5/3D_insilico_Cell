function [projected_face_areas] = calculate_projected_area_of_faces(Normal_vect,face_normals,face_areas)
num_faces = size(face_normals,1);
cosine_angle_to_each_face = sum((Normal_vect'*ones(1,num_faces))'.*face_normals,2);
projected_face_areas = face_areas.*cosine_angle_to_each_face;
end