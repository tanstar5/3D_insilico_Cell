function [distance_from_plane] = measure_distance_from_plane(point_in_plane,Normal_vect,points)
num_points = size(points,1);
% distance_from_plane = NaN(size(points,1),1);
distance_vect_from_point_in_plane = (point_in_plane'*ones(1,num_points))'-points;
distance_from_plane = distance_vect_from_point_in_plane*(Normal_vect');

end