function [median_coor] = median_faces(triangulation)

ConnectivityList = triangulation.ConnectivityList;
x_cor_points = triangulation.Points(:,1);
y_cor_points = triangulation.Points(:,2);
z_cor_points = triangulation.Points(:,3);

x_mat = x_cor_points(ConnectivityList);
y_mat = y_cor_points(ConnectivityList);
z_mat = z_cor_points(ConnectivityList);

x_mean = mean(x_mat,2);
y_mean = mean(y_mat,2);
z_mean = mean(z_mat,2);

median_coor = [x_mean y_mean z_mean];

    
end