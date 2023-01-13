function [Points_new] = cluster_points(Points,cluster_rad)
% size(Points)
num_particles = size(Points,1);
if num_particles~=0
    x_coor = Points(:,1);
    y_coor = Points(:,2);
    z_coor = Points(:,3);
    
    %% alculating the distance matrix
    x_mat = x_coor*ones(1,num_particles);
    y_mat = y_coor*ones(1,num_particles);
    z_mat = z_coor*ones(1,num_particles);
    
    distance_mat = ((x_mat - x_mat').^2 + (y_mat - y_mat').^2 + (z_mat - z_mat').^2).^(.5);
    distance_mat_mask = zeros(size(distance_mat));
    distance_mat_mask(distance_mat<=cluster_rad) = 1;
    identity_mat = diag(ones(num_particles,1));
    distance_mat_mask(identity_mat==1) = 0;
    
    
    %% Updating the points based on clustering
    %% Calculate all the edges based on links
    %% take the center of mass of each edge as a new point and delete the coors of old points
    cluster_links = [];
    for pt = 1:num_particles
        [close_neighbours,~] = find(distance_mat_mask(:,pt) == 1 );
        edges_current = [pt*ones(length(close_neighbours),1),close_neighbours];
        cluster_links = [ cluster_links;edges_current ];
    end
    cluster_links = unique(sort(cluster_links,2),'rows');
    unique_cluster_particles = unique(cluster_links(:));
    unique_charater_notClusterMask = ones(size(unique_cluster_particles));
    % cluster = [];
    
    if isempty(cluster_links)==0
        current_particle = unique_cluster_particles(1);
        loop = 1;
        while (sum(unique_charater_notClusterMask)~=0)
            connected_with_current_particle_mask = ismember(cluster_links(:,1),current_particle)|...
                ismember(cluster_links(:,2),current_particle);
            interesting_edges = cluster_links(connected_with_current_particle_mask,:);
            connected_with_current_particle = unique(interesting_edges(:));
            
            %% deleting also edges containing connected particles only edges
            cluster(loop,1) = {[connected_with_current_particle(:)]};
            connected_with_current_particle_mask = connected_with_current_particle_mask|(ismember(cluster_links(:,1),connected_with_current_particle)&...
                ismember(cluster_links(:,2),connected_with_current_particle));            
            
            cluster_links(connected_with_current_particle_mask,:) = [];
            unique_charater_notClusterMask(loop) =  0;
            loop = loop + 1;
            current_particle = unique_cluster_particles(2);
        end
        clus_centroids = [];
        particles_in_cluster = [];
        for clus = 1:size(cluster,1)
            clus_current = cell2mat(cluster(clus,1));
            if isempty(clus_current)==0
                current_centroid = mean(Points(clus_current,:),1);
                clus_centroids = [clus_centroids;current_centroid];
                particles_in_cluster = unique([reshape(particles_in_cluster(:),1,[]),...
                    reshape(clus_current(:),1,[])]);
            end
        end
        Points(particles_in_cluster,:) = [];
        Points_new = unique([clus_centroids;Points],'rows');
    else
        Points_new = unique(Points,'rows');
    end
else
    Points_new = [];
end
end