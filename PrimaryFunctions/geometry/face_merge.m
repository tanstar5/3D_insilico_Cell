function [merged_faces,new_face_pairs] = face_merge(faces,face_pairs,common_edge,common_vertex)
old_faces = faces;
common_edge_marker = sum(double(ismember(old_faces,common_edge)),2)==2;
face_idx = 1:size(faces,1);

merged_faces = old_faces;
merged_faces(max(face_idx(common_edge_marker)),:) = NaN(1,3);

face_pairs_to_merge = old_faces(common_edge_marker,:);

face1 = face_pairs_to_merge(1,:);
[~,id]=  find(face1==common_vertex);
face1_circshifted = circshift(face1,1-id);

face2 = face_pairs_to_merge(2,:);
[~,id]=  find(face2==common_vertex);
face2_circshifted = circshift(face2,1-id);

face_circ = [face1_circshifted;face2_circshifted];

% determining the first face to be in order
edge_ordered = [common_vertex, common_edge(common_edge~=common_vertex)];
first_edge_facecirc = face_circ(:,1:2);
face2_determined = face_circ(ismember(first_edge_facecirc,edge_ordered,'rows'),: );
face1_determined = face_circ( not(ismember(first_edge_facecirc,edge_ordered,'rows')),: );

new_face = [face1_determined(1:2),face2_determined(3)];
merged_faces(min(face_idx(common_edge_marker)),:) = new_face;

% modifying the face pair
new_face_pairs = face_pairs;
merged_face_pair_id = sort(face_idx(common_edge_marker),2);
face_pairs_sorted = sort(face_pairs,2);

face_pair_to_del_mask = ismember(face_pairs_sorted,merged_face_pair_id,'rows');
new_face_pairs(face_pair_to_del_mask,:) = NaN(1,2);

face_id_to_del_from_face_pair = max(face_idx(common_edge_marker));
new_face_pairs(new_face_pairs==face_id_to_del_from_face_pair) = min(face_idx(common_edge_marker));

% correcting for NaNs
NaN_idx_in_face_list = face_id_to_del_from_face_pair;
new_face_pairs(new_face_pairs>NaN_idx_in_face_list) = new_face_pairs(new_face_pairs>NaN_idx_in_face_list)-1;
new_face_pairs(isnan(new_face_pairs(:,1)),:) = [];
merged_faces(isnan(merged_faces(:,1)),:) = [];
end