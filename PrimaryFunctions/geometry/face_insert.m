function [new_faces,new_face_pairs] = face_insert(faces,face_pairs,edge_to_del,common_vertex,opposite_vertex)
new_faces = faces;
face_idx = 1:size(faces,1);



to_del_face_mask = sum(double(ismember(faces,edge_to_del)),2)==2;

face_id_to_del = face_idx(to_del_face_mask);
if(isempty(face_id_to_del)==1)
   disp('eroor') 
end

new_faces(face_id_to_del,:) = NaN(1,3);

face_deleted = faces(face_id_to_del,:);

[~,id]=find(face_deleted==common_vertex);
face_deleted_circshift = circshift(face_deleted,1-id);

face1_formed = [face_deleted_circshift(1:2),opposite_vertex];
face2_formed = [face_deleted_circshift(1),opposite_vertex,face_deleted_circshift(3)]; 

new_faces(face_id_to_del,:) = face1_formed;
new_faces(end+1,:) = face2_formed;

new_face_formed = [face1_formed;face2_formed];
corresponding_indices = [face_id_to_del;size(new_faces,1)];

% adjacent_face_ids = (ismember(face_pairs,face_id_to_del));



% correcting facePairs
face_id_deleted_containingPairMask = sum(double(ismember(face_pairs,face_id_to_del)),2)==1;
facepairs_to_correct = face_pairs(face_id_deleted_containingPairMask,:);

adjacent_face_ids = unique(facepairs_to_correct(facepairs_to_correct~=face_id_to_del));
adjacent_faces = new_faces(adjacent_face_ids,:);

new_face_formed_id = [face_id_to_del,size(new_faces,1)];
corresponding_face = NaN(2,1);
for ad_fc_id = 1:2
    current_adjacent_face = adjacent_faces(ad_fc_id,:);
    if length(intersect( current_adjacent_face, new_face_formed(1,:) ))==2
        corresponding_face(ad_fc_id) = new_face_formed_id(1);
    elseif length(intersect( current_adjacent_face, new_face_formed(2,:) ))==2
        corresponding_face(ad_fc_id) = new_face_formed_id(2);
    else
        fprintf('    >>>ERROR IN face_insert\n');
    end
end

new_face_pairs_formed = NaN(2,2);
new_face_pairs_formed(:,1) = adjacent_face_ids;
new_face_pairs_formed(:,2) = corresponding_face;
new_face_pairs_formed(3,:) = [corresponding_face(1),corresponding_face(2)];

face_pairs(face_id_deleted_containingPairMask,:) = [];
face_pairs = [face_pairs;new_face_pairs_formed];



% facepairs_to_corrected =facepairs_to_correct;
% for f_id = 1:2
%    fc_pair_current = facepairs_to_correct(f_id,:);
%    adjacent_face_id = fc_pair_current(fc_pair_current~=face_id_to_del);
%    if length(intersect( new_faces(adjacent_face_id,:),new_face_formed(f_id,:)))==2
%        fc_pair_current(fc_pair_current~=adjacent_face_id)= corresponding_indices(f_id);
%    else
%        fc_pair_current(fc_pair_current~=adjacent_face_id)= corresponding_indices( mod(f_id,2)+1 );
%    end
%    facepairs_to_corrected(f_id,:) = fc_pair_current;
% end
% face_pairs(face_id_deleted_containingPairMask,:) =  facepairs_to_corrected;
% face_pairs(end+1,:) = [face_id_to_del;size(new_faces,1)];
new_face_pairs = face_pairs;

end