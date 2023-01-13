function [masked_mat] = mask_mat_to_sequence(mat_to_mask,id_sequence)
masked_mat = NaN(size(mat_to_mask));
for i = 1:length(id_sequence)
    masked_mat( ismember(mat_to_mask,id_sequence(i)) ) = i;
end
end