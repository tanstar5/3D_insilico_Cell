obj = obj_list(3,1);
obj.neighbours_lipid_compositions_up = lipid_comp(obj.neighbours,:);
pos_objs(i,:) = obj.Pos;
[grad_vect,grad_vect_along_surf] = calculate_gradient_green_gauss_method(obj)