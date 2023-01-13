function [barbed_end] = initialize_barbed_ends(pointed_end_base,actin_filament_direction,length)
barbed_end_coorx = pointed_end_base(:,1) + length.*actin_filament_direction(:,1);
barbed_end_coory = pointed_end_base(:,2) + length.*actin_filament_direction(:,2);
barbed_end_coorz = pointed_end_base(:,3) + length.*actin_filament_direction(:,3);
barbed_end = [barbed_end_coorx,barbed_end_coory,barbed_end_coorz];
end