%%% Protocol for simulating Cell

%% Initialize
figure(10)
cell1 = active_cell_object(1);
actin_length = 10;
cell1 = create_initial_actin_mesh_pointed_ends_from_stl(cell1,1,actin_length); %% function active_cell_bject
hold on;
unit_membrane_mesh_length = .5; 

cell1 = create_initial_membrane_mesh_barbed_ends(cell1,1,unit_membrane_mesh_length);
hold off
