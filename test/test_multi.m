input_file_name = 'circle_rad_20.stl';
destinationFolder = 'C:\Users\tanum\OneDrive\Documents\PhD work\workfolder_surface\Self_Organization_project\active_cell_version41_two_component_gradient_term\outputs\test_multi';

kB = 1;
lipid_head_area = 0.5; % nm^2
kappa = 60*kB*300;
surface_modulous = kappa/4*12; %% relation between kappa and surface modulous is kappa = surface_molous*t^2 where t is bilayer thickness (wikipedia: search Lipid bilayer mechanics)
% or in popular literature kappa =
% surface_modulous*t^2/12
temperature = 100;
penalty_factor = 0.05*kB*temperature;

radius = 300;
spont_curv_lipids = [-0.2+1/radius,1/radius+0.2];

time_steps = 1000;
step_size = 0.3;
particle_packet_size = 30;
patch_min_edge = 1.4;
patch_max_edge = 1.4;



obj = ...
    parameters(kB,lipid_head_area,kappa,surface_modulous,penalty_factor, ...
    temperature,spont_curv_lipids,radius,time_steps,step_size,particle_packet_size,patch_min_edge,patch_max_edge);


delete(gcp('nocreate')); %delete the current pool

parpool(9)
parfor K = 1:2
    if K == 1
        parameters_obj1 = obj;
        parameters_obj1.temperature = 100;
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulation(parameters_obj1,input_file_name,destinationFolder,1,1);
    end
    if K == 2
        parameters_obj2 = obj;
        parameters_obj2.temperature = 150;
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulation(parameters_obj2,input_file_name,destinationFolder,1,2);
    end
    
    if K == 3
        parameters_obj3 = obj;
        parameters_obj3.temperature = 200;
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulation(parameters_obj3,input_file_name,destinationFolder,1,3);
    end
    
    if K == 4
        parameters_obj4 = obj;
        parameters_obj4.temperature = 250;
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulation(parameters_obj4,input_file_name,destinationFolder,1,4);
    end
    
    if K == 5
        parameters_obj5 = obj;
        parameters_obj5.temperature = 280;
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulation(parameters_obj5,input_file_name,destinationFolder,1,5);
    end
    
    if K == 6
        parameters_obj6 = obj;
        parameters_obj6.temperature = 300;
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulation(parameters_obj6,input_file_name,destinationFolder,1,6);
    end
    
end
