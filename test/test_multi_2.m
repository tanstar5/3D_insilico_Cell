input_file_name = 'circle_rad_20_and_5.stl';
destinationFolder = 'T:\1Data\Tanumoy\selforganization\active_cell_version42_two_component_gradient_term\outputs';

kB = 1;
lipid_head_area = 0.5; % nm^2
kappa = 27*kB*300;
surface_modulous = kappa/4*12; %% relation between kappa and surface modulous is kappa = surface_molous*t^2 where t is bilayer thickness (wikipedia: search Lipid bilayer mechanics)
% or in popular literature kappa =
% surface_modulous*t^2/12
kappa_coefficient = 1;
vanderwaals = 2*kB*300;

temperature = 300;
penalty_factor = 8*kB*300;

radius = 300;
spont_curv_lipids = [-0.2+1/radius,1/radius+0.2];

time_steps = 1000;
step_size = 0.3;
particle_packet_size = 30;
patch_min_edge = 1.4;
patch_max_edge = 1.4;



obj = ...
    parameters(kB,lipid_head_area,kappa,surface_modulous,kappa_coefficient,vanderwaals,penalty_factor, ...
    temperature,spont_curv_lipids,radius,time_steps,step_size,particle_packet_size,patch_min_edge,patch_max_edge);

%% Create parameters for all experimental conditions
parameters_obj1 = obj;
% parameters_obj1.temperature = 274;
% parameters_obj1.penalty_factor = 2*parameters_obj1.penalty_factor;
parameters_obj1.spont_curv_lipids = [-0.05+1/radius,1/radius+0.05];

parameters_obj2 = obj;
% parameters_obj2.temperature = 280;
% parameters_obj2.penalty_factor = 4*parameters_obj2.penalty_factor;
parameters_obj2.spont_curv_lipids = [-0.08+1/radius,1/radius+0.08];

parameters_obj3 = obj;
% parameters_obj3.temperature = 300;
% parameters_obj3.penalty_factor = 8*parameters_obj3.penalty_factor;
parameters_obj3.spont_curv_lipids = [-0.1+1/radius,1/radius+0.1];

parameters_obj4 = obj;
% parameters_obj4.temperature = 350;
% parameters_obj4.penalty_factor = 16*parameters_obj4.penalty_factor;
parameters_obj4.spont_curv_lipids = [-0.15+1/radius,1/radius+0.15];

parameters_obj5 = obj;
% parameters_obj5.temperature = 400;
% parameters_obj5.penalty_factor = 32*parameters_obj5.penalty_factor;
parameters_obj5.spont_curv_lipids = [-0.18+1/radius,1/radius+0.18];

parameters_obj6 = obj;
% parameters_obj6.temperature = 500;
% parameters_obj6.penalty_factor = 64*parameters_obj6.penalty_factor;
parameters_obj6.spont_curv_lipids = [-0.2+1/radius,1/radius+0.2];


delete(gcp('nocreate')); %delete the current pool

parpool(6)
parfor K = 1:6
    if K == 1
        test_simulation_parameters(parameters_obj1,input_file_name,destinationFolder,1);
    end
    if K == 2
        test_simulation_parameters(parameters_obj2,input_file_name,destinationFolder,2);
    end
    
    if K == 3
        test_simulation_parameters(parameters_obj3,input_file_name,destinationFolder,3);
    end
    
    if K == 4
        test_simulation_parameters(parameters_obj4,input_file_name,destinationFolder,4);
    end
    
    if K == 5
        test_simulation_parameters(parameters_obj5,input_file_name,destinationFolder,5);
    end
    
    if K == 6
        test_simulation_parameters(parameters_obj6,input_file_name,destinationFolder,6);
    end
    
end


%% Run experiment
delete(gcp('nocreate')); %delete the current pool

parpool(6)
parfor K = 1:6
    if K == 1
        
        
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulation(parameters_obj1,input_file_name,destinationFolder,30,1);
    end
    if K == 2
        
        
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulation(parameters_obj2,input_file_name,destinationFolder,30,2);
    end
    
    if K == 3
        
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulation(parameters_obj3,input_file_name,destinationFolder,30,3);
    end
    
    if K == 4
        
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulation(parameters_obj4,input_file_name,destinationFolder,30,4);
    end
    
    if K == 5
        
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulation(parameters_obj5,input_file_name,destinationFolder,30,5);
    end
    
    if K == 6
        
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulation(parameters_obj6,input_file_name,destinationFolder,30,6);
    end
    
end
