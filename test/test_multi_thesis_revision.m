input_file_name = 'circle_rad_20_wedges.stl';
destinationFolder = 'T:\1Data\Tanumoy\selforganization\active_cell_version42_two_component_gradient_term\outputThesis';

%% Physical parameters
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

time_steps = 300;
step_size = 0.3;
particle_packet_size = 30;
patch_min_edge = 1.4;
patch_max_edge = 1.4;


obj = ...
    parameters(kB,lipid_head_area,kappa,surface_modulous,kappa_coefficient,vanderwaals,penalty_factor, ...
    temperature,spont_curv_lipids,radius,time_steps,step_size,particle_packet_size,patch_min_edge,patch_max_edge);

%% Create parameters for all experimental conditions
parameters_obj = obj;
parameters_obj1 = obj;
parameters_obj2 = obj;
parameters_obj3 = obj;
parameters_obj4 = obj;
parameters_obj5 = obj;
parameters_obj6 = obj;
% [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj,input_file_name,destinationFolder,5,1,1,1);
% % [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj1,input_file_name,destinationFolder,3,1,1);
% [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj,input_file_name,destinationFolder,5,1,1,0);
delete(gcp('nocreate')); %delete the current pool

%% Curvature based accumulation with and without membrane fluctuation
% parpool(2)
% parfor K = 1:6
%     if K == 1      
%         [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj,input_file_name,destinationFolder,5,1,1,1);
%     end
%     if K == 2
%         [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj,input_file_name,destinationFolder,5,2,1,0);
%     end
%     
%     
% end

%% Curvature based accumulation  without membrane fluctuation at different temperature
% parpool(6)
% parameters_obj1.temperature = 278;
% parameters_obj2.temperature = 290; 
% parameters_obj3.temperature = 300; 
% parameters_obj4.temperature = 310;
% parameters_obj5.temperature = 320;
% parameters_obj6.temperature = 350;
% parfor K = 1:6
%     if K == 1
%          
%         [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj1,input_file_name,destinationFolder,5,1,1,0);
%     end
%     if K == 2
%         
%         [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj2,input_file_name,destinationFolder,5,2,1,0);
%     end
%     if K == 3
%         
%         [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj3,input_file_name,destinationFolder,5,3,1,0);
%     end
%     if K == 4
%         
%         [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj4,input_file_name,destinationFolder,5,4,1,0);
%     end
%     if K == 5
%          
%         [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj5,input_file_name,destinationFolder,5,5,1,0);
%     end
%     if K == 6
%          
%         [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj6,input_file_name,destinationFolder,5,6,1,0);
%     end
%     
%     
% end

%% Curvature based accumulation  without membrane fluctuation at different temperature with no vanderwaal
destinationFolder = 'T:\1Data\Tanumoy\selforganization\active_cell_version42_two_component_gradient_term\outputThesis\Diff_temperature_withoutVanderwaals';
parpool(6)
parameters_obj1.temperature = 278;
parameters_obj2.temperature = 290; 
parameters_obj3.temperature = 300; 
parameters_obj4.temperature = 310;
parameters_obj5.temperature = 320;
parameters_obj6.temperature = 350;
parameters_obj1.vanderwaals  = 0;
parameters_obj2.vanderwaals  = 0; 
parameters_obj3.vanderwaals  = 0; 
parameters_obj4.vanderwaals  = 0;
parameters_obj5.vanderwaals  = 0;
parameters_obj6.vanderwaals  = 0;
parfor K = 1:6
    if K == 1
         
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj1,input_file_name,destinationFolder,5,1,1,0);
    end
    if K == 2
        
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj2,input_file_name,destinationFolder,5,2,1,0);
    end
    if K == 3
        
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj3,input_file_name,destinationFolder,5,3,1,0);
    end
    if K == 4
        
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj4,input_file_name,destinationFolder,5,4,1,0);
    end
    if K == 5
         
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj5,input_file_name,destinationFolder,5,5,1,0);
    end
    if K == 6
         
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj6,input_file_name,destinationFolder,5,6,1,0);
    end
    
    
end
parameters_obj1 = obj;
parameters_obj2 = obj;
parameters_obj3 = obj;
parameters_obj4 = obj;
parameters_obj5 = obj;
parameters_obj6 = obj;
delete(gcp('nocreate'));


%% Curvature based accumulation  without membrane fluctuation at same temperature but different kappa
destinationFolder = 'T:\1Data\Tanumoy\selforganization\active_cell_version42_two_component_gradient_term\outputThesis\Room_temperature_diffKappa';
parpool(6)

parameters_obj1.kappa  = 0.2*27*kB*300;
parameters_obj2.kappa  = 0.5*27*kB*300; 
parameters_obj3.kappa  = 0.8*27*kB*300; 
parameters_obj4.kappa  = 27*kB*300;
parameters_obj5.kappa  = 1.2*27*kB*300;
parameters_obj6.kappa  = 1.5*27*kB*300;
parfor K = 1:6
    if K == 1
         
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj1,input_file_name,destinationFolder,5,1,1,0);
    end
    if K == 2
        
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj2,input_file_name,destinationFolder,5,2,1,0);
    end
    if K == 3
        
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj3,input_file_name,destinationFolder,5,3,1,0);
    end
    if K == 4
        
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj4,input_file_name,destinationFolder,5,4,1,0);
    end
    if K == 5
         
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj5,input_file_name,destinationFolder,5,5,1,0);
    end
    if K == 6
         
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj6,input_file_name,destinationFolder,5,6,1,0);
    end
    
    
end
parameters_obj1 = obj;
parameters_obj2 = obj;
parameters_obj3 = obj;
parameters_obj4 = obj;
parameters_obj5 = obj;
parameters_obj6 = obj;
delete(gcp('nocreate'));


%% Curvature based accumulation  without membrane fluctuation at same temperature but different vanderwaals
destinationFolder = 'T:\1Data\Tanumoy\selforganization\active_cell_version42_two_component_gradient_term\outputThesis\Room_temperature_diffvanderwaala';
parpool(6)

parameters_obj1.vanderwaals  = 0.2*2*kB*300;
parameters_obj2.vanderwaals  = 0.5*2*kB*300; 
parameters_obj3.vanderwaals  = 0.8*2*kB*300; 
parameters_obj4.vanderwaals  = 2*kB*300;
parameters_obj5.vanderwaals  = 1.2*2*kB*300;
parameters_obj6.vanderwaals  = 1.5*2*kB*300;
parfor K = 1:6
    if K == 1
         
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj1,input_file_name,destinationFolder,5,1,1,0);
    end
    if K == 2
        
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj2,input_file_name,destinationFolder,5,2,1,0);
    end
    if K == 3
        
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj3,input_file_name,destinationFolder,5,3,1,0);
    end
    if K == 4
        
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj4,input_file_name,destinationFolder,5,4,1,0);
    end
    if K == 5
         
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj5,input_file_name,destinationFolder,5,5,1,0);
    end
    if K == 6
         
        [patch_min_edge_scaled,patch_max_edge_scaled] = run_simulationv2(parameters_obj6,input_file_name,destinationFolder,5,6,1,0);
    end
    
    
end
parameters_obj1 = obj;
parameters_obj2 = obj;
parameters_obj3 = obj;
parameters_obj4 = obj;
parameters_obj5 = obj;
parameters_obj6 = obj;
delete(gcp('nocreate'));





