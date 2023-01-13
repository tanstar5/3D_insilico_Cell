classdef parameters
    
    properties
        % constants
        kB
        lipid_head_area
        
        % physical parameters
        kappa
        surface_modulous
        kappa_coefficient
        vanderwaals
        penalty_factor
        temperature
        spont_curv_lipids
        radius
        
        % simulation parameters
        time_steps
        step_size
        particle_packet_size
        patch_min_edge 
        patch_max_edge
        
    end
    
    methods
        
        function obj = ...
            parameters(kB,lipid_head_area,kappa,surface_modulous,kappa_coefficient,vanderwaals,penalty_factor,temperature,spont_curv_lipids,radius,time_steps,step_size,particle_packet_size,patch_min_edge,patch_max_edge)
            if nargin ~= 0
                obj.kB = kB;
                obj.lipid_head_area = lipid_head_area;
                obj.kappa = kappa;
                obj.surface_modulous = surface_modulous;
                obj.kappa_coefficient = kappa_coefficient;
                obj.vanderwaals = vanderwaals;
                obj.penalty_factor= penalty_factor;
                obj.temperature= temperature;
                obj.spont_curv_lipids= spont_curv_lipids;
                obj.radius= radius;
                obj.time_steps= time_steps;
                obj.step_size= step_size;
                obj.particle_packet_size= particle_packet_size;
                obj.patch_min_edge= patch_min_edge;
                obj.patch_max_edge= patch_max_edge;
                
            end
        end
        
    end
end