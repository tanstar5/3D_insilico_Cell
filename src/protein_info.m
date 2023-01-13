classdef protein_info
    
    properties
        type
        size_radius
        shape_function
    end
    
    methods        
        function obj = protein_info(type,size_rad,shape_func)
            if nargin ~= 0
                obj.type = type;
                obj.size_radius = size_rad;
                obj.shape_function = shape_func;
            end
        end    
    end    
end