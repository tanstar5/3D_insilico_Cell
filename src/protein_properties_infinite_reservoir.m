classdef protein_properties_infinite_reservoir
   properties
       types
       concentration_bulk_cytosol       
       concentration_bulk_extracellular
       structure_factor 
       c1_preffered 
       c2_preffered
       mem_binding_energy_per_mol
       
   end
   
   methods       
       function [protein_properties_infinite_reservoir] = ...
               protein_properties_infinite_reservoir(concentration_bulk_cytosol,concentration_bulk_extracellular,structure_factor,c1_preffered,c2_preffered)
           protein_properties_infinite_reservoir.types  = types;
           protein_properties_infinite_reservoir.concentration_bulk_cytosol  = concentration_bulk_cytosol;
           protein_properties_infinite_reservoir.concentration_bulk_extracellular  = concentration_bulk_extracellular;
           protein_properties_infinite_reservoir.structure_factor = structure_factor;
           protein_properties_infinite_reservoir.c1_preffered = c1_preffered;
           protein_properties_infinite_reservoir.c2_preffered = c2_preffered;
       end
       
       function [protein_type_obj] = get_protein_props(protein_obj_all_type,index)
           protein_type_obj.type = protein_obj_all_type.type(index);
           protein_type_obj.concentration_bulk_cytosol = protein_obj_all_type.concentration_bulk_cytosol(index);
           protein_type_obj.concentration_bulk_extracellular = protein_obj_all_type.concentration_bulk_extracellular(index);
           protein_type_obj.structure_factor = protein_obj_all_type.structure_factor(index);
           protein_type_obj.c1_preffered = protein_obj_all_type.c1_preffered(index);
           protein_type_obj.c2_preffered = protein_obj_all_type.c2_preffered(index);
       end
       
   end
end