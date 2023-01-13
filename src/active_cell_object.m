classdef active_cell_object
    properties
        ID
        actin
        membrane
        cytosol
    end
    
    methods
        %% Initilising and structure defining functions
        function obj = active_cell_object(arg)
            obj.ID = arg;
        end
        
        function [obj] =  create_initial_actin_mesh_pointed_ends_from_stl(obj,obj_fig_handle,actin_length,actin_mesh_max_length,actin_mesh_min_length)
            %% Reading an STL file object from disk
            filename = uigetfile( '.stl' );
            TR_shape = stlread(filename);
            points = TR_shape.Points;
            elements = TR_shape.ConnectivityList;
            h1 = figure(obj_fig_handle);
            h1_sub1 = subplot(1,1,1);
            %             assignin('base',"TR_actin",TR_shape);
            model = createpde();
            geometryFromMesh(model,points',elements');
            generateMesh(model,'GeometricOrder', 'linear','Hmax',actin_mesh_max_length,'Hmin',actin_mesh_min_length);
            assignin('base',"model",model);
            
            
            %% generating the surface mesh created by the pointed end for normal vector determination
            fprintf('Creating Cytoskeleton\n')
            [T_surf_pointed,pointed_end_meshpoints] = freeBoundary(triangulation(model.Mesh.Elements',model.Mesh.Nodes'));
            fprintf('   Plotting in figure of the surface formed by the pointed ends of the actin mesh %d points \n',...
                size(model.Mesh.Nodes',1));
            
            trisurf(T_surf_pointed,pointed_end_meshpoints(:,1),pointed_end_meshpoints(:,2),pointed_end_meshpoints(:,3),ones(size(pointed_end_meshpoints,1),1),'FaceAlpha',.4); daspect(h1_sub1,[1,1,1]);
            FEMmeshed_pointed_end_surf_triangulation = triangulation(T_surf_pointed,pointed_end_meshpoints);
            fprintf('   Calculating base of the actin filaments from pointed mesh\n');
            actin_obj.pointed_ends =  median_faces(triangulation(T_surf_pointed,pointed_end_meshpoints));
            actin_obj.pointed_end_mesh = FEMmeshed_pointed_end_surf_triangulation;
            hold(h1_sub1,'on'); scatter3(h1_sub1,actin_obj.pointed_ends(:,1),actin_obj.pointed_ends(:,2),actin_obj.pointed_ends(:,3),'filled');hold(h1_sub1,'off');
            
            %% Generating the barbed end surface based on the pointed end surface
            F = faceNormal(FEMmeshed_pointed_end_surf_triangulation,[1:size(T_surf_pointed,1)]');
            hold(h1_sub1,'on');quiver3(actin_obj.pointed_ends(:,1),actin_obj.pointed_ends(:,2),actin_obj.pointed_ends(:,3),F(:,1),F(:,2),F(:,3),1,'color','r');hold(h1_sub1,'off');
            actin_obj.actin_direction_vects = F;
            actin_obj.barbed_end = initialize_barbed_ends(actin_obj.pointed_ends,actin_obj.actin_direction_vects,actin_length);
            hold(h1_sub1,'on'); scatter3(h1_sub1,actin_obj.barbed_end(:,1),actin_obj.barbed_end(:,2),actin_obj.barbed_end(:,3),'filled');hold(h1_sub1,'off');
            %% saving variables to objects
            obj.actin = actin_obj;
            
            %% Exporting variables to workspace
            assignin('base',"FEMmeshed_pointed_end_surf_triangulation",triangulation(T_surf_pointed,pointed_end_meshpoints));
            assignin('base',"pointed_ends",actin_obj.pointed_ends);
        end
        
        function [obj] = create_initial_membrane_mesh_barbed_ends(obj,obj_fig_handle,unit_membrane_mesh_length)
            barbed_ends_coor = obj.actin.barbed_end;
            fprintf('Creating membrane mesh\n');
            %             [TR_connectivity,volume] = convhull(barbed_ends_coor(:,1),barbed_ends_coor(:,2),barbed_ends_coor(:,3));
            %% Calculating the average normal vector for projection into barbed_end_mesh
            pointed_end_mesh = obj.actin.pointed_end_mesh;
            length_actin_filaments = vecnorm(-obj.actin.pointed_ends + obj.actin.barbed_end,2,2);
            [barbed_end_mesh] = initialize_barbed_end_mesh_from_pointed_end(pointed_end_mesh,length_actin_filaments);
            assignin('base',"barbed_end_mesh",barbed_end_mesh);
            fprintf('    Plotting membrane mesh\n');
            h1 = figure(obj_fig_handle);
            h1_sub1 = subplot(1,1,1);
%             trisurf(barbed_end_mesh,'FaceColor','red','FaceAlpha',.3);
            trimesh(barbed_end_mesh,'EdgeColor','b','LineWidth',2,'FaceAlpha',0); daspect([1,1,1])
            
            %% Creating the global membrane mesh from the barbed end surface
            [cytoskeleton_mesh_list] = initialize_cytoskeleton_mesh(barbed_end_mesh);
            d = unit_membrane_mesh_length;
%             [cytoskeleton_mesh_list,new_edgePoints_inserted_points] = segment_all_edge(cytoskeleton_mesh_list,d);            
%             [cytoskeleton_mesh_list,~] = insert_surface_node_points(cytoskeleton_mesh_list,d);            
%             [cytoskeleton_mesh_list,new_surfacePoints_inserted_points] = make_membrane_patch(cytoskeleton_mesh_list,d);
%             [cytoskeleton_mesh_list] = refine_fine_membrane_mesh(cytoskeleton_mesh_list);
%             [global_membrane_mesh,local_membrane_patches] =  ...
%                 concatenate_membrane_patches(cytoskeleton_mesh_list,barbed_end_mesh.Points,new_edgePoints_inserted_points,new_surfacePoints_inserted_points);
            [global_membrane_mesh,cytoskeleton_mesh_list] = create_global_membrane_mesh(cytoskeleton_mesh_list,barbed_end_mesh.Points,d);
%             %% Refining mesh from staright zero area triangles
%             triangles_global_mesh = global_membrane_mesh.ConnectivityList;
%             all_points = global_membrane_mesh.Points;
%             area_zero_noter_logical = zeros(size(triangles_global_mesh,1),1);
%             for tri = 1:size(triangles_global_mesh,1)
%                 tri_ind = triangles_global_mesh(tri,:);
%                 Points_current = all_points(tri_ind,:);
%                 if abs(real(area_triangle(Points_current)))<=10e-10
%                     area_zero_noter_logical(tri) = 1;
%                     fprintf('Zero area triangle detected with area %d in triangle %d\n',area_triangle(Points_current), tri );
%                 else
%                     area_zero_noter_logical(tri) = 0;
%                     fprintf('Zero area triangle not detected with area %d in triangle %d\n',area_triangle(Points_current), tri );
%                 end
%             end
%             triangles_global_mesh(area_zero_noter_logical==1,:)=[];
%             global_membrane_mesh = triangulation(triangles_global_mesh,all_points);
            
            %%
            obj.membrane.GlobalMesh = global_membrane_mesh;
%             obj.membrane.LocalMesh = local_membrane_patches;
            trisurf(global_membrane_mesh,'FaceColor','r','FaceAlpha',.5); daspect([1,1,1])
            hold off
            
        end
        
        %% Dynamic functions
        
        
        
    end
    
end