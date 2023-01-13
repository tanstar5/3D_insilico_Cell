function [membrane_mesh] = import_surface_mesh_from_stl(filename,hmin,hmax)
model = createpde;
importGeometry(model,filename);
mesh_default = generateMesh(model,'Hmax',hmax,'Hmin',hmin,'GeometricOrder','linear');
points = (mesh_default.Nodes)';
connectivity = (mesh_default.Elements)';

tetrahedron_mesh = triangulation(connectivity,points);
[F,P] = freeBoundary(tetrahedron_mesh);

membrane_mesh = triangulation(F,P);


end