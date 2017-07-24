function h = plotMesh(mesh, f, az, el, flag)
    
    mesh1 = mesh; 
    mesh1.F = mesh.surface.TRIV'; 
    mesh1.V = [mesh.surface.X'; mesh.surface.Y'; mesh.surface.Z'];
    mesh1.Nf = mesh.Nf'; 
    mesh1.Nv = mesh.Nv'; 
    
    %Default: Rotate the shape by pi/2
    if nargin < 5
        flag = 'Rotate'; 
    end
    
    if strcmp(flag, 'Rotate')
        mesh1.V = [mesh1.surface.Y'; mesh1.surface.Z'; mesh1.surface.X'];     
        per = [2, 3, 1];
    else
        per = [1, 2, 3];
    end
    
    if var(f) > 1E-10
        h = trimesh(mesh1.F', mesh1.V(1,:)', mesh1.V(2,:)' ,mesh1.V(3,:)', f, 'FaceColor', 'interp', 'EdgeColor', 'interp', ...
            'AmbientStrength', 0.35, 'DiffuseStrength', 0.65, 'SpecularStrength', 0.0, 'FaceLighting', 'gouraud', ...
            'SpecularExponent', 10, 'VertexNormals', -mesh1.Nv(per,:)', 'BackFaceLighting', 'reverselit', 'LineStyle', 'none');
    else
        h = trimesh(mesh1.F', mesh1.V(1,:)', mesh1.V(2,:)' ,mesh1.V(3,:)', 'FaceColor', 'w', 'EdgeColor', 'interp', ...
            'AmbientStrength', 0.35, 'DiffuseStrength', 0.65, 'SpecularStrength', 0.0, 'FaceLighting', 'gouraud', ...
            'SpecularExponent', 10, 'VertexNormals', -mesh1.Nv(per,:)', 'BackFaceLighting', 'reverselit', 'LineStyle', 'none');
    end
    
    colormap('hot'); colormap(flipud(colormap));  
    set(gcf, 'Color', 'w', 'Renderer', 'OpenGL');
    set(gca, 'Projection', 'perspective');    
    axis equal;
    axis off;
    view(az,el);
    camlight(30,40);
end