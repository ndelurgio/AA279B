function visualize_planet(planet,X,Y,Z)
% Creates plot of planet sourced from image file
% Input:
%   planet - string name of planet

    if planet == "Mars"
        R = 3390; % Mars radius [km]
        img_file = "8k_mars.jpg";
    else
        R = 6371; % Earth radius [km]
        img_file = "Earth2.jpg";
    end
    I = imread(img_file);
    [x,y,z] = sphere;              % create a sphere 
    % Move sphere to specified location
    if nargin > 1
        x = x+X; y = y+Y; z = z+Z;
    end
    s = surface(R*x,R*y,R*z);            % plot spherical surface
    s.FaceColor = 'texturemap';    % use texture mapping
    s.CData = flipud(I);           % set color data to topographic data
    s.EdgeColor = 'none';          % remove edges
    s.FaceLighting = 'gouraud';    % preferred lighting for curved surfaces
    s.SpecularStrength = 0.4;      % change the strength of the reflected light
    
    light('Position',[1 0 1])     % add a light
    axis equal; grid on;
end