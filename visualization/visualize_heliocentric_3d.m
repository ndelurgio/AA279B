function visualize_heliocentric_3d()
% Visualizes a body of specified radius
% Outputs: Plot of body
    scatter3(0,0,0,72,[227,193,2]/255,'filled','o')
    hold on; grid on; %axis equal;
    title('Heliocentric Inertial Frame')
    xlabel('X [AU]');
    ylabel('Y [AU]');
    zlabel('Z [AU]')
    % view(3);
end