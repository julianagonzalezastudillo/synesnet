%% create left hemisphere
clear

sides = {'lh', 'rh'};

% Load brain
% get vertex data (brain)
load('brain_surface.mat');

% Define node generation constants
range_shading_brain = [0.5 1];
range_shading_nodes = [0 1];
distance_factor = 1.1;
sphere_maxRadius = 5;
sphere_minRadius = 1; % 0.4;
binsize = false;

% load nodes
addpath('/Users/juliana.gonzalez/ownCloud/github/synesnet/plots/glb/')
file_name = 'coreness_thr_mean_ctr_selection';

% for each hemisphere
for i = 1:length(sides)
    side = sides{i};
    % get hemisphere data (hemisphereSurface)
    load(append('hemisphere_surface_', side, '.mat'));
    
    % load nodes, node names for each hemispher (Xnet, xyz, color, names, names_idx)
    load(append(file_name, '_', side, '.mat'));
    
    % significant nodes names
    names = string(names);
    names_ = regexprep(names, '\s', ''); % Remove spaces using regular expression
    names_idx_ = names_idx +1; % beacuse it comes from python
    
    % generate spheres
    nodes_h = generate_nodes(Xnet, xyz, sphere_maxRadius, sphere_minRadius, color, distance_factor, binsize);
    
    % plot hemisphere
    hemisphereVerticesIdx = hemisphereSurface.hemisphereVerticesIdx;
    hemisphere.faces = hemisphereSurface.hemisphereFaces;
    
    hemisphere.vertices = brain.vertices(hemisphereVerticesIdx , :);
    hemisphere.vertices = hemisphere.vertices - ones(size(hemisphere.vertices , 1) ,1)* mean(brain.vertices);
    
    tmp_shading_colorHem = brain.shading_pre(hemisphereVerticesIdx) * diff(range_shading_brain);
    tmp_shading_colorHem = tmp_shading_colorHem - min(tmp_shading_colorHem) + range_shading_brain(1);
    shading_colorHem = repmat(tmp_shading_colorHem,1,3);
    hemisphere.color = shading_colorHem;
    
    % load nodes
    nodes_h.colors = nodes_h.colors * diff(range_shading_nodes);
    
    % to center brain and nodes
    load(append(file_name, '.mat'));
    nodes = generate_nodes(Xnet, xyz, sphere_maxRadius, sphere_minRadius, color, distance_factor, binsize);
    
    nodes_h.vertices = nodes_h.vertices - ones(size(nodes_h.vertices , 1) ,1)* mean(nodes.vertices);
    
    nodes_h.xyz = nodes_h.xyz - ones(size(nodes_h.xyz , 1) ,1)* mean(nodes.vertices);
    nodes_h.vertices(:,2) = nodes_h.vertices(:,2) + 6;
    nodes_h.vertices(:,3) = nodes_h.vertices(:,3) - 3;
    
    % Plot the spheres
    allVertices = vertcat(hemisphere.vertices, nodes_h.vertices); 
    allIndices = vertcat(hemisphere.faces, (nodes_h.indices + size(hemisphere.vertices,1)));
    allColors = vertcat(hemisphere.color, nodes_h.colors);
    figure;
    fig = gcf; % Get the current figure handle
    
    trisurf(hemisphere.faces,hemisphere.vertices(:,1),hemisphere.vertices(:,2), ...
        hemisphere.vertices(:,3),'edgecolor','none', 'FaceLighting', 'gouraud', ...
        'AmbientStrength', 0.5, 'FaceVertexCData', hemisphere.color, ...
        'FaceAlpha',.75);
    hold on
    trisurf(nodes_h.indices , ...
        nodes_h.vertices(:,1),nodes_h.vertices(:,2),nodes_h.vertices(:,3), ...
        'edgecolor','none', 'FaceLighting', 'gouraud', 'AmbientStrength', 0.5, ...
        'FaceVertexCData', nodes_h.colors);
    
    % Add label names
    for i = 1:size(names_,1)
        % offset = i * 2 * maxRadius;
        label = sprintf(names_(i));
        % text(nodes_h.xyz(names_idx_(i), 1) -5, ...
        %      nodes_h.xyz(names_idx_(i), 2) +5, ...
        %      nodes_h.xyz(names_idx_(i), 3) +5, label, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 10, 'FontWeight', 'bold');
    end
    
    % Set the desired figure size
    figureWidth = 600;  % Width in pixels
    figureHeight = 600; % Height in pixels
    set(fig, 'Units', 'pixels', 'Position', [100, 100, figureWidth, figureHeight]);
    
    % shading interp
    shading interp
    axis equal;
    grid off
    axis off
    % view(90, 0);  % lateral view
    % view(0, 90);  % top view
    % view(270, 0); % lateral view
end