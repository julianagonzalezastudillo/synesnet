clear all

% Parameters 
aa = 0.1;
range = [0 1];
distance_factor = 1.15;

% load brain
load('simple_brain_surface.mat');

% load nodes
addpath('/Users/juliana.gonzalez/ownCloud/github/synesnet/plots/glb/')
% nodes_file = 'plots/glb/strength_thr_t-val.mat';
% nodes_file = 'plots/glb/coreness_norm_thr_t-val.mat';
% nodes_file = 'plots/glb/coreness_norm_by_rand_conserve_strenght_distribution_thr_t-val_selection';
% nodes_file = 'plots/glb/coreness_norm_thr_mean_syn.mat';
% nodes_file = 'plots/glb/coreness_norm_thr_mean_ctr.mat';
% nodes_file = 'plots/glb/coreness_norm_strength_selection_thr_mean_syn.mat';
% nodes_file = 'plots/glb/coreness_norm_strength_selection_thr_t-val.mat';
% nodes_file = 'plots/glb/coreness_norm_strength_selection_thr_mean_syn.mat';
% nodes_file = 'plots/glb/literature.mat';
nodes_file = 'coreness_thr_mean_ctr_selection.mat';
load(nodes_file);

% load significant nodes names

names = string(names);
names = regexprep(names, '\s', ''); % Remove spaces using regular expression
names_idx = names_idx +1; % beacuse it comes from python

% generate spheres
sphere_maxRadius = 4;
sphere_minRadius = 0.4;
nodes = generate_nodes(Xnet, xyz, sphere_maxRadius, sphere_minRadius, color, distance_factor);
nodes.vertices = nodes.vertices - ones(size(nodes.vertices , 1) ,1)* mean(nodes.vertices);
nodes.vertices(:,2) = nodes.vertices(:,2) + 5;
nodes.vertices(:,1) = nodes.vertices(:,1) - 2;

% center
brain.vertices = brain.vertices - ones(size(brain.vertices , 1) ,1)* mean(brain.vertices);

tmp_shading_color = brain.shading_pre * diff(range);
tmp_shading_color = tmp_shading_color - min(tmp_shading_color) + range(1);
shading_color = repmat(tmp_shading_color,1,3);
brain.color = shading_color;

n_color = [0.9961 0.7461 0.2461];
numVertices = size(brain.vertices, 1);
color_v = repmat(n_color, numVertices, 1);  
color_v = repmat(tmp_shading_color, 1, 3) .* color_v;  % Apply shading to the color vector
brain.color = color_v;


nodes.colors = nodes.colors * diff(range);
nodes.xyz = nodes.xyz - ones(size(nodes.xyz , 1) ,1)* mean(nodes.vertices);

%% Plot the spheres
% allVertices = vertcat(hemisphere.vertices, nodes_h.vertices); 
% allIndices = vertcat(hemisphere.faces, (nodes_h.indices + size(hemisphere.vertices,1)));
% allColors = vertcat(hemisphere.color, nodes_h.colors);
figure;
fig = gcf; % Get the current figure handle
% trisurf(allIndices,allVertices(:,1),allVertices(:,2),allVertices(:,3), ...
%     'edgecolor','none', 'FaceLighting', 'gouraud', 'AmbientStrength', 0.5, ...
%     'FaceVertexCData', allColors);
trisurf(brain.faces,brain.vertices(:,1),brain.vertices(:,2), ...
    brain.vertices(:,3),'edgecolor','none', 'FaceLighting', 'gouraud', ...
    'AmbientStrength', 0.5, 'FaceVertexCData', brain.color, ...
    'FaceAlpha',.75);
hold on
trisurf(nodes.indices , ...
    nodes.vertices(:,1),nodes.vertices(:,2),nodes.vertices(:,3), ...
    'edgecolor','none', 'FaceLighting', 'gouraud', 'AmbientStrength', 0.5, ...
    'FaceVertexCData', nodes.colors);

% Add label names
for i = 1:size(names,1)
    % offset = i * 2 * maxRadius;
    label = sprintf(names(i));
    % text(nodes.xyz(names_idx(i), 1) -5, ...
    %      nodes.xyz(names_idx(i), 2) +28, ...
    %      nodes.xyz(names_idx(i), 3) - 9, label, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 10, 'FontWeight', 'bold');
end

% Set the desired figure size
figureWidth = 600;  % Width in pixels
figureHeight = 600; % Height in pixels
set(fig, 'Units', 'pixels', 'Position', [100, 100, figureWidth, figureHeight]);

% trisurf(allIndices, allVertices(:, 1), allVertices(:, 2), allVertices(:, 3), C, 'edgecolor','none');
% light
% shading interp
shading interp
axis equal;

xlabel('X');
ylabel('Y');
zlabel('Z');
grid off
axis off
% view(0, 90);

