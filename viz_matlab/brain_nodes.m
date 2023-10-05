clear

% Parameters 
aa = 0.1;
range_shading_brain = [0.6 1];
range_shading_nodes = [0 1];
distance_factor = 1.05;
sizeType = 'unique';

% load brain
load('brain_surface.mat');

% load nodes
addpath('/Users/juliana.gonzalez/ownCloud/github/synesnet/plots/glb/')
nodes_file = 'coreness_thr_mean_ctr_selection.mat';

load(nodes_file);

% load significant nodes names

names = string(names);
names = regexprep(names, '\s', ''); % Remove spaces using regular expression
names_idx = names_idx +1; % beacuse it comes from python

% generate spheres
sphere_maxRadius = 5.1;
sphere_minRadius = 5;
nodes = generate_nodes(Xnet, xyz, sphere_maxRadius, sphere_minRadius, color, distance_factor, sizeType);
nodes.vertices = nodes.vertices - ones(size(nodes.vertices , 1) ,1) * mean(nodes.vertices);
nodes.vertices(:,2) = nodes.vertices(:,2) + 6;
nodes.vertices(:,3) = nodes.vertices(:,3) - 3;

% center brain
brain.vertices = brain.vertices - ones(size(brain.vertices , 1) ,1) * mean(brain.vertices);

% brain shading
tmp_shading_color = brain.shading_pre * diff(range_shading_brain);
tmp_shading_color = tmp_shading_color - min(tmp_shading_color) + range_shading_brain(1);
shading_color = repmat(tmp_shading_color,1,3);
brain.color = shading_color;

% brain corlor
n_color = [1 1 1]; % grey
numVertices = size(brain.vertices, 1);
color_v = repmat(n_color, numVertices, 1);  
color_v = repmat(tmp_shading_color, 1, 3) .* color_v;  % Apply shading to the color vector
brain.color = color_v;

% nodes shading
nodes.colors = nodes.colors * diff(range_shading_nodes);

% nodes position
nodes.xyz = nodes.xyz - ones(size(nodes.xyz , 1) ,1) * mean(nodes.vertices);

% Plot the spheres
figure;
fig = gcf; % Get the current figure handle
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

shading interp
axis equal;

xlabel('X');
ylabel('Y');
zlabel('Z');
grid off
axis off
% view(0, 90);

