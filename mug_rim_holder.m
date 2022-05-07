clear all
close all
clc

num_reps = 12;
gap_size = 3.6;
resolution = 0.25;
clearance = 0.8;
mug_inner_radius = 51/2;
mug_rim_radius = 75.5/2;
mug_rim_minor_radius_h = 1.85;
mug_rim_minor_radius_v = 1.5;
cutout_radius = 7.25;
outer_descent = 10.5;
outer_radius_increase = 2.5;

polygon = csvread('border_subregion.csv');
polygon(:,2) = 500-polygon(:,2); % compensate for the original being in image coords

polygon = flipud(polygon); % reverse the order (makes it make more sense visually)

x_shift = mug_rim_radius - polygon(287,1) - 0.25;
polygon(:,1) = polygon(:,1) + x_shift;

y_shift = mug_rim_minor_radius_v - polygon(287,2);
polygon(:,2) = polygon(:,2) + y_shift;

% Rotate the polygon ever so slightly to move it away from the mug wall
x_offset_from_rot = -0.6;
rotation_center = [mug_inner_radius+clearance polygon(17,2)];
theta = 0.035;

rot_matrix = [cos(theta) -sin(theta); sin(theta) cos(theta)];
polygon = polygon - repmat(rotation_center, size(polygon,1), 1);
polygon = polygon * rot_matrix';
polygon = polygon + repmat(rotation_center, size(polygon,1), 1);

n_divs = ceil((pi*mug_rim_minor_radius_v)/resolution);
n_divs2 = ceil((pi/2*cutout_radius)/resolution); % larger circle, but we need number of divs now
theta = linspace(0.9*pi,pi/2,n_divs + 1)';
theta2 = linspace(pi/2,0,n_divs2 + 1)';
theta = [theta; theta2(2:end)];

polygon = polygon(17:240,:); % Keep only the part of the polygon needed to join with the rim
polygon(1,1) = mug_inner_radius + clearance;

polygon = [polygon;...
mug_rim_radius + mug_rim_minor_radius_h * cos(theta) mug_rim_minor_radius_v * sin(theta)];

n_divs = ceil(sqrt(outer_descent^2+outer_radius_increase^2)/resolution);
y_values = linspace(0,-outer_descent,n_divs + 1)';
x_values = linspace(mug_rim_radius + mug_rim_minor_radius_h,...
                    mug_rim_radius + mug_rim_minor_radius_h + outer_radius_increase,n_divs + 1)';
x_values = x_values(2:end);
y_values = y_values(2:end);

polygon = [polygon;...
x_values y_values];

% Divide the polygon into pieces for this experiment
crit_pt = find(polygon(:,2)==max(polygon(:,2)),1,'first');
polygonA = polygon(1:crit_pt,:);
polygonB = polygon(crit_pt:crit_pt+numel(theta2)-1,:); % note that polygons A and B share a vertex (intentional)
polygonC = polygon(crit_pt+numel(theta2)-1:end,:); % note that polygons A and B share a vertex (intentional)

polygonC = [polygonC; mug_rim_radius+10 polygonC(end,2)];

x_sample_pts = linspace(0,10,16)';
polygonD = [polygonC(end,:);...
  [polygonC(end,1)+x_sample_pts 10-(2/3)*2.^(-x_sample_pts/5)]; ...
  polygonC(end,1)+10 10; ...
  polygonC(end,1) 10; ... 
  polygonA(1,1) 10; ...
  polygonA(1,:) ...
  ];

plot(polygonA(:,1),polygonA(:,2), 'r-', 'linewidth', 2);
hold on
plot(polygonB(:,1),polygonB(:,2), 'c-', 'linewidth', 2);
plot(polygonC(:,1),polygonC(:,2), 'g-', 'linewidth', 2);
plot(polygonD(:,1),polygonD(:,2), 'b-', 'linewidth', 2);
axis equal

% Single combined polygon used for extrusion
polygon = [polygonA; polygonB(2:end,:); polygonC(2:end,:); polygonD(2:end,:)];
polygon = flipud(polygon);
polygon = [polygon; polygon(1,:)]; % repeat starting vertex for the surface generation below

% Alternative polygon (cross section when one of the spacers is present)
polygonBPrime = [mug_rim_radius + cutout_radius * cos(theta2) cutout_radius * sin(theta2)];
altPolygon = [...
mug_inner_radius + clearance cutout_radius; ...
mug_rim_radius cutout_radius; ...
polygonBPrime; ...
mug_rim_radius + cutout_radius -outer_descent; ...
mug_rim_radius + 10 -outer_descent; ...
[mug_rim_radius+10+x_sample_pts 10-(2/3)*2.^(-x_sample_pts/5)]; ...
mug_rim_radius+10+10 10; ...
mug_rim_radius + 10 10; ...
mug_inner_radius + clearance 10 ...
];

figure
plot([polygon(:,1); polygon(1,1)], [polygon(:,2); polygon(1,2)], 'r-', 'linewidth', 2);
hold on
plot([altPolygon(:,1); altPolygon(1,1)], [altPolygon(:,2); altPolygon(1,2)], 'c-', 'linewidth', 2);
axis equal

outerSurface = @(t, angle_param) [...
    polygon(t,1) .* cos(angle_param .* (pi/num_reps - asin(0.5 * gap_size ./ polygon(t,1)))) ...
    polygon(t,2) ...
    polygon(t,1) .* sin(angle_param .* (pi/num_reps - asin(0.5 * gap_size ./ polygon(t,1))))];

vertices = generate_cylindrical_surface(outerSurface, ...
  1, size(polygon,1), size(polygon,1)-1, -1, 1, 24);
% Orientation flip
temp = vertices(1:3:end,:);
vertices(1:3:end,:) = vertices(2:3:end,:);
vertices(2:3:end,:) = temp;
  
verticesA = vertices;

altPolygon = [altPolygon; altPolygon(1,:)]; % repeat starting vertex for the surface generation below

outerAltSurface = @(t, angle_param) [...
    altPolygon(t,1) .* cos(pi/num_reps + angle_param .* asin(0.5 * gap_size ./ altPolygon(t,1))) ...
    altPolygon(t,2) ...
    altPolygon(t,1) .* sin(pi/num_reps + angle_param .* asin(0.5 * gap_size ./ altPolygon(t,1)))];

theta_min = pi/num_reps - gap_size/(2*mug_rim_radius);
theta_max = pi/num_reps + gap_size/(2*mug_rim_radius);
vertices = generate_cylindrical_surface(outerAltSurface, ...
  1, size(altPolygon,1), size(altPolygon,1)-1, -1, 1, 12);
  
verticesB = vertices;

polygonsBAndBPrime = [polygonB polygonBPrime];
bToBPrimeSurface = @(t, theta) [...
    polygonsBAndBPrime(sub2ind(size(polygonsBAndBPrime),theta,2*t-1)) ...
    polygonsBAndBPrime(sub2ind(size(polygonsBAndBPrime),theta,2*t)) ...
    zeros(numel(theta),1)];
    
verticesFan1 = generate_quad_surface(bToBPrimeSurface, ...
  1, 2, 1, 1, size(polygonsBAndBPrime,1), size(polygonsBAndBPrime,1)-1);

polygonAExt = [polygonA; polygonA(end,1) cutout_radius];
fanSurface = @(t) [polygonAExt(t,:) zeros(numel(t),1)];
verticesFan2 = generate_fan(fanSurface, 1, size(polygonAExt,1), size(polygonAExt,1)-1, ...
  [polygonAExt(1,1) cutout_radius 0]);
  
verticesRemainingFill = [...
  mug_rim_radius + cutout_radius 0 0; ...
  mug_rim_radius + mug_rim_minor_radius_h + outer_radius_increase -outer_descent 0; ...
  mug_rim_radius + cutout_radius -outer_descent 0; ...
  mug_rim_radius + cutout_radius 0 0; ...
  mug_rim_radius + mug_rim_minor_radius_h 0 0; ...
  mug_rim_radius + mug_rim_minor_radius_h + outer_radius_increase -outer_descent 0; ...
  ];

% we must apply a custom rotation angle to each point that is dependent on x
% pi/num_reps - asin(0.5 * gap_size ./ verticesFan1(:,1))
% TODO: There is perhaps some way to vectorize this for-loop? Using tensor math, maybe?
verticesDivider = [verticesFan1; verticesFan2; verticesRemainingFill];
verticesDividerA = verticesDivider;
verticesDividerB = verticesDivider;
for i=1:size(verticesDivider,1)
  rotation_angle = pi/num_reps - asin(0.5 * gap_size ./ verticesDivider(i,1));
  rotation_matrix = [...
    cos(rotation_angle) 0 -sin(rotation_angle); ...
    0                   1 0; ...
    sin(rotation_angle) 0  cos(rotation_angle)];
    verticesDividerA(i,:) = verticesDivider(i,:) * rotation_matrix';
    verticesDividerB(i,:) = verticesDivider(i,:) * rotation_matrix;
endfor
verticesDividerB = orientation_flip(verticesDividerB);

vertices = [verticesA; verticesB; verticesDividerA; verticesDividerB];

fid = fopen('mug_rim_holder.stl','wt');
fprintf(fid, 'solid mug_rim_holder\n');
for i=0:num_reps-1
  theta_extra = 2*pi*i/num_reps;
  newVertices = vertices;
  rotation_matrix = [...
    cos(theta_extra) 0 -sin(theta_extra); ...
    0                1 0; ...
    sin(theta_extra) 0  cos(theta_extra)];
  newVertices = newVertices * rotation_matrix';
  for f = 1:(size(newVertices,1)/3)
    fprintf(fid, 'facet normal 0.0 0.0 1.0\n');
    fprintf(fid, '    outer loop\n');
    fprintf(fid, '        vertex %.3f %.3f %.3f\n', newVertices(3*f-2,:));
    fprintf(fid, '        vertex %.3f %.3f %.3f\n', newVertices(3*f-1,:));
    fprintf(fid, '        vertex %.3f %.3f %.3f\n', newVertices(3*f,:));
    fprintf(fid, '    endloop\n');
    fprintf(fid, 'endfacet\n');
  endfor
endfor
fprintf(fid, 'endsolid mug_rim_holder\n');
fclose(fid)



