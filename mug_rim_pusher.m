clear all
close all
clc

num_reps = 12;
gap_size = 3.5;
small_radius = 4.5;
large_radius = 7.5;
xy_radius = 45/2 - small_radius;
extra = xy_radius * tan(pi / num_reps) - 3*gap_size/2;
offset = xy_radius - large_radius - (gap_size/2) / tan(pi / num_reps);

resolution = 0.25;

innerCircle = @(theta) [...
    xy_radius + small_radius * cos(theta) ...
    zeros(size(theta)) ...
    small_radius * sin(theta)];

innerCylinder = @(theta, y) [...
    xy_radius + small_radius * cos(theta) ...
    y ...
    small_radius * sin(theta)];
    
outerCylinder = @(theta, y) [...
    xy_radius + large_radius * cos(theta) ...
    y ...
    large_radius * sin(theta)];
    
annulus = @(theta, r) [...
    xy_radius + r .* cos(theta) ...
    zeros(size(theta)) ...
    r .* sin(theta)];
    
xy_rect = @(x, y) [...
    x ...
    y ...
    zeros(size(x))];
    
xz_rect = @(x, z) [...
    x ...
    zeros(size(x)) ...
    z];
    
yz_rect = @(y, z) [...
    zeros(size(y)) ...
    y ...
    z];

fid = fopen('mug_rim_pusher.stl','wt');
fprintf(fid, 'solid mug_rim_pusher\n');

% inner cylinder
n_theta_divs = 4 * ceil(0.5 * pi * large_radius / resolution);
n_y_divs = 2 * ceil(0.5 * gap_size / resolution);
verticesA = generate_quad_surface(innerCylinder, ...
  0, pi/2, n_theta_divs/4, ...
  -gap_size/2, gap_size/2, n_y_divs);
verticesA = [verticesA;...
  xy_radius-large_radius-offset  gap_size/2 large_radius;...
  xy_radius-large_radius-offset -gap_size/2 large_radius;...
  xy_radius                      gap_size/2 small_radius;...
  xy_radius-large_radius-offset -gap_size/2 large_radius;...
  xy_radius                     -gap_size/2 small_radius;...
  xy_radius                      gap_size/2 small_radius;...
  xy_radius-large_radius-offset  gap_size/2 large_radius;...
  xy_radius                      gap_size/2 small_radius;...
  xy_radius                      gap_size/2 large_radius;...
  xy_radius-large_radius-offset -gap_size/2 large_radius;...
  xy_radius                     -gap_size/2 large_radius;...
  xy_radius                     -gap_size/2 small_radius;...
  xy_radius-large_radius-offset  gap_size/2 large_radius;...
  xy_radius                      gap_size/2 large_radius;...
  xy_radius                      (3*gap_size/2+extra) large_radius;...
  xy_radius                     -gap_size/2 large_radius;...
  xy_radius-large_radius-offset -gap_size/2 large_radius;...
  xy_radius                     -(3*gap_size/2+extra) large_radius;...
  xy_radius-large_radius-offset -gap_size/2 large_radius;...
  xy_radius-large_radius-offset  gap_size/2 large_radius;...
  0                              0          large_radius;...
  ];
verticesA2 = verticesA;
verticesA2(:,3) = -verticesA2(:,3);
verticesA2 = orientation_flip(verticesA2);
verticesA = [verticesA; verticesA2];

% outer cylinder
verticesB = generate_quad_surface(outerCylinder, ...
  -pi/2, pi/2, n_theta_divs/2, ...
  gap_size/2, (3*gap_size/2+extra), n_y_divs);
verticesB2 = verticesB;
verticesB2(:,2) = -verticesB2(:,2);
verticesB2 = orientation_flip(verticesB2);

verticesB = [verticesB; verticesB2];

% annuli
n_r_divs = ceil((large_radius - small_radius) / resolution);
verticesC = generate_quad_surface(annulus, ...
  -pi/2, pi/2, n_theta_divs/2, ...
  small_radius, large_radius, n_r_divs);
verticesC1 = verticesC;
verticesC1(:,2) = gap_size/2;
##verticesC2 = verticesC;
##verticesC2(:,2) = -(3*gap_size/2+extra);
verticesC3 = verticesC;
verticesC3(:,2) = -gap_size/2;
verticesC3 = orientation_flip(verticesC3);
##verticesC4 = verticesC;
##verticesC4(:,2) = (3*gap_size/2+extra);
##verticesC4 = orientation_flip(verticesC4);

verticesC = [verticesC1; verticesC3];

##% half-circles
##verticesD = generate_fan(innerCircle, -pi/2, pi/2, n_theta_divs/2, [xy_radius 0 0]);
##verticesD2 = verticesD;
##verticesD2(:,2) = -(3*gap_size/2+extra);
##verticesD3 = verticesD;
##verticesD4 = verticesD;
##verticesD4(:,2) = (3*gap_size/2+extra);
##verticesD4 = orientation_flip(verticesD4);
##verticesD = [verticesD2; verticesD4];

theta = linspace(-pi/2, pi/2, n_theta_divs/2 + 1)';
pts1 = [...
    xy_radius + large_radius * cos(theta) ...
    (3*gap_size/2+extra) * ones(size(theta)) ...
    large_radius * sin(theta)];
pts2 = pts1;
pts2(:,2) = -(3*gap_size/2+extra);
pts2(:,1:2) = pts2(:,1:2) * ...
    [ cos(2*pi/num_reps) sin(2*pi/num_reps); ...
     -sin(2*pi/num_reps) cos(2*pi/num_reps)];
pts1A = pts1(1:end-1,:);
pts1B = pts1(2:end,:);
pts2A = pts2(1:end-1,:);
pts2B = pts2(2:end,:);
verticesD = zeros(6 * size(pts1A,1),3);
verticesD(1:6:end,:) = pts1A;
verticesD(2:6:end,:) = pts1B;
verticesD(3:6:end,:) = pts2B;
verticesD(4:6:end,:) = pts2B;
verticesD(5:6:end,:) = pts2A;
verticesD(6:6:end,:) = pts1A;
verticesD = orientation_flip(verticesD);

##% fan fillers
##verticesE = generate_fan(innerCircle, pi/2, pi, n_theta_divs/2, [xy_radius-small_radius 0 small_radius]);
##verticesE1 = verticesE;
##verticesE1(:,2) = gap_size/2;
##verticesE1 = orientation_flip(verticesE1);
##verticesE2 = verticesE;
##verticesE2(:,2) = -gap_size/2;
##verticesE1 = [verticesE1; verticesE2];
##verticesE2 = verticesE1;
##verticesE2(:,3) = -verticesE2(:,3);
##verticesE2 = orientation_flip(verticesE2);
##
##verticesE = [verticesE1; verticesE2];

##% rectangles
##verticesF = generate_quad_surface(xz_rect, ...
##  xy_radius-small_radius, xy_radius, 1, ...
##  small_radius, large_radius, 1);
####verticesF = [verticesF; ...
####  generate_quad_surface(xz_rect, ...
####  xy_radius-large_radius-offset, xy_radius-small_radius, 1, ...
####  0, large_radius, 1)];
##verticesF(:,2) = -gap_size/2;
##
##verticesF = [verticesF;
##  xy_radius-large_radius-offset gap_size/2 large_radius;
##  xy_radius gap_size/2 large_radius
##  xy_radius (3*gap_size/2+extra) large_radius
##];
##
####temp_rect = generate_quad_surface(yz_rect, ...
####  gap_size/2, (3*gap_size/2+extra), 1, ...
####  0, large_radius, 1);
####temp_rect(:,1) = xy_radius-large_radius-offset;
####verticesF = [verticesF; temp_rect];
##
##verticesF1 = verticesF;
##verticesF2 = verticesF;
##
##verticesF2 = verticesF1;
##verticesF2(:,2) = -verticesF2(:,2);
##verticesF2 = orientation_flip(verticesF2);
##
##verticesF1 = [verticesF1; verticesF2];
##verticesF2 = verticesF1;
##verticesF2(:,3) = -verticesF2(:,3);
##verticesF2 = orientation_flip(verticesF2);
##
##verticesF = [verticesF1; verticesF2];

vertices = [verticesA; verticesB; verticesC; verticesD];

for i=0:num_reps-1
  theta_extra = 2*pi*i/num_reps;
  newVertices = vertices;
  newVertices(:,1:2) = newVertices(:,1:2) * ...
    [cos(theta_extra) -sin(theta_extra); ...
     sin(theta_extra) cos(theta_extra)];
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

fprintf(fid, 'endsolid mug_rim_pusher\n');
fclose(fid)