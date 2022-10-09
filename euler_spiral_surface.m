clear all
close all
clc

offset = 0.9; % minimum wall thickness recommended by PrusaSlicer is 0.86
subdivs = 512;

polygon = csvread('euler_spiral_points.csv');

% Find the first location where the tangent is horizontal
tangents = diff(polygon,1,1);
idxH = find(diff(sign(tangents(:,2))) ~= 0,1,'last') + 1
y_min = polygon(idxH,2) - offset
x_max = polygon(end,1) + offset

polygon2 = [-1 polygon(1,2); polygon; polygon(end,1) 1];
tangents = diff(polygon2,1,1);
normals = [tangents(:,2) -tangents(:,1)];
normals = bsxfun(@rdivide, normals, sqrt(sum(normals.^2,2))); % make unit length
normals = 0.5 * (normals(1:end-1,:) + normals(2:end,:)); % make symmetric
normals = bsxfun(@rdivide, normals, sqrt(sum(normals.^2,2))); % make unit length

normals(1,:) = [0 -1];
normals(end,:) = [1 0];

polygon2 = polygon2(2:end-1,:) + offset * normals;

polygon3 = polygon2;
multipliersA = x_max ./ polygon3(idxH:end,1);
multipliersB = y_min ./ polygon3(idxH:end,2);
multipliers = min(multipliersA, multipliersB);

idxOuterShell = idxH + find(multipliersA < multipliersB,1,'first') - 1

size(multipliers)
size(polygon3(idxH:end,:))
polygon3(idxH:end,:) = bsxfun(@times, multipliers, polygon3(idxH:end,:));
%polygon3(idxH:end-5,:) = y_min * bsxfun(@rdivide, polygon3(idxH:end-5,:), )

figure
hold on
plot(polygon(:,1),polygon(:,2), 'r-', 'linewidth', 2);
plot(polygon(1,1),polygon(1,2), 'bo', 'markersize', 10);
plot(polygon2(:,1),polygon2(:,2), 'b-', 'linewidth', 2);
plot(polygon(idxH,1),polygon(idxH,2), 'go', 'markersize', 10);
plot(polygon3(:,1),polygon3(:,2), 'g-', 'linewidth', 2);
axis equal

size(polygon)

polygonSurface = @(t,u) [...
    polygon(mod(u+size(polygon,1)-1,size(polygon,1))+1,1) .* sin(t) ...
    polygon(mod(u+size(polygon,1)-1,size(polygon,1))+1,2) ...
    polygon(mod(u+size(polygon,1)-1,size(polygon,1))+1,1) .* cos(t) ...
    ];
fanSurface = @(t) [...
    polygon(2,1) .* sin(t) ...
    polygon(2,2) .* ones(size(t)) ...
    polygon(2,1) .* cos(t) ...
    ];
    
polygon2Surface = @(t,u) [...
    polygon2(mod(u+size(polygon2,1)-1,size(polygon2,1))+1,1) .* sin(t) .* ((t <= 3 * pi / 4) | (t >= 5 * pi / 4)) + ...
    polygon3(mod(u+size(polygon3,1)-1,size(polygon3,1))+1,1) .* sin(t) .* ((t > 3 * pi / 4) & (t < 5 * pi / 4)) ...
    polygon2(mod(u+size(polygon2,1)-1,size(polygon2,1))+1,2) .* ((t <= 3 * pi / 4) | (t >= 5 * pi / 4)) + ...
    polygon3(mod(u+size(polygon3,1)-1,size(polygon3,1))+1,2) .* ((t > 3 * pi / 4) & (t < 5 * pi / 4)) ...
    polygon2(mod(u+size(polygon2,1)-1,size(polygon2,1))+1,1) .* cos(t) .* ((t <= 3 * pi / 4) | (t >= 5 * pi / 4)) + ...
    (polygon3(mod(u+size(polygon3,1)-1,size(polygon3,1))+1,1) .* cos(t) .* ((mod(u+size(polygon3,1)-1,size(polygon3,1))+1)<idxOuterShell) - ...
    polygon3(mod(u+size(polygon3,1)-1,size(polygon3,1))+1,1) .* ((mod(u+size(polygon3,1)-1,size(polygon3,1))+1)>=idxOuterShell)) .* ((t > 3 * pi / 4) & (t < 5 * pi / 4)) ...
    ];
fan2Surface = @(t) [...
    polygon2(2,1) .* sin(t) ...
    polygon2(2,2) .* ones(size(t)) ...
    polygon2(2,1) .* cos(t) ...
    ];

vertices = [generate_cylindrical_surface(polygonSurface, ...
  0, 2 * pi, subdivs, ...
  2, size(polygon,1), size(polygon,1)-2)];
verticesF = generate_fan(fanSurface, 0, 2 * pi, subdivs, [0 polygon(1,2) 0]);
verticesF = orientation_flip(verticesF);
vertices = [vertices; verticesF];
  
vertices2 = [generate_cylindrical_surface(polygon2Surface, ...
  0, 2 * pi, subdivs, ...
  2, size(polygon2,1), size(polygon2,1)-2)];
%vertices2 = orientation_flip(vertices2);
  
verticesF = generate_fan(fan2Surface, 0, 2 * pi, subdivs, [0 polygon2(1,2) 0]);
verticesF = orientation_flip(verticesF);
vertices2 = [vertices2; verticesF];

##vertices2B = [generate_quad_surface(polygon3Surface, ...
##  3 * pi / 4, 5 * pi / 4, subdivs / 4, ...
##  2, size(polygon3,1), size(polygon3,1)-2)];
##vertices2B = orientation_flip(vertices2B);
##vertices2 = [vertices2; vertices2B];

vertices = orientation_flip(vertices);
vertices = [vertices; vertices2];
  
verticesB = vertices;
verticesB(:,2) = -verticesB(:,2);
verticesB = orientation_flip(verticesB);

vertices = [vertices; verticesB];

fid = fopen('euler_spiral_surface.stl','wt');
fprintf(fid, 'solid euler_spiral_surface\n');

for f = 1:(size(vertices,1)/3)
  fprintf(fid, 'facet normal 0.0 0.0 1.0\n');
  fprintf(fid, '    outer loop\n');
  fprintf(fid, '        vertex %.3f %.3f %.3f\n', vertices(3*f-2,:));
  fprintf(fid, '        vertex %.3f %.3f %.3f\n', vertices(3*f-1,:));
  fprintf(fid, '        vertex %.3f %.3f %.3f\n', vertices(3*f,:));
  fprintf(fid, '    endloop\n');
  fprintf(fid, 'endfacet\n');
endfor

fprintf(fid, 'endsolid euler_spiral_surface\n');
fclose(fid)