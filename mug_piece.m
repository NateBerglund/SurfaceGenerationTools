clear all
close all
clc

extrusion_height = 15.4;

polygon = csvread('mug_piece_border_polygon.csv');
polygon(:,2) = 20-polygon(:,2); % compensate for the original being in image coords

polygon = flipud(polygon); % reverse the order (makes it make more sense visually)
polygon = [polygon; polygon(1,:)]; % repeat the first vertex

figure
hold on
plot(polygon(:,1),polygon(:,2), 'r-', 'linewidth', 2);
plot(polygon(1,1),polygon(1,2), 'bo', 'markersize', 10);
axis equal

size(polygon)

polygonFan = @(t) [...
    polygon(mod(t+size(polygon,1)-1,size(polygon,1))+1,:) ...
    zeros(numel(t),1)];
firstFanVertex = [13 5 0];
%verticesA = generate_fan(polygonFan, 1, 1, 1-1, firstFanVertex);
secondFanVertex = [13 16 0];
verticesB = generate_fan(polygonFan, 1, 370, 370-1, secondFanVertex);
verticesC = generate_fan(polygonFan, 370, 613, 613-370, firstFanVertex);

vertices = [verticesB; verticesC; ...
secondFanVertex; firstFanVertex; [polygon(1,:) 0]; ...
firstFanVertex; secondFanVertex; [polygon(370,:) 0]; ...
];

vertices2 = vertices;
vertices2(:,3) = extrusion_height;
vertices = orientation_flip(vertices);

vertices = [vertices; vertices2];

fid = fopen('mug_piece.stl','wt');
fprintf(fid, 'solid extruded_polygon\n');

for f = 1:(size(vertices,1)/3)
  fprintf(fid, 'facet normal 0.0 0.0 1.0\n');
  fprintf(fid, '    outer loop\n');
  fprintf(fid, '        vertex %.3f %.3f %.3f\n', vertices(3*f-2,:));
  fprintf(fid, '        vertex %.3f %.3f %.3f\n', vertices(3*f-1,:));
  fprintf(fid, '        vertex %.3f %.3f %.3f\n', vertices(3*f,:));
  fprintf(fid, '    endloop\n');
  fprintf(fid, 'endfacet\n');
endfor

fprintf(fid, 'endsolid extruded_polygon\n');
fclose(fid)


