num_divisions = 12;
gap_size = 3.5;
large_radius = 38;
small_radius = 1.5;
edge_offset = 10;

fid = fopen('mug_rim.stl','wt');
fprintf(fid, 'solid mug_rim\n');

% half torus
gap_theta = gap_size / (large_radius - 3*small_radius);
theta = linspace(3*gap_theta/2, 2*pi/num_divisions - 3*gap_theta/2, 960/num_divisions);

thetaCC = theta(1:end-1);
theta = theta(2:end);
p = linspace(0,pi,21);
r_values = [...
  linspace(large_radius + 3*small_radius, large_radius + 6*small_radius/5, 5) ...
  large_radius + small_radius * cos(p) ...
  linspace(large_radius - 6*small_radius/5, large_radius - 3*small_radius, 5) ...
  large_radius + 3*small_radius * cos(pi - p(1:end-1)) ...
  ];
z_values = [...
  zeros(1,5) ...
  small_radius * sin(p) ...
  zeros(1,5) ...
  3*small_radius * sin(p(1:end-1))
  ];
r_valuesCC = [r_values(end) r_values(1:end-1)];
z_valuesCC = [z_values(end) z_values(1:end-1)];

[thetaGrid1, rGrid1] = ndgrid(theta, r_values);
[thetaGrid1, zGrid1] = ndgrid(theta, z_values);
[thetaGrid2, rGrid2] = ndgrid(theta, r_valuesCC);
[thetaGrid2, zGrid2] = ndgrid(theta, z_valuesCC);
[thetaGrid3, rGrid3] = ndgrid(thetaCC, r_valuesCC);
[thetaGrid3, zGrid3] = ndgrid(thetaCC, z_valuesCC);
[thetaGrid4, rGrid4] = ndgrid(thetaCC, r_values);
[thetaGrid4, zGrid4] = ndgrid(thetaCC, z_values);

verticesA = zeros(6 * numel(thetaGrid1), 3);
verticesA(1:6:end,:) = [...
    rGrid1(:) .* cos(thetaGrid1(:)) ...
    rGrid1(:) .* sin(thetaGrid1(:)) ...
    zGrid1(:)];
verticesA(2:6:end,:) = [...
    rGrid2(:) .* cos(thetaGrid2(:)) ...
    rGrid2(:) .* sin(thetaGrid2(:)) ...
    zGrid2(:)];
verticesA(3:6:end,:) = [...
    rGrid3(:) .* cos(thetaGrid3(:)) ...
    rGrid3(:) .* sin(thetaGrid3(:)) ...
    zGrid3(:)];
verticesA(4:6:end,:) = [...
    rGrid3(:) .* cos(thetaGrid3(:)) ...
    rGrid3(:) .* sin(thetaGrid3(:)) ...
    zGrid3(:)];
verticesA(5:6:end,:) = [...
    rGrid4(:) .* cos(thetaGrid4(:)) ...
    rGrid4(:) .* sin(thetaGrid4(:)) ...
    zGrid4(:)];
verticesA(6:6:end,:) = [...
    rGrid1(:) .* cos(thetaGrid1(:)) ...
    rGrid1(:) .* sin(thetaGrid1(:)) ...
    zGrid1(:)];
    
phiA = p(1:end-1);
phiB = p(2:end);
r_B = linspace(6*small_radius/5, 3*small_radius, 5);
r_A = [small_radius r_B(1:end-1)];

[phiGrid1, rGrid1] = ndgrid(phiA, r_A);
[phiGrid2, rGrid2] = ndgrid(phiA, r_B);
[phiGrid3, rGrid3] = ndgrid(phiB, r_B);
[phiGrid4, rGrid4] = ndgrid(phiB, r_A);

verticesB = zeros(6 * numel(phiGrid1), 3);
verticesB(1:6:end,:) = [...
    (large_radius + rGrid1(:) .* cos(phiGrid1(:))) ...
    zeros(numel(rGrid1), 1) ...
    rGrid1(:) .* sin(phiGrid1(:))];
verticesB(2:6:end,:) = [...
    (large_radius + rGrid2(:) .* cos(phiGrid2(:))) ...
    zeros(numel(rGrid2), 1) ...
    rGrid2(:) .* sin(phiGrid2(:))];
verticesB(3:6:end,:) = [...
    (large_radius + rGrid3(:) .* cos(phiGrid3(:))) ...
    zeros(numel(rGrid3), 1) ...
    rGrid3(:) .* sin(phiGrid3(:))];
verticesB(4:6:end,:) = [...
    (large_radius + rGrid3(:) .* cos(phiGrid3(:))) ...
    zeros(numel(rGrid3), 1) ...
    rGrid3(:) .* sin(phiGrid3(:))];
verticesB(5:6:end,:) = [...
    (large_radius + rGrid4(:) .* cos(phiGrid4(:))) ...
    zeros(numel(rGrid4), 1) ...
    rGrid4(:) .* sin(phiGrid4(:))];
verticesB(6:6:end,:) = [...
    (large_radius + rGrid1(:) .* cos(phiGrid1(:))) ...
    zeros(numel(rGrid1), 1) ...
    rGrid1(:) .* sin(phiGrid1(:))];
    
verticesB1 = verticesB;
verticesB1(:,1:2) = verticesB1(:,1:2) * ...
    [ cos(gap_theta/2) -sin(gap_theta/2); ...
      sin(gap_theta/2) cos(gap_theta/2)];
% Orientation flip
temp = verticesB1(1:3:end,:);
verticesB1(1:3:end,:) = verticesB1(2:3:end,:);
verticesB1(2:3:end,:) = temp;
verticesB2 = verticesB;
verticesB2(:,1:2) = verticesB2(:,1:2) * ...
    [ cos(gap_theta/2) sin(gap_theta/2); ...
      -sin(gap_theta/2) cos(gap_theta/2)];
      
theta = linspace(gap_theta/2, 3*gap_theta/2, 6);
thetaCC = theta(1:end-1);
theta = theta(2:end);
p = linspace(0,pi,21);
r_values = [...
  linspace(large_radius + 3*small_radius, large_radius + 6*small_radius/5, 5) ...
  large_radius + small_radius * cos(p) ...
  linspace(large_radius - 6*small_radius/5, large_radius - 3*small_radius, 5) ...
  large_radius + 3*small_radius * cos(pi - p(1:11)) ...
  ];
z_values = [...
  zeros(1,5) ...
  small_radius * sin(p) ...
  zeros(1,5) ...
  3*small_radius * sin(p(1:11))
  ];
r_valuesCC = r_values(1:end-1);
z_valuesCC = z_values(1:end-1);
r_values = r_values(2:end);
z_values = z_values(2:end);

[thetaGrid1, rGrid1] = ndgrid(theta, r_values);
[thetaGrid1, zGrid1] = ndgrid(theta, z_values);
[thetaGrid2, rGrid2] = ndgrid(theta, r_valuesCC);
[thetaGrid2, zGrid2] = ndgrid(theta, z_valuesCC);
[thetaGrid3, rGrid3] = ndgrid(thetaCC, r_valuesCC);
[thetaGrid3, zGrid3] = ndgrid(thetaCC, z_valuesCC);
[thetaGrid4, rGrid4] = ndgrid(thetaCC, r_values);
[thetaGrid4, zGrid4] = ndgrid(thetaCC, z_values);

verticesC = zeros(6 * numel(thetaGrid1), 3);
verticesC(1:6:end,:) = [...
    rGrid1(:) .* cos(thetaGrid1(:)) ...
    rGrid1(:) .* sin(thetaGrid1(:)) ...
    zGrid1(:)];
verticesC(2:6:end,:) = [...
    rGrid2(:) .* cos(thetaGrid2(:)) ...
    rGrid2(:) .* sin(thetaGrid2(:)) ...
    zGrid2(:)];
verticesC(3:6:end,:) = [...
    rGrid3(:) .* cos(thetaGrid3(:)) ...
    rGrid3(:) .* sin(thetaGrid3(:)) ...
    zGrid3(:)];
verticesC(4:6:end,:) = [...
    rGrid3(:) .* cos(thetaGrid3(:)) ...
    rGrid3(:) .* sin(thetaGrid3(:)) ...
    zGrid3(:)];
verticesC(5:6:end,:) = [...
    rGrid4(:) .* cos(thetaGrid4(:)) ...
    rGrid4(:) .* sin(thetaGrid4(:)) ...
    zGrid4(:)];
verticesC(6:6:end,:) = [...
    rGrid1(:) .* cos(thetaGrid1(:)) ...
    rGrid1(:) .* sin(thetaGrid1(:)) ...
    zGrid1(:)];
    
verticesC2 = verticesC;
verticesC2(:,1:2) = verticesC2(:,1:2) * ...
    [ cos(2*gap_theta) -sin(2*gap_theta); ...
      sin(2*gap_theta) cos(2*gap_theta)];
      
r_values = large_radius + 3*small_radius * cos(pi - p(11:end));
r_valuesA = r_values(1:end-1);
r_valuesB = r_values(2:end);
z_values = 3*small_radius * sin(p(11:end));
z_valuesA = z_values(1:end-1);
z_valuesB = z_values(2:end);
r_values2 = large_radius + edge_offset * ones(11,1);
r_values2A = r_values2(1:end-1);
r_values2B = r_values2(2:end);

verticesD = zeros(60, 3);
verticesD(1:6:end,:) = [...
    r_valuesA(:) ...
    zeros(numel(r_valuesA), 1) ...
    z_valuesA(:)];
verticesD(2:6:end,:) = [...
    r_values2A(:) ...
    zeros(numel(r_values2A), 1) ...
    z_valuesA(:)];
verticesD(3:6:end,:) = [...
    r_values2B(:) ...
    zeros(numel(r_values2B), 1) ...
    z_valuesB(:)];
verticesD(4:6:end,:) = [...
    r_values2B(:) ...
    zeros(numel(r_values2B), 1) ...
    z_valuesB(:)];
verticesD(5:6:end,:) = [...
    r_valuesB(:) ...
    zeros(numel(r_valuesB), 1) ...
    z_valuesB(:)];
verticesD(6:6:end,:) = [...
    r_valuesA(:) ...
    zeros(numel(r_valuesA), 1) ...
    z_valuesA(:)];
    
verticesD1 = verticesD;
verticesD1(:,1:2) = verticesD1(:,1:2) * ...
    [ cos(gap_theta/2) -sin(gap_theta/2); ...
      sin(gap_theta/2) cos(gap_theta/2)];

verticesD2 = verticesD;
verticesD2(:,1:2) = verticesD2(:,1:2) * ...
    [ cos(gap_theta/2) sin(gap_theta/2); ...
      -sin(gap_theta/2) cos(gap_theta/2)];
      
% Orientation flip
temp = verticesD2(1:3:end,:);
verticesD2(1:3:end,:) = verticesD2(2:3:end,:);
verticesD2(2:3:end,:) = temp;
      
verticesD3 = verticesD;
verticesD3(:,1:2) = verticesD3(:,1:2) * ...
    [ cos(3*gap_theta/2) -sin(3*gap_theta/2); ...
      sin(3*gap_theta/2) cos(3*gap_theta/2)];
      
% Orientation flip
temp = verticesD3(1:3:end,:);
verticesD3(1:3:end,:) = verticesD3(2:3:end,:);
verticesD3(2:3:end,:) = temp;
      
verticesD4 = verticesD;
verticesD4(:,1:2) = verticesD4(:,1:2) * ...
    [ cos(3*gap_theta/2) sin(3*gap_theta/2); ...
      -sin(3*gap_theta/2) cos(3*gap_theta/2)];

radius = (large_radius + edge_offset) * gap_theta/2;
p = linspace(-pi/2,pi/2,21);
theta_values = gap_theta/2 * sin(p);
theta_valuesA = theta_values(1:end-1);
theta_valuesB = theta_values(2:end);
r_values = large_radius + edge_offset + radius * cos(p);
r_valuesA = r_values(1:end-1);
r_valuesB = r_values(2:end);
z_values = 3*small_radius * sin(p(11:end));
z_valuesA = z_values(1:end-1);
z_valuesB = z_values(2:end);

[rGrid1, zGrid1] = ndgrid(r_valuesA, z_valuesA);
[thetaGrid1, zGrid1] = ndgrid(theta_valuesA, z_valuesA);
[rGrid2, zGrid2] = ndgrid(r_valuesA, z_valuesB);
[thetaGrid2, zGrid2] = ndgrid(theta_valuesA, z_valuesB);
[rGrid3, zGrid3] = ndgrid(r_valuesB, z_valuesB);
[thetaGrid3, zGrid3] = ndgrid(theta_valuesB, z_valuesB);
[rGrid4, zGrid4] = ndgrid(r_valuesB, z_valuesA);
[thetaGrid4, zGrid4] = ndgrid(theta_valuesB, z_valuesA);

verticesE = zeros(6 * numel(thetaGrid1), 3);
verticesE(1:6:end,:) = [...
    rGrid1(:) .* cos(thetaGrid1(:)) ...
    rGrid1(:) .* sin(thetaGrid1(:)) ...
    zGrid1(:)];
verticesE(2:6:end,:) = [...
    rGrid2(:) .* cos(thetaGrid2(:)) ...
    rGrid2(:) .* sin(thetaGrid2(:)) ...
    zGrid2(:)];
verticesE(3:6:end,:) = [...
    rGrid3(:) .* cos(thetaGrid3(:)) ...
    rGrid3(:) .* sin(thetaGrid3(:)) ...
    zGrid3(:)];
verticesE(4:6:end,:) = [...
    rGrid3(:) .* cos(thetaGrid3(:)) ...
    rGrid3(:) .* sin(thetaGrid3(:)) ...
    zGrid3(:)];
verticesE(5:6:end,:) = [...
    rGrid4(:) .* cos(thetaGrid4(:)) ...
    rGrid4(:) .* sin(thetaGrid4(:)) ...
    zGrid4(:)];
verticesE(6:6:end,:) = [...
    rGrid1(:) .* cos(thetaGrid1(:)) ...
    rGrid1(:) .* sin(thetaGrid1(:)) ...
    zGrid1(:)];
    
radius = (large_radius + edge_offset) * 3 * gap_theta/2;
p = linspace(-pi/2,pi/2,21);
theta_values = 3 * gap_theta/2 * sin(p);
theta_valuesA = theta_values(1:end-1);
theta_valuesB = theta_values(2:end);
r_values = large_radius + edge_offset + radius * cos(p);
r_valuesA = r_values(1:end-1);
r_valuesB = r_values(2:end);
z_values = 3*small_radius * sin(p(11:end));
z_valuesA = z_values(1:end-1);
z_valuesB = z_values(2:end);

[rGrid1, zGrid1] = ndgrid(r_valuesA, z_valuesA);
[thetaGrid1, zGrid1] = ndgrid(theta_valuesA, z_valuesA);
[rGrid2, zGrid2] = ndgrid(r_valuesA, z_valuesB);
[thetaGrid2, zGrid2] = ndgrid(theta_valuesA, z_valuesB);
[rGrid3, zGrid3] = ndgrid(r_valuesB, z_valuesB);
[thetaGrid3, zGrid3] = ndgrid(theta_valuesB, z_valuesB);
[rGrid4, zGrid4] = ndgrid(r_valuesB, z_valuesA);
[thetaGrid4, zGrid4] = ndgrid(theta_valuesB, z_valuesA);

verticesF = zeros(6 * numel(thetaGrid1), 3);
verticesF(1:6:end,:) = [...
    rGrid1(:) .* cos(thetaGrid1(:)) ...
    rGrid1(:) .* sin(thetaGrid1(:)) ...
    zGrid1(:)];
verticesF(2:6:end,:) = [...
    rGrid2(:) .* cos(thetaGrid2(:)) ...
    rGrid2(:) .* sin(thetaGrid2(:)) ...
    zGrid2(:)];
verticesF(3:6:end,:) = [...
    rGrid3(:) .* cos(thetaGrid3(:)) ...
    rGrid3(:) .* sin(thetaGrid3(:)) ...
    zGrid3(:)];
verticesF(4:6:end,:) = [...
    rGrid3(:) .* cos(thetaGrid3(:)) ...
    rGrid3(:) .* sin(thetaGrid3(:)) ...
    zGrid3(:)];
verticesF(5:6:end,:) = [...
    rGrid4(:) .* cos(thetaGrid4(:)) ...
    rGrid4(:) .* sin(thetaGrid4(:)) ...
    zGrid4(:)];
verticesF(6:6:end,:) = [...
    rGrid1(:) .* cos(thetaGrid1(:)) ...
    rGrid1(:) .* sin(thetaGrid1(:)) ...
    zGrid1(:)];

% Orientation flip
temp = verticesF(1:3:end,:);
verticesF(1:3:end,:) = verticesF(2:3:end,:);
verticesF(2:3:end,:) = temp;

radiusA = gap_theta/2;
radiusB = 3*gap_theta/2;
phi = linspace(-pi/2,pi/2,21);
phiA = phi(1:end-1);
phiB = phi(2:end);

[rGrid1, pGrid1] = ndgrid(radiusA, phiA);
[rGrid2, pGrid2] = ndgrid(radiusA, phiB);
[rGrid3, pGrid3] = ndgrid(radiusB, phiB);
[rGrid4, pGrid4] = ndgrid(radiusB, phiA);

verticesG = zeros(6 * numel(rGrid1), 3);
verticesG(1:6:end,:) = [...
    (large_radius + edge_offset + (large_radius+edge_offset)*rGrid1(:) .* cos(pGrid1(:))) .* cos(rGrid1(:) .* sin(pGrid1(:))) ...
    (large_radius + edge_offset + (large_radius+edge_offset)*rGrid1(:) .* cos(pGrid1(:))) .* sin(rGrid1(:) .* sin(pGrid1(:))) ...
    zeros(numel(pGrid1),1)];
verticesG(2:6:end,:) = [...
    (large_radius + edge_offset + (large_radius+edge_offset)*rGrid2(:) .* cos(pGrid2(:))) .* cos(rGrid2(:) .* sin(pGrid2(:))) ...
    (large_radius + edge_offset + (large_radius+edge_offset)*rGrid2(:) .* cos(pGrid2(:))) .* sin(rGrid2(:) .* sin(pGrid2(:))) ...
    zeros(numel(pGrid2),1)];
verticesG(3:6:end,:) = [...
    (large_radius + edge_offset + (large_radius+edge_offset)*rGrid3(:) .* cos(pGrid3(:))) .* cos(rGrid3(:) .* sin(pGrid3(:))) ...
    (large_radius + edge_offset + (large_radius+edge_offset)*rGrid3(:) .* cos(pGrid3(:))) .* sin(rGrid3(:) .* sin(pGrid3(:))) ...
    zeros(numel(pGrid3),1)];
verticesG(4:6:end,:) = [...
    (large_radius + edge_offset + (large_radius+edge_offset)*rGrid3(:) .* cos(pGrid3(:))) .* cos(rGrid3(:) .* sin(pGrid3(:))) ...
    (large_radius + edge_offset + (large_radius+edge_offset)*rGrid3(:) .* cos(pGrid3(:))) .* sin(rGrid3(:) .* sin(pGrid3(:))) ...
    zeros(numel(pGrid3),1)];
verticesG(5:6:end,:) = [...
    (large_radius + edge_offset + (large_radius+edge_offset)*rGrid4(:) .* cos(pGrid4(:))) .* cos(rGrid4(:) .* sin(pGrid4(:))) ...
    (large_radius + edge_offset + (large_radius+edge_offset)*rGrid4(:) .* cos(pGrid4(:))) .* sin(rGrid4(:) .* sin(pGrid4(:))) ...
    zeros(numel(pGrid4),1)];
verticesG(6:6:end,:) = [...
    (large_radius + edge_offset + (large_radius+edge_offset)*rGrid1(:) .* cos(pGrid1(:))) .* cos(rGrid1(:) .* sin(pGrid1(:))) ...
    (large_radius + edge_offset + (large_radius+edge_offset)*rGrid1(:) .* cos(pGrid1(:))) .* sin(rGrid1(:) .* sin(pGrid1(:))) ...
    zeros(numel(pGrid1),1)];
    
verticesG2 = verticesG;
verticesG2(:,3) = 3*small_radius;
% Orientation flip
temp = verticesG2(1:3:end,:);
verticesG2(1:3:end,:) = verticesG2(2:3:end,:);
verticesG2(2:3:end,:) = temp;
    
vertices = [verticesA; verticesB1; verticesB2; verticesC; verticesC2; ...
            verticesD1; verticesD2; verticesD3; verticesD4; verticesE; ...
            verticesF; verticesG; verticesG2];

for i=0:num_divisions-1
  theta_extra = 2*pi*i/num_divisions;
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

fprintf(fid, 'endsolid mug_rim\n');
fclose(fid)