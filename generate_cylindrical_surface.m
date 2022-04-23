function vertices = generate_cylindrical_surface(f, t_min, t_max, n_t, u_min, u_max, n_u)
  
  t_params = linspace(t_min, t_max, n_t + 1)';
  t_params = t_params(1:end-1);
  t_params_shifted = [t_params(end); t_params(1:end-1)];
  u_params = linspace(u_min, u_max, n_u + 1)';
  u_params_shifted = u_params(2:end);
  u_params = u_params(1:end-1);
  
  [tGrid1, uGrid1] = ndgrid(t_params, u_params);
  [tGrid2, uGrid2] = ndgrid(t_params, u_params_shifted);
  [tGrid3, uGrid3] = ndgrid(t_params_shifted, u_params_shifted);
  [tGrid4, uGrid4] = ndgrid(t_params_shifted, u_params);
  
  vertices = zeros(6 * numel(tGrid1), 3);
  vertices(1:6:end,:) = f(tGrid1(:), uGrid1(:));
  vertices(2:6:end,:) = f(tGrid2(:), uGrid2(:));
  vertices(3:6:end,:) = f(tGrid3(:), uGrid3(:));
  vertices(4:6:end,:) = f(tGrid3(:), uGrid3(:));
  vertices(5:6:end,:) = f(tGrid4(:), uGrid4(:));
  vertices(6:6:end,:) = f(tGrid1(:), uGrid1(:));

endfunction