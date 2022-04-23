function vertices = generate_fan(f, t_min, t_max, n_t, center_pt)
  
  t_params = linspace(t_min, t_max, n_t + 1)';
  t_params_shifted = t_params(2:end);
  t_params = t_params(1:end-1);
  
  vertices = zeros(3 * numel(t_params), 3);
  
  vertices(1:3:end,:) = repmat(center_pt, numel(t_params), 1);
  vertices(2:3:end,:) = f(t_params);
  vertices(3:3:end,:) = f(t_params_shifted);

endfunction