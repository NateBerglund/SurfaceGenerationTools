function vertices = orientation_flip(vertices)
  temp = vertices(1:3:end,:);
  vertices(1:3:end,:) = vertices(2:3:end,:);
  vertices(2:3:end,:) = temp;
endfunction