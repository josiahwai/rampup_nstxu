function vec = pop(vec, val, idx)

if size(vec,1) == 1, vec = vec'; end

vec = [vec(1:idx-1,:); val; vec(idx:end,:)];


