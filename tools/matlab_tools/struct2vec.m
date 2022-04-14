% makes a vector given a struct and fields to concatenate
function v = struct2vec(s, fields)

v = [];

for i = 1:length(fields)
  fd = fields{i};
  v = [v; s.(fd)(:)];
end



