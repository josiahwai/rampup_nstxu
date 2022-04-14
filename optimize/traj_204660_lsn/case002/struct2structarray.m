% takes a struct with fields each of size [N x A x B x etc]  and rearranges 
% to an [N x 1] struct array with fields each of size [A x B x etc]


function s2 = struct2structarray(s)

fds = fieldnames(s);
N = size(s.(fds{1}), 1);

for j = 1:length(fds)
  fd = fds{j};
  for i = 1:N  
    s2(i).(fd) = squeeze(s.(fd)(i,:,:,:));
  end
end




