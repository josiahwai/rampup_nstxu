% Copies the fields of B to the fields of A. If overwrite=1, the fields of 
% B will overwrite contents of corresponding field in A. If overwrite=0, 
% only fields of B that do not exist in A will be copied. 

function A = copyfields(A,B,overwrite)

if nargin == 2, overwrite = false; end

fns = fieldnames(B);
for i = 1:length(fns)
  fn = fns{i};  
  if ~isfield(A,fn) || overwrite
    A.(fn) = B.(fn);
  end
end
