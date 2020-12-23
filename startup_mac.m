restoredefaultpath
root = '/Users/jwai/Research/rampup_nstxu/';
addpath(genpath(root));


% Remove from path, any files in a directory labeled 'old'
% (Convenient for development to avoid path conflicts -- one can backup a
% copy with the same filename by saving it in any directory called old.) 
d = dir('**/old*');
for i = 1:length(d)
  rmpath([d(i).folder '/' d(i).name]);
end



