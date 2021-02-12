restoredefaultpath
RAMPROOT = '/u/jwai/rampup_nstxu/';
addpath(genpath(RAMPROOT));
setenv('RAMPROOT', RAMPROOT);


% Remove from path, any files in a directory labeled 'old'
% (Convenient for development to avoid path conflicts -- one can backup a
% copy with the same filename by saving it in any directory called old.) 
cd(RAMPROOT)
d = dir('**/old*');
for i = 1:length(d)
  rmpath([d(i).folder '/' d(i).name]);
end



