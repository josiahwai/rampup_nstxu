restoredefaultpath
RAMPROOT = '/Users/jwai/Research/rampup_nstxu/';
addpath(genpath(RAMPROOT));
setenv('RAMPROOT', RAMPROOT);

mdsplus_dir = getenv('MDSPLUS_DIR');
addpath(genpath(mdsplus_dir));
addpath(genpath('/Users/jwai/Research/coneqt/data_access'))

addpath('/Users/jwai/Research/coneqt/model_dev')

% Remove from path, any files in a directory labeled 'old'
% (Convenient for development to avoid path conflicts -- one can backup a
% copy with the same filename by saving it in any directory called old.) 
d = dir('**/old*');
for i = 1:length(d)
  rmpath([d(i).folder '/' d(i).name]);
end



