restoredefaultpath
RAMPROOT = '/Users/jwai/Research/rampup_nstxu/';
addpath(genpath(RAMPROOT));
setenv('RAMPROOT', RAMPROOT);


MDSPLUS_DIR = getenv('MDSPLUS_DIR');
if isempty(MDSPLUS_DIR)
  MDSPLUS_DIR = '/usr/local/mdsplus';
  setenv('MDSPLUS_DIR', MDSPLUS_DIR);
end
addpath(genpath(MDSPLUS_DIR));


addpath(genpath('/Users/jwai/Research/coneqt/data_access'))
addpath('/Users/jwai/Research/coneqt/model_dev')

% Remove from path, any files in a directory labeled 'old'
% (Convenient for development to avoid path conflicts -- one can backup a
% copy with the SAME filename by saving it in any directory called old.) 
d = dir('**/old*');
for i = 1:length(d)
  rmpath([d(i).folder '/' d(i).name]);
end



