function create_noise_link(filename)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:
%       create_noise_link(filename)  
%
%  PURPOSE: Creates a soft link in the working directory called "noise.mat"
%           to an existing file named by filename.  The file specified by
%           filename must exist or the fucntion will error out.
%
%  INPUTS:
%     filename = the full path to the existing file that is to be linked to 
%                "noise.mat"
%
%  OUTPUTS:
%     None.
%
%  RESTRICTIONS:
%     None.
%
%  WRITTEN BY:  Brian Sammuli on 6/9/08
%
%  MODIFICATION HISTORY:
%
%
%  SEE ALSO: make_noise_file, make_noise 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create soft link: noise.mat -> filename

disp('Linking noise.mat to new noise file...')
gatools_root = getenv('GATOOLS_ROOT');
link_command = [gatools_root,'/matlab/tok_models/set_noise_file ', filename]
[success, comm_output] = system(link_command)
if (success ~= 0)
    error('Failed to create noise.mat link.')
end
disp('Done with noise')

