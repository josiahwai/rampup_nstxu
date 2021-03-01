function make_noise(filename, tstart, tend, dt, num_diags, vars)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:
%    make_noise(filename, tstart, tend, dt, num_diags)
%    make_noise(filename, tstart, tend, dt, num_diags, vars)
%  
%
%  PURPOSE: Generate a file containing a noise matrix that can be read into 
%           a "from file" simulink block and create a soft link to that
%           file called "noise.mat" in the working directory.  
%
%           The matrix is formatted such that the first row is a vector of
%           times.  Each subsequent row is a time series of noise that can
%           be added to a diagnostic output signal.  Row i of the matrix 
%           maps to the i-1 diagnostic.
%
%  INPUTS:
%     filename = the full path to the .mat file that will be created
%
%     tstart = the starting time of the time vector
%
%     tend = the end time of the time vector
%
%     dt = the time step of the time vector
%
%     num_diags = the number of diagnostic outputs
%
%     vars = (optional) a vector of variances (sigma^2).  Must be of length
%               num_diags. If vars is not used, then a zero-noise matrix is
%               generated.
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
%  SEE ALSO: make_noise_file, create_noise_link 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin < 5)
    error('To few arguments.')
elseif  (nargin == 5)
    make_noise_file(filename, tstart, tend, dt, num_diags)
elseif (nargin == 6)
    make_noise_file(filename, tstart, tend, dt, num_diags, vars)
end

create_noise_link(filename)
