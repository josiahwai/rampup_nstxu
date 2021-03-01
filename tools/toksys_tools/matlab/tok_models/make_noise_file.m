function make_noise_file(filename, tstart, tend, dt, num_diags, vars, noise_type)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:
%    make_noise_file(filename, tstart, tend, dt, num_diags)
%    make_noise_file(filename, tstart, tend, dt, num_diags, vars)
%    make_noise_file(filename, tstart, tend, dt, num_diags, vars, noise_type)
%  
%
%  PURPOSE: Generate a file containing a noise matrix that can be read into 
%           a "from file" simulink block and create a soft link to that
%           file called "noise.mat".  
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
%               num_diags.
%
%     noise_type = (optional) the type of noise to be generated. If vars is
%                  specified, defaults to 'gaussian'.  Otherwise defaults
%                  to 'zeros' (ie no noise).
%                    Allowed Values:
%                    - 'gaussian' - Noise will be gaussian 
%                    - 'zeros' - No noise 
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
%  SEE ALSO: make_noise, create_noise_link 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (nargin < 5)
    error('To few arguments')
elseif(nargin == 5)
    noise_type = 'zeros';
    vars = zeros(1, num_diags);
elseif (nargin == 6)
    noise_type = 'gaussian';
elseif (nargin == 7)
    if(strcmpi(noise_type, 'gaussian') ~= 1)
        disp('Non-Gaussian distributions not yet supported. Setting type to Gaussian')
        noise_type = 'gaussian'; 
    end
end

t = [tstart:dt:tend];

[mt, nt] = size(t);
noise = zeros(num_diags+1, nt); %num_diags+1 because we must include the time vector
noise(1,:) = t;

if (strcmpi(noise_type, 'gaussian'))
    disp('Creating Gaussian white noise...')
    [mvars, nvars] = size(vars);
    if (nvars > mvars) %Make sure vars is a column vector
        vars = vars';
        [mvars, nvars] = size(vars);
    end
    if(nvars > 1)
        error('vars must be a vector')    
    end

    if (~(mvars == num_diags)) 
        error('num_diags must equal the number of points in vars vector')
    end


    for i=1:mvars
        noise(i+1,:) = sqrt(vars(i))*randn(1, nt);     
    end
else %Default to zero noise
    disp('Defaulting to zero noise...')
end
disp('Saving noise file...')
save(filename, 'noise');

%%% Create soft link: noise.mat -> filename
%disp('Linking noise.mat to new noise file...')
%link_command = ['set_noise_file ', filename];
%[success, comm_output] = system(link_command);
%if (success ~= 0)
%    error('Failed to create noise.mat link.')
%end
%disp('Done with noise')

