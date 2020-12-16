%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gsdesign_demo_kstar_ISS
%
%  PURPOSE: Design ITER-similar shape for KSTAR
%
%  INPUTS:  None (just run the code)
%
%  OUTPUTS: Plots of and tables of feedforward currents
%
%  METHOD:  
	
%
%  WRITTEN BY:  Anders Welander  ON	4/8/14
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Target separatrix from ITER_data_2010-v3.3.xls
iter.rsep = [5.5513, 5.5171, 5.5083, 5.4646, 5.4205, 5.4195, 5.376, 5.3312, 5.322, 5.2864, 5.2417, 5.1171, 5.4195, 5.5171, 5.539, 5.6146, 5.6989, 5.7122, 5.8098, 5.8531, 5.9073, 5.9967, 6.0049, 6.1024, 6.1404, 6.2, 6.2758, 6.2976, 6.3951, 6.4112, 6.4927, 6.5376, 6.5902, 6.6588, 6.6878, 6.7762, 6.7854, 6.8829, 6.8924, 6.9805, 6.9992, 7.078, 7.0994, 7.1756, 7.1938, 7.2732, 7.282, 7.3621, 7.3707, 7.4413, 7.4683, 7.5167, 7.5659, 7.5874, 7.6507, 7.6634, 7.7129, 7.761, 7.7713, 7.8228, 7.8585, 7.8727, 7.9168, 7.9561, 7.9594, 7.9959, 8.0291, 8.0537, 8.0613, 8.088, 8.1112, 8.1313, 8.1493, 8.1512, 8.1672, 8.1802, 8.1897, 8.1963, 8.2001, 8.2011, 8.1993, 8.1948, 8.1876, 8.1775, 8.1642, 8.1512, 8.1461, 8.1282, 8.1082, 8.0855, 8.0596, 8.0537, 8.0285, 7.9968, 7.9622, 7.9561, 7.9222, 7.881, 7.8585, 7.8349, 7.7876, 7.761, 7.7353, 7.6814, 7.6634, 7.6212, 7.5659, 7.5572, 7.4919, 7.4683, 7.4191, 7.3707, 7.3407, 7.2732, 7.2572, 7.1756, 7.1684, 7.078, 7.074, 6.9805, 6.9735, 6.8829, 6.8662, 6.7854, 6.7503, 6.6878, 6.6235, 6.5902, 6.4927, 6.4785, 6.3951, 6.3201, 6.2976, 6.2, 6.1293, 6.1024, 6.0049, 5.9073, 5.8761, 5.8098, 5.7122, 5.6146, 5.5171, 5.4195, 5.322, 5.2244, 5.1268, 5.0875, 5.0293, 4.9317, 4.9213, 4.8341, 4.8101, 4.7366, 4.7258, 4.6561, 4.639, 4.5984, 4.5491, 4.5415, 4.5063, 4.4686, 4.4439, 4.4359, 4.4048, 4.3773, 4.353, 4.3463, 4.3316, 4.3114, 4.2932, 4.2771, 4.263, 4.2506, 4.2488, 4.24, 4.2297, 4.2207, 4.2133, 4.2073, 4.2025, 4.1989, 4.1964, 4.1948, 4.1943, 4.1948, 4.1962, 4.1986, 4.202, 4.2064, 4.2117, 4.2179, 4.225, 4.2332, 4.2423, 4.2488, 4.252, 4.2631, 4.2755, 4.289, 4.3035, 4.3188, 4.335, 4.3463, 4.3519, 4.3698, 4.3888, 4.4087, 4.4295, 4.4439, 4.4507, 4.4727, 4.4955, 4.5192, 4.5415, 4.5434, 4.5679, 4.593, 4.6188, 4.639, 4.6447, 4.6706, 4.6971, 4.7239, 4.7366, 4.75, 4.7762, 4.8028, 4.8301, 4.8341, 4.856, 4.8825, 4.9104, 4.9317, 4.9388, 4.9681, 4.9997, 5.1171, 4.8341, 4.7366, 4.6783, 4.639, 4.5415, 4.4492, 4.4439, 4.3463, 4.2488, 4.2266];
iter.zsep = [-4.3754, -4.2976, -4.2778, -4.1802, -4.0827, -4.0805, -3.9851, -3.8876, -3.8674, -3.79, -3.6924, -3.4139, -3.2712, -3.2175, -3.2046, -3.1597, -3.1071, -3.0987, -3.0382, -3.0095, -2.9731, -2.912, -2.9064, -2.8411, -2.8144, -2.7718, -2.7168, -2.701, -2.6313, -2.6193, -2.557, -2.5217, -2.4798, -2.4241, -2.4004, -2.3266, -2.3189, -2.2374, -2.229, -2.149, -2.1315, -2.0552, -2.0339, -1.9556, -1.9363, -1.8489, -1.8388, -1.7412, -1.7308, -1.6437, -1.6093, -1.5461, -1.4791, -1.4485, -1.351, -1.3314, -1.2534, -1.1738, -1.1559, -1.0583, -0.9892, -0.9607, -0.8632, -0.7735, -0.7656, -0.6681, -0.5705, -0.4972, -0.4729, -0.3754, -0.2778, -0.1802, -0.0827, -0.0729, 0.0149, 0.1124, 0.21, 0.3076, 0.4051, 0.5027, 0.6002, 0.6978, 0.7954, 0.8929, 0.9905, 1.0619, 1.088, 1.1856, 1.2832, 1.3807, 1.4783, 1.4976, 1.5759, 1.6734, 1.771, 1.7864, 1.8685, 1.9661, 2.0143, 2.0637, 2.1612, 2.2116, 2.2588, 2.3563, 2.3863, 2.4539, 2.5385, 2.5515, 2.649, 2.6816, 2.7466, 2.8075, 2.8441, 2.9235, 2.9417, 3.0315, 3.0393, 3.1328, 3.1368, 3.2279, 3.2344, 3.3172, 3.332, 3.4009, 3.4295, 3.4788, 3.5271, 3.551, 3.6154, 3.6246, 3.6776, 3.7222, 3.735, 3.7854, 3.8198, 3.8322, 3.8722, 3.9068, 3.9173, 3.9379, 3.9611, 3.9774, 3.9867, 3.9882, 3.981, 3.9635, 3.9334, 3.9173, 3.8892, 3.8274, 3.8198, 3.746, 3.7222, 3.6383, 3.6246, 3.5271, 3.4998, 3.4295, 3.332, 3.3154, 3.2344, 3.1368, 3.0639, 3.0393, 2.9417, 2.8441, 2.7466, 2.7169, 2.649, 2.5515, 2.4539, 2.3563, 2.2588, 2.1612, 2.1449, 2.0637, 1.9661, 1.8685, 1.771, 1.6734, 1.5759, 1.4783, 1.3807, 1.2832, 1.1856, 1.088, 0.9905, 0.8929, 0.7954, 0.6978, 0.6002, 0.5027, 0.4051, 0.3076, 0.21, 0.1438, 0.1124, 0.0149, -0.0827, -0.1802, -0.2778, -0.3754, -0.4729, -0.5389, -0.5705, -0.6681, -0.7656, -0.8632, -0.9607, -1.0273, -1.0583, -1.1559, -1.2534, -1.351, -1.4408, -1.4485, -1.5461, -1.6437, -1.7412, -1.8174, -1.8388, -1.9363, -2.0339, -2.1315, -2.1784, -2.229, -2.3266, -2.4241, -2.5217, -2.5366, -2.6193, -2.7168, -2.8144, -2.8876, -2.912, -3.0095, -3.1071, -3.4139, -3.5303, -3.5705, -3.5949, -3.6114, -3.6528, -3.6924, -3.6947, -3.7372, -3.7801, -3.79];
iter.rx = +5.1171;
iter.zx = -3.4139;
ix = min(find(iter.zsep == iter.zx));

% Create the ITER-similar-shape
% Settings to make outer boundary run parallel to limiter
scale_factor = 0.24;
r_offset = 0.3;
z_offset = -0.13;
iss.rsep = r_offset + scale_factor*iter.rsep;
iss.zsep = z_offset + scale_factor*iter.zsep;
iss.rx = r_offset + scale_factor*iter.rx;
iss.zx = z_offset + scale_factor*iter.zx;

% Load KSTAR tokamak configuration
if exist('kstar_obj_filename','var')
  try
    clear tok_data_struct
    load(kstar_obj_filename)
    tok_data_struct.rg;
  catch
    error('kstar_obj_filename does not hold name of a matfile containing tok_data_struct')
  end
else
  kstar_obj_filename = [getenv('GATOOLS_ROOT'), ...
    '/tokamaks/kstar/make/kstar_obj_mks_struct_6565.mat'];
  if exist(kstar_obj_filename,'file')
    load(kstar_obj_filename)
  else
    kstar_obj_filename = ...
      '/m/GAtools/tokamaks/kstar/make/kstar_obj_mks_struct_6565.mat';
    if exist(kstar_obj_filename,'file')
      load(kstar_obj_filename)
    else
      disp('Could not find kstar objects. Please set the variable:')
      disp('kstar_obj_filename')
      disp('to name of matfile containing tok_data_struct for kstar')
    end
  end
end

config = tok_data_struct;
config.constraints = 1;
config.psikn = [0 0.40 0.70 1];
config.no_edge_current = true;
config.no_edge_gradient = true;
init = [];

clear spec gsdesign

% Specify where the poloidal field should equal zero
spec.targets.rx = iss.rx;
spec.targets.zx = iss.zx;
spec.weights.x = 1000;

% Specify what point defines the boundary
spec.targets.rbdef = iss.rx;
spec.targets.zbdef = iss.zx;

% Specify where the flux should equal the boundary flux
spec.targets.rsep = iss.rsep;
spec.targets.zsep = iss.zsep;

% Specify plasma current, li, betap, boundary flux
spec.targets.cpasma = round(150*scale_factor^2/2)*1e5;
spec.weights.cpasma = 1e-3;
spec.targets.li = 1;
spec.weights.li = 10;
spec.targets.betap = 0.65;
spec.weights.betap = 5;
spec.targets.psibry = 0;
spec.weights.psibry = 5;

% Specify connections
spec.cccirc = [1 2 3 4 5 6 7 1 2 8 9 10 11 7 13 12 -13 12];

% Lock radial and vertical control coils to 0
spec.locks.ic = nan(18,1); % nan means don't lock
spec.locks.ic(15:18) = 0;

% Call gsdesign with these specs
eq = gsdesign(spec, init, config);

return

disp('Generating table of solutions for different psibry...')
psibrys = -2:0.5:2;
report = 'psibry  ';
for i = 1:config.nc
  report = [report config.ccnames(i,:) '   '];
end
for i = 1:length(psibrys)
  spec.targets.psibry = psibrys(i);
  eq = gsdesign(spec);
  report(i+1,1:7) = sprintf('%7.4f',psibrys(i));
  for j = 1:config.nc
     report(i+1,j*7-1:j*7+5) = sprintf('%7d',round(eq.ic(j)));
  end
end
disp(report)
