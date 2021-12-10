% Semi-Free-Boundary Grad-Shafranov solver
% function eq = semifreegs(target, tok_data_struct, opts)

ccc
load('matlab.mat')
opts.max_iterations = 20;

% options arguments for each update step
eq = nan;
pla_opts.cold_start = 1;
pla_opts.time = target.time;
app_opts.use_vessel_currents = 0;
eq_opts.plotit = 0;

if ~exist('opts','var'), opts = struct; end
if ~isfield(opts, 'max_iterations'), opts.max_iterations = 10; end


% iterate for solution
for iter = 1:opts.max_iterations
  pla = update_psi_pla(eq, target, tok_data_struct, pla_opts);
  app = update_psi_app(pla, target, tok_data_struct, app_opts);
  eq = eq_analysis(app.psizr, pla, tok_data_struct, eq_opts);
  pla_opts.cold_start = 0;    
end

eq = copyfields(eq, pla, []);
eq = copyfields(eq, app, []);

if opts.plotit
  semifreegs_plots(eq, opts.init, target, tok_data_struct)
end






























