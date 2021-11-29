function [targets, targets_array, efit01_eqs] = read_target(shot, times, tok_data_struct)

ROOT = getenv('RAMPROOT');

tree = 'EFIT01';
tokamak = 'nstxu';
server = 'skylark.pppl.gov:8501';
opts.cache_dir = [ROOT '/fetch/cache/'];
efit01_eqs = fetch_eqs_nstxu(shot, times, tree, tokamak, server, opts);

circ = nstxu2016_circ(tok_data_struct);
nr = tok_data_struct.nr;
nz = tok_data_struct.nz;
rg = tok_data_struct.rg;
zg = tok_data_struct.zg;


% targets are: desired boundary, boundary-defining pt, and Ip
N = length(efit01_eqs.time);

% target boundary
gap_opts.plotit = 0;
for i = 1:N
  gaps(i) = get_nstxu_gaps(efit01_eqs.gdata(i), gap_opts);
end
targets.rcp = [gaps(:).r]';
targets.zcp = [gaps(:).z]';
ngaps = size(targets.rcp, 2);


% boundary defining point
bry_opts.plotit = 0;
for i = 1:N
  bry(i) = eq_bdef_analysis(efit01_eqs.gdata(i), tok_data_struct, bry_opts);
end
targets.rbdef = [bry(:).rbdef]';
targets.zbdef = [bry(:).zbdef]';
targets.islimited = [bry(:).islimited]';

% coil and vessel currents 
targets.icx = [efit01_eqs.gdata(:).icx]';
targets.ivx = [efit01_eqs.gdata(:).ivx]';
targets.icx_select = targets.icx(:,circ.ikeep);

% Ip
targets.ip = [efit01_eqs.gdata(:).cpasma]';
targets.pres = [efit01_eqs.gdata(:).pres]';
targets.fpol = [efit01_eqs.gdata(:).fpol]';
targets.pprime = [efit01_eqs.gdata(:).pprime]';
targets.ffprim = [efit01_eqs.gdata(:).ffprim]';


% Wmhd
targets.wmhd = read_wmhd(efit01_eqs, tok_data_struct);

% li
for i = 1:N
  [~,~,~,targets.li(i,1)] = inductance(efit01_eqs.gdata(i), tok_data_struct);
end

% put onto original requested timebase
N = length(times); 
fns = fieldnames(targets);
for i = 1:length(fns)
  targets.(fns{i})(1:N,:) = interp1(efit01_eqs.time, targets.(fns{i}), times);
end
[~,i] = min(abs(efit01_eqs.time(:) - times(:)'));
efit01_eqs.time = efit01_eqs.time(i);
efit01_eqs.adata = efit01_eqs.adata(i);
efit01_eqs.gdata = efit01_eqs.gdata(i);
efit01_eqs.tms = efit01_eqs.tms(i);


% rearrange data for easier access later
for i = 1:N
  for j = 1:length(fns)
    targets_array(i).(fns{j}) = targets.(fns{j})(i,:);
  end
end

targets.time = times;







