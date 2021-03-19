ccc

load('train_shots_times.mat')
shots = double(shots_times.shots);

save_name = 'train_responses';
save_dir = '/Users/jwai/Research/rampup_nstxu/gspert_response/builds/';


tree = 'EFIT01';
tokamak = 'nstxu';
server = 'skylark.pppl.gov:8501';

tok_data_struct = load('nstxu_obj_config2016_6565.mat').tok_data_struct;
circ = nstxu2016_circ(tok_data_struct);


X = [];
validshots = [];
validtimes = [];
validshots_all = [];
validtimes_all = [];
j = 0;
nbatch = 0;

for ishot = 1:length(shots)
  
  shot = shots(ishot);
  true_times = double(shots_times.times{ishot});
  
  % nominal_times = (20:20:200) / 1e3;
  % nominal_times = (250:50:2000) / 1e3;
  nominal_times = [20:20:200  250:50:2000] / 1e3;
  
  k = nominal_times > min(true_times) & nominal_times < max(true_times);
  nominal_times(~k) = [];
  
  ntimes = length(nominal_times);
  
  for itime = 1:ntimes
    time = nominal_times(itime);
    
    str = ['Fetching shot ' num2str(ishot) ' of ' num2str(length(shots)) ': '];
    str = [str [num2str(shot) ' @ ' num2str(floor(time*1000)) 'ms']];
    disp(str);
    
    try
      xi = get_xmat_target(shot, time, tree, tokamak, server, tok_data_struct, circ);
      j = j + 1;
      
      X(:,j) = xi;
      validshots(j) = shot;
      validtimes(j) = time;  
      validshots_all(end+1) = shot;
      validtimes_all(end+1) = shot;
            
      % Save data every 200 samples
      if j >= 200
        nbatch = nbatch + 1;
        
        response = variables2struct(X, validshots, validtimes);
        save([save_dir save_name num2str(nbatch) '.mat'], 'response')
        
        X = [];
        validshots = [];
        validtimes = [];
        j = 0;
      end
      
    catch
      warning(['Could not obtain ' num2str(shot) ' ' num2str(time)])
    end
  end      
end

% save the last partial batch
nbatch = nbatch + 1;
response = variables2struct(X, validshots, validtimes, validshots_all, validtimes_all);
save([save_dir save_name num2str(nbatch) '.mat'], 'response')





%%
function y = get_xmat_target(shot, time, tree, tokamak, server, tok_data_struct, circ)

eq = read_eq(shot,time,tree,tokamak,server);

response = gspert(eq, tok_data_struct);

mpc = tok_data_struct.mpc;
mpv = tok_data_struct.mpv;

dcphidis = [response.dcphidis response.dcphidip(:)];
dcphidix = dcphidis * circ.Pxx;

% xmat_cv = circ.Pxx(1:end-1, 1:end-1)' * [mpc'; mpv'] * dcphidix;
y = dcphidix(:);

%   xmats = response.xmats;
%   xmatp = [mpc'; mpv'] * reshape(response.dcphidip,[], 1);
%   Ldum = 0;
%   xmatx = [xmats xmatp; xmatp' Ldum];
%   xmat = circ.Pxx' * xmatx * circ.Pxx;
%   xmat_cv = xmat(1:end-1, :);
%   y = reshape(xmat_cv, [], 1);

end

