
% Smooth a timeseries

function ts2 = smooth_ts(ts, opts)

if ~exist('opts', 'var'), opts = struct; end
if ~isfield(opts, 'p'), opts.p = 0.999; end
if ~isfield(opts, 'remove_outliers'), opts.remove_outliers = 0; end
if ~isfield(opts, 'plotit'), opts.plotit = 0; end
if ~isfield(opts, 'force_start_point'), opts.force_start_point = 0; end
if ~isfield(opts, 'force_end_point'), opts.force_end_point = 0; end


sz = size(ts.Data);
Data = squeeze(reshape(ts.Data, [], length(ts.Time)));
Time = ts.Time;

Data2 = 0*Data;
ts2 = ts;

for i = 1:size(Data,1)
  
  y = Data(i,:);
  t = Time;
  
  if opts.remove_outliers
    [y, idx] = rmoutliers(y, 'movmean', 11);
    t(idx) = [];
  end
  
  % make start/end point more important
  if opts.force_start_point
    N = length(ts.Time);
    t = [t(1)*ones(N,1); t];
    y = [y(1)*ones(N,1); y];
  end 
  if opts.force_end_point
    N = length(ts.Time);
    t = [t; t(end)*ones(N,1)];
    y = [y; y(end)*ones(N,1)];
  end
  
  f = fit(t(:), y(:) , 'smoothingspline', 'SmoothingParam', opts.p);
  Data2(i,:) = f(ts.Time);
end

ts2.Data = reshape(Data2, sz);



if opts.plotit
  plot(ts, '--r')
  hold on
  plot(ts2, 'b')
end
















