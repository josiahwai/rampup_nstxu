% SUMMARY OF TAGS:
% circuit names := {'OH', 'PF1AU', 'PF2U', 'PF3U', 'PF5', 'PF3L', 'PF2L', 'PF1AL'};
%
% (+) OH:     {'OH_P1S_1PV', 'OH_P1S_2PV', ... 
%              'OH_P2S_1PV', 'OH_P2S_2PV', ...
%              'OH_P3S_1PV', 'OH_P3S_2PV', ...
%              'OH_P4S_1PV', 'OH_P4S_2PV', ...
%              'OH_P5S_1PV', 'OH_P5S_2PV', ...
%              'OH_P6S_1PV', 'OH_P6S_2PV', ...
%              'OH_P7S_1PV', 'OH_P7S_2PV', ...
%              'OH_P8S_1PV', 'OH_P8S_2PV'};
%               && IOH
%               k1 = [1 3 9 11];     % same sign as IOH
%               k2 = [6 8 10 12];    % opposite IOH and k1
%                
% (+) PF1AU: PF1AU_P1SV, PF1AU_P2SV, PF1AU_P2SI && IPF1AU
% (+) PF2U: PF2U_P1SV, PF2U_P2SV,  PF2U_CUR2 && IPF2U
% (-) PF3U: PF3U_P1S_1PV, PF3U_P1S_2PV, PF3U_P2S_1PV, PF3U_P2S_2PV, PF3U_CUR2 && IPF3U 
% (-) PF5: PF5_P1SV, PF5_P2SV, PF5_P3SV, PF5_CUR1 && IPF5
% (-) PF3L: PF3L_P1S_1PV, PF3L_P1S_2PV, PF3L_P2S_1PV, PF3L_P2S_2PV, PF3L_CUR2 && IPF3L
% (+) PF2L: PF2L_P1SV, PF2L_P2SV,  PF2L_CUR2 && IPF2L
% (+) PF1AL: PF1AL_P1SV, PF1AL_P2SV, PF1AL_CUR1 && IPF1AL

% note for PF3U/L: the 2PV voltages are opposite 1PV voltages. sign 2PV goes
%  with sign IPF3U/L

function coils = mds_fetch_current_voltage(shot, plotit)

% setup connection
MDSPLUS=getenv('MDSPLUS_DIR');
addpath(genpath(MDSPLUS))
mdshost = 'skylark.pppl.gov:8501';
mdsclose;
mdsdisconnect;
mdsconnect(mdshost);

% Load sim
sim = load('sim_inputs204660.mat');
sim_t = sim.sim_inputs.tspan(2:end);
sim_v = sim.sim_inputs.Uhat;

% ===============================================
% Load current and voltages for each coil circuit
% ===============================================
coil.names = {'OH', 'PF1AU', 'PF2U', 'PF3U', 'PF5', 'PF3L', 'PF2L', 'PF1AL'};

dig = '.EPICS.FCPC.DIGITIZERS:';
tree = 'ENGINEERING';

% =======
% OH coil
% =======
ohtags = strcat(dig, ...
         {'OH_P1S_1PV', 'OH_P1S_2PV', ...
          'OH_P2S_1PV', 'OH_P2S_2PV', ...
          'OH_P3S_1PV', 'OH_P3S_2PV', ...
          'OH_P4S_1PV', 'OH_P4S_2PV', ...
          'OH_P5S_1PV', 'OH_P5S_2PV', ...
          'OH_P6S_1PV', 'OH_P6S_2PV', ...
          'OH_P7S_1PV', 'OH_P7S_2PV', ...
          'OH_P8S_1PV', 'OH_P8S_2PV'});
 
k1 = [1 3 9 11];     % same sign as IOH
k2 = [6 8 10 12];    % opposite IOH and k1
[voh1, voha, t] = fetch_and_prep(tree, shot, ohtags(k1), []);
[voh2, vohb] = fetch_and_prep(tree, shot, ohtags(k2), t);
voh = (voh1 - voh2) / 2;

coils.t = t;
coils.v(1,:) = voh;
coils.ic(1,:) = fetch_and_prep(tree, shot, '.ANALYSIS:IOH', t);

if plotit
  figure
  ax(1) = subplot(211);
  hold on  
  plot(t, voha, 'k')
  plot(t, -vohb, 'k')
  plot(t, voh, 'r', 'linewidth', 3)
  plot(sim_t, sim_v(1,:), 'b', 'linewidth', 2)
  ylabel('Voltage')
  
  ax(2) = subplot(212);
  hold on
  plot(t, coils.ic(1,:))
  ylabel('Current')
  sgtitle('OH')  
  linkaxes(ax, 'x')
  xlim([0 0.5])
end

% ==========
% PF1AU coil
% ==========

tags = strcat(dig, {'PF1AU_P2SV', 'PF1AU_P2SV'});
[v_avg, v] = fetch_and_prep(tree, shot, tags, t);
I = fetch_and_prep(tree, shot, '.ANALYSIS:IPF1AU', t);

coils.v(2,:) = v_avg;
coils.ic(2,:) = I;

if plotit
  figure
  ax(1) = subplot(211);
  hold on  
  plot(t, v_avg, 'r', 'linewidth', 2)
  plot(sim_t, sim_v(2,:), 'b', 'linewidth', 2)  
  plot(t, v, 'k')    
  plot(sim_t, sim_v', 'color', [1 1 1] * 0.8)
  plot(t, v_avg, 'r', 'linewidth', 2)
  plot(sim_t, sim_v(2,:), 'b', 'linewidth', 2)
  legend('Power supply', 'Simulation')
  xlabel('Time')
  ylabel('Voltage')
  title('PF1AU')
  ax(2) = subplot(212);
  hold on
  plot(t,I)
  ylabel('Current')
  xlabel('Time')
  linkaxes(ax,'x')
  xlim([0 0.5])
end


% ==========
% PF2U coil
% ==========
tags = strcat(dig, {'PF2U_P1SV', 'PF2U_P2SV'});
[v_avg, v] = fetch_and_prep(tree, shot, tags, t);
I = fetch_and_prep(tree, shot, '.ANALYSIS:IPF2U', t);

coils.v(3,:) = v_avg;
coils.ic(3,:) = I;

if plotit
  figure
  ax(1) = subplot(211);
  hold on  
  plot(t, v_avg, 'r', 'linewidth', 2)
  plot(sim_t, sim_v(3,:), 'b', 'linewidth', 2)  
  plot(t, v, 'k')    
  plot(sim_t, sim_v', 'color', [1 1 1] * 0.8)
  plot(t, v_avg, 'r', 'linewidth', 2)
  plot(sim_t, sim_v(3,:), 'b', 'linewidth', 2)
  legend('Power supply', 'Simulation')
  xlabel('Time')
  ylabel('Voltage')
  title('PF2U')
  ax(2) = subplot(212);
  hold on
  plot(t,I)
  ylabel('Current')
  xlabel('Time')
  linkaxes(ax,'x')
  xlim([0 0.5])
end

%%
% ==========
% PF3U coil
% ==========
tags = strcat(dig, {'PF3U_P1S_1PV', 'PF3U_P1S_2PV', 'PF3U_P2S_1PV', 'PF3U_P2S_2PV'});
[v_avg1, v1] = fetch_and_prep(tree, shot, tags([1 3]), t);
[v_avg2, v2] = fetch_and_prep(tree, shot, tags([2 4]), t);
v_avg = (v_avg2 - v_avg1) / 2;

I = fetch_and_prep(tree, shot, '.ANALYSIS:IPF3U', t);

coils.v(4,:) = v_avg;
coils.ic(4,:) = I;

if plotit
  figure
  ax(1) = subplot(211);
  hold on  
  plot(t, v_avg, 'r', 'linewidth', 2)
  plot(sim_t, sim_v(4,:), 'b', 'linewidth', 2)  
  plot(t, v2, 'k')
  plot(sim_t, sim_v', 'color', [1 1 1] * 0.8)
  plot(t, v_avg, 'r', 'linewidth', 2)
  plot(sim_t, sim_v(4,:), 'b', 'linewidth', 2)
  legend('Power supply', 'Simulation')
  xlabel('Time')
  ylabel('Voltage')
  title('PF3U')
  ax(2) = subplot(212);
  hold on
  plot(t,I)
  ylabel('Current')
  xlabel('Time')
  linkaxes(ax,'x')
  xlim([0 0.5])
end

%% =========
% PF5 coil
% ==========
tags = strcat(dig, {'PF5_P1SV', 'PF5_P2SV', 'PF5_P3SV'});
[v_avg, v] = fetch_and_prep(tree, shot, tags, t);
I = fetch_and_prep(tree, shot, '.ANALYSIS:IPF5', t);
v_avg = -v_avg;

coils.v(5,:) = v_avg;
coils.ic(5,:) = I;

if plotit
  figure
  ax(1) = subplot(211);
  hold on  
  plot(t, v_avg, 'r', 'linewidth', 2)
  plot(sim_t, sim_v(5,:), 'b', 'linewidth', 2)  
  plot(t, -v, 'k')
  plot(sim_t, sim_v', 'color', [1 1 1] * 0.8)
  plot(t, v_avg, 'r', 'linewidth', 2)
  plot(sim_t, sim_v(5,:), 'b', 'linewidth', 2)
  legend('Power supply', 'Simulation')
  xlabel('Time')
  ylabel('Voltage')
  title('PF5')
  ax(2) = subplot(212);
  hold on
  plot(t,I)
  ylabel('Current')
  xlabel('Time')
  linkaxes(ax,'x')
  xlim([0 0.5])
end

%% =========
% PF3L coil
% ==========
tags = strcat(dig, {'PF3L_P1S_1PV', 'PF3L_P1S_2PV', 'PF3L_P2S_1PV', 'PF3L_P2S_2PV'});
[v_avg1, v1] = fetch_and_prep(tree, shot, tags([1 3]), t);
[v_avg2, v2] = fetch_and_prep(tree, shot, tags([2 4]), t);
v_avg = (v_avg2 - v_avg1) / 2;

I = fetch_and_prep(tree, shot, '.ANALYSIS:IPF3L', t);

coils.v(6,:) = v_avg;
coils.ic(6,:) = I;

if plotit
  figure
  ax(1) = subplot(211);
  hold on  
  plot(t, v_avg, 'r', 'linewidth', 2)
  plot(sim_t, sim_v(6,:), 'b', 'linewidth', 2)  
  plot(sim_t, sim_v', 'color', [1 1 1] * 0.8)
  plot(t, v_avg, 'r', 'linewidth', 2)
  plot(sim_t, sim_v(6,:), 'b', 'linewidth', 2)
  legend('Power supply', 'Simulation')
  xlabel('Time')
  ylabel('Voltage')
  title('PF3L')
  ax(2) = subplot(212);
  hold on
  plot(t,I)
  ylabel('Current')
  xlabel('Time')
  linkaxes(ax,'x')
  xlim([0 0.5])
end

%% =========
% PF2L coil
% ==========
tags = strcat(dig, {'PF2L_P1SV', 'PF2L_P2SV'});
[v_avg, v] = fetch_and_prep(tree, shot, tags, t);
I = fetch_and_prep(tree, shot, '.ANALYSIS:IPF2L', t);

coils.v(7,:) = v_avg;
coils.ic(7,:) = I;

if plotit
  figure
  ax(1) = subplot(211);
  hold on  
  plot(t, v_avg, 'r', 'linewidth', 2)
  plot(sim_t, sim_v(7,:), 'b', 'linewidth', 2)  
  plot(t, v, 'k')    
  plot(sim_t, sim_v', 'color', [1 1 1] * 0.8)
  plot(t, v_avg, 'r', 'linewidth', 2)
  plot(sim_t, sim_v(7,:), 'b', 'linewidth', 2)
  legend('Power supply', 'Simulation')
  xlabel('Time')
  ylabel('Voltage')
  title('PF2L')
  ax(2) = subplot(212);
  hold on
  plot(t,I)
  ylabel('Current')
  xlabel('Time')
  linkaxes(ax,'x')
  xlim([0 0.5])
end


%% =========
% PF1AL coil
% ==========
tags = strcat(dig, {'PF1AL_P2SV', 'PF1AL_P2SV'});
[v_avg, v] = fetch_and_prep(tree, shot, tags, t);
I = fetch_and_prep(tree, shot, '.ANALYSIS:IPF1AL', t);

coils.v(8,:) = v_avg;
coils.ic(8,:) = I;

if plotit
  figure
  ax(1) = subplot(211);
  hold on  
  plot(t, v_avg, 'r', 'linewidth', 2)
  plot(sim_t, sim_v(8,:), 'b', 'linewidth', 2)  
  plot(t, v, 'k')    
  plot(sim_t, sim_v', 'color', [1 1 1] * 0.8)
  plot(t, v_avg, 'r', 'linewidth', 2)
  plot(sim_t, sim_v(8,:), 'b', 'linewidth', 2)
  legend('Power supply', 'Simulation')
  xlabel('Time')
  ylabel('Voltage')
  title('PF1AL')
  ax(2) = subplot(212);
  hold on
  plot(t,I)
  ylabel('Current')
  xlabel('Time')
  linkaxes(ax,'x')
  xlim([0 0.5])
end

%% ==============
% Vessel currents
% ===============
tree = 'efit01';
tag = '.RESULTS.AEQDSK:CCBRSP';
mdsopen(tree, shot);
t_iv = mdsvalue(strcat('dim_of(', tag, ')'));
ccefit = mdsvalue(tag);
ccefit = interp1(t_iv, ccefit', t)';
nvess = 40;
iv = ccefit(end-nvess+1:end, :);
coils.iv = iv;
mdsclose;


%% ==============
% Plasma current
% ===============
tag = '.RESULTS.AEQDSK:IPMEAS';
coils.ip = fetch_and_prep(tree, shot, tag, t);


coils.t = double(coils.t);
coils.v = double(coils.v);
coils.ic = double(coils.ic);
coils.iv = double(coils.iv);
coils.ip = double(coils.ip);
end

function [avg_signal, indiv_signals, timebase] = fetch_and_prep(tree, shot, tags, timebase, plotit)

  mdsopen(tree, shot);
  
  if ~exist('plotit','var'), plotit = 0; end
  smoothit = 1;      % applies smoothing filter
  smoothtime = 0.05; % [sec]
  shiftit = 1;       % shift signal so that for time<0, signal=0
  if ~iscell(tags), tags = {tags}; end
  
  if isempty(timebase)
     timebase = mdsvalue(strcat('dim_of(', tags{1}, ')'));
  end  
    
  for i = 1:length(tags)
    t = mdsvalue(strcat('dim_of(', tags{i}, ')'));
    signal = mdsvalue(tags{i});    
    if smoothit
      dt = mean(diff(t));
      signal = smoothdata(signal', 'movmean', floor(smoothtime / dt));
    end    
    if shiftit
      signal = signal - median(signal(t < 0));
    end    
    signal = interp1(t, signal, timebase);    
    indiv_signals(i,:) = signal;
  end
  
  mdsclose;
  
  if length(tags) > 1
    avg_signal = mean(indiv_signals); 
  else
    avg_signal = signal;
  end
  
  if plotit
    figure
    hold on
    plot(timebase, indiv_signals, 'color', [1 1 1] * 0.8)
    plot(timebase, avg_signal, 'b', 'linewidth', 2)
  end
end

















