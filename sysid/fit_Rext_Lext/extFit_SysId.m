% External circuit L/R Fitting Algorithm

% Fits resistance (Rext) and inductance (Lext) of power supplies using data 
% from a dedicated set of experiments. 
%
% This is edited from the version in Dan Boyer's "coneqt" repository. 
%
% Author: Josiah Wai 11/2021

ROOT = getenv('RAMPROOT');
saveit = 1;
savedir = [ROOT '/sysid/fit_Rext_Lext'];


load('nstxu_obj_config2016_6565.mat');
circ = nstxu2016_circ(tok_data_struct);
Pxx = circ.Pxx(1:end-1,1:end-1); % remove Ip circuit from transition map
Rxx = Pxx' * diag([tok_data_struct.resc; tok_data_struct.resv]) * Pxx;
Mxx = Pxx' * [tok_data_struct.mcc tok_data_struct.mcv; tok_data_struct.mcv' tok_data_struct.mvv] * Pxx;
ncx = circ.ncx;


%......................................................
% Generate the grey-box model for system identification
odefun = 'SystemDynamics';


%..............
% Load the data

% Coil Currents
include_coils = {'OH', 'PF1aU', 'PF1bU', 'PF1cU', 'PF2U', 'PF3U', 'PF4', ...
        'PF5', 'PF3L', 'PF2L', 'PF1cL', 'PF1bL', 'PF1aL'};
    
%%   
shotlist = [202135,202145,202142,202125,202118,202283,202085,202073,202092,202114];
starttimes = [0,0,0,0,0,0,0,0,0,0];
endtimes = [0.6,1.5,1.5,0.5,0.5,1.5,0.5,0.5,0.5,0.5];
fitcoillist = {[10],[9],[11],[5],[8],[1],[13],[2],[6],[7]};
i=0;

for shot=shotlist
    i=i+1;
    icsignals = get_icsignals(shot, [], [], include_coils);
    vobjcsignals = get_vobjcsignals(shot, [], [], include_coils);

    icsigs = icsignals.sigs;
    ictime = icsignals.times;
    
    vobjcsigs = vobjcsignals.sigs;
    vobjctime = vobjcsignals.times;
    
    Ts = mean(diff(ictime)); 
    Fc = 1000;
    Fs = 1 / Ts;
     
    [b,a] = butter(4, Fc/(Fs/2));
     
    icsigs_filtered = zeros(size(icsigs,1), size(icsigs,2));
    
    for ii = 1:size(icsigs,2)
         
        icsigs_filtered(:,ii) = filtfilt(b, a, icsigs(:,ii));
         
    end
    
    Fs = 1/mean(diff(vobjctime));     
    [b,a] = butter(4, Fc/(Fs/2));     
    vobjcsigs_filtered = zeros(size(vobjcsigs,1), size(vobjcsigs,2));
    
    for ii = 1:size(vobjcsigs,2)         
        vobjcsigs_filtered(:,ii) = filtfilt(b, a, vobjcsigs(:,ii));         
    end
    
    tstart = starttimes(i);
    tend = endtimes(i);
    tsample = tstart:Ts:tend;

    
    icts = timeseries(icsigs_filtered,ictime);
    vobjcts = timeseries(vobjcsigs,vobjctime);
    icts = resample(icts,tsample);
    vobjcts = resample(vobjcts,tsample);
    
    
    %...................
    % Identify the model
    fit_coils = fitcoillist{i};
    shotdata{i} = iddata(double(icts.Data(:,fit_coils)), double(vobjcts.Data(:,fit_coils)), Ts);
    if i==1
       all_data = shotdata{i};
    else
       all_data = merge(all_data,shotdata{i});
    end

    maxy = 0;
    miny = 0;    
    maxy = max(maxy,max(shotdata{i}.OutputData));
    miny = min(miny,min(shotdata{i}.OutputData));    
    rangey = maxy-miny;
    
    Rxx_use = diag(diag(Rxx));
    Rxx_use(1:ncx,1:ncx) = 10*eye(ncx);
    Rxx_use(fit_coils,fit_coils) = Rxx(fit_coils,fit_coils);
    
    Rext_mOhm = 40;
    Lext_mH = 4;

    parameters = {'Mxx', Mxx; 'Rxx', Rxx_use;...
        'fit_coils',fit_coils; 'ncx', ncx; 'Rext_mOhm',Rext_mOhm; ...
        'Lext_mH',Lext_mH; };

    fcn_type = 'd';

    sys = idgrey(odefun, parameters, fcn_type,[],Ts,'InputDelay',3);

    sys.Structure.Parameters(1).Free = false;
    sys.Structure.Parameters(2).Free = false;
    sys.Structure.Parameters(3).Free = false;
    sys.Structure.Parameters(4).Free = false;
            
    sys.Structure.Parameters(5).Minimum = -10*ones(length(fit_coils),1);
    sys.Structure.Parameters(5).Maximum = 500*ones(length(fit_coils),1);
    sys.Structure.Parameters(6).Minimum = -1*ones(length(fit_coils),1);
    sys.Structure.Parameters(6).Maximum = 100*ones(length(fit_coils),1);    

    opt = greyestOptions('Display', 'on', 'InitialState', 'zero', 'DisturbanceModel', 'none', ...
    'Focus', 'simulation', 'SearchMethod', 'fmincon','OutputWeight',diag(1./rangey));
 
    sys_est = greyest(shotdata{i}, sys, opt);
    [y,t,x] = lsim(sys_est,shotdata{i}.InputData,tsample);
    figure
    plot(y)
    hold on
    plot(shotdata{i}.OutputData)
    Lext_fit(i) = sys_est.Structure.Parameters(1,6).Value/1000;
    Rext_fit(i) = sys_est.Structure.Parameters(1,5).Value/1000;
end
%%

% ORGANIZE AND SAVE FITS

% reorder the fits to correspond with standard ordering and naming
idx_fits = [fitcoillist{:}];
% load('nstxu_obj_config2016_6565.mat')
% circ = nstxu2016_circ(tok_data_struct);

for i = 1:length(circ.ikeep)
  [~,k] = find(idx_fits == circ.ikeep(i));
  if isempty(k), k = nan; end
  isort(i) = k;
end

% coils that were used during campaign, see circ.keep_coils for names
Rext_keep_coils = Rext_fit(isort)'; 
Lext_keep_coils = Lext_fit(isort)'; 

% assign high resistance to unused coils
Rext = 1000*ones(circ.ncx,1); 
Rext(circ.ikeep) = Rext_keep_coils;

Lext = median(Lext_keep_coils)*ones(circ.ncx,1); 
Lext(circ.ikeep) = Lext_keep_coils;

disp(' ')
disp('Fitted values for resistance R [Ohm] and inductance L [Henry]')
disp('Note: PF1BU/L, PF1CU/L, PF4 not used, given high resistance.') 
disp(' ')
for i = 1:circ.ncx
  fprintf('%s: R=%.04f, L=%.04f \n', circ.ccnames{i}, Rext(i), Lext(i))
end

ccnames = circ.ccnames;
info = 'Power supply resistances [Ohms] and inductances [Henry] fit via system id (extFit_SysId.m).'; 
ext_fit = variables2struct(Rext, Lext, ccnames, info);

if saveit  
  save([savedir '/ext_fit.mat'], 'ext_fit')  
end




















