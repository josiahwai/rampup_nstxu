
% External circuit L/R Fitting Algorithm
vacuum_system = load('NSTXU_vaccum_system.mat').NSTXU_vacuum_system;

Rxx = vacuum_system.Rxx;
Mxx = vacuum_system.Mxx;

%Rvv = diag(vacuum_system.Rxx(vacuum_system.ncx+1:end,vacuum_system.ncx+1:end));

ncx = vacuum_system.ncx;


%...............................................................................
% Generate the grey-box model for system identification

fit_coils = [13];
odefun = 'coil_dynamics';




%...............................................................................
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
    
    %icts = timeseries(icsignals.sigs,icsignals.times);
    %vobjcts = timeseries(vobjcsignals.sigs,vobjcsignals.times);
    
    % Flux Loops

%     flsignals = get_flsignals(shot, [], [], flsignals_names);
% 
%     flts = timeseries(flsignals.sigs,flsignals.times);
    % fltime = flsignals.times;
    % flsigs = flsignals.sigs;
    % 
    % % Loop Voltages
    % 
    % lvsignals = get_lvsignals(shot, [], [], lvsignals_names);
    % 
    % lvtime = lvsignals.times;
    % lvsigs = lvsignals.sigs;
    % 
    % % Mirnov Probes
    % 
    % bpsignals = get_bpsignals(shot, [], [], bpsignals_names);
    % 
    % bptime = bpsignals.times;
    % bpsigs = bpsignals.sigs;
    % 
    % % Concatenate diagnostic signals
    % 
    % diagsigs = [flsigs bpsigs];
    % 
    % diagtime = (fltime + bptime)/2;
    Ts = mean(diff(ictime));
    Fc = 1000;
    Fs = 1/mean(diff(ictime));
     
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
    %icsigs_filtered = icsigs;
    %vobjcsigs_filtered = vobjcsigs;
%     lpfreq = 180; %Hz
%     ictspass = idealfilter(icts,[lpfreq,inf],'notch');
%     vobjctspass = idealfilter(vobjcts,[0,lpfreq],'pass');
    %icdot = [diff(ictspass.Data); zeros(1,ncx)]/Ts;
    
    % Trim coil data to length of magnetics signals

    %tstart = min(flsignals.times);
    %tend   = max(flsignals.times);
    tstart = starttimes(i);
    tend = endtimes(i);
    tsample = tstart:Ts:tend;

    %flts = resample(flts,tsample);
    icts = timeseries(icsigs_filtered,ictime);
    vobjcts = timeseries(vobjcsigs,vobjctime);
    icts = resample(icts,tsample);
    vobjcts = resample(vobjcts,tsample);
    % idx = find((icsignals.times >= tstart) & (icsignals.times <= tend));
    % 
    % ictime  = icsignals.times([idx; idx(end)+1]);
    % icsigs  = icsignals.sigs( [idx; idx(end)+1],:);
    % 
    % times = ictime;
    % 
    % dt = mean(diff(times));

    %...............................................................................
    % Filter the data

    % Design Butterworth filter

    % Fc = 60;
    % Fs = 1/mean(diff(diagtime));
    % 
    % [b,a] = butter(4, Fc/(Fs/2));
    % 
    % icsigs_filtered = zeros(size(icsigs,1), size(icsigs,2));
    % 
    % for ii = 1:size(icsigs,2)
    %     
    %     icsigs_filtered(:,ii) = filtfilt(b, a, icsigs(:,ii));
    %     
    % end
    % 
    % diagsigs_filtered = zeros(size(diagsigs,1), size(diagsigs,2));
    % 
    % for ii = 1:size(diagsigs,2)
    %     
    %     diagsigs_filtered(:,ii) = filtfilt(b, a, diagsigs(:,ii));
    %     
    % end



    %...............................................................................
    % Determine vessel current contributions to the diagnostic signals

    %yc = Cc*ictspass.Data';

    %yv = fltspass.Data' - yc;

    %...............................................................................
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
    %for i=1:length(shotlist)
    maxy = max(maxy,max(shotdata{i}.OutputData));
    miny = min(miny,min(shotdata{i}.OutputData));
    %end
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
%     Ts = 0.0002;
    sys = idgrey(odefun, parameters, fcn_type,[],Ts,'InputDelay',3);

    sys.Structure.Parameters(1).Free = false;
    sys.Structure.Parameters(2).Free = false;
    sys.Structure.Parameters(3).Free = false;
    sys.Structure.Parameters(4).Free = false;
            
    sys.Structure.Parameters(5).Minimum = -10*ones(length(fit_coils),1);
    sys.Structure.Parameters(5).Maximum = 500*ones(length(fit_coils),1);
    sys.Structure.Parameters(6).Minimum = -1*ones(length(fit_coils),1);
    sys.Structure.Parameters(6).Maximum = 100*ones(length(fit_coils),1);
    
%     sys.Structure.Parameters(5).Free = true;
%     sys.Structure.Parameters(5).Minimum = -0.1*ones(length(fit_coils),1);
%     sys.Structure.Parameters(5).Maximum = 0.1*ones(length(fit_coils),1);    
%     
%     sys.Structure.Parameters(6).Free = true;
%     sys.Structure.Parameters(6).Minimum = -0.1*ones(length(fit_coils),1);
%     sys.Structure.Parameters(6).Maximum = 0.1*ones(length(fit_coils),1);    

    opt = greyestOptions('Display', 'on', 'InitialState', 'zero', 'DisturbanceModel', 'none', ...
    'Focus', 'simulation', 'SearchMethod', 'fmincon','OutputWeight',diag(1./rangey));
 
    sys_est = greyest(shotdata{i}, sys, opt)
    [y,t,x] = lsim(sys_est,shotdata{i}.InputData,tsample);
    figure
    plot(y)
    hold on
    plot(shotdata{i}.OutputData)
    Lext_fit{i} = sys_est.Structure.Parameters(1,6).Value/1000;
    Rext_fit{i} = sys_est.Structure.Parameters(1,5).Value/1000;
end
%%



%%
%Rvv_fit = sys_est.Structure.Parameters(1,3).Value/1000;

%save('Rvv_fit.mat','Rvv_fit');

nplotrows = 3;
nplotcols = 2;
nplotspertab = nplotrows*nplotcols;
ntabs = ceil(size(Cin,1)/nplotspertab);

for exp = 1:length(shotlist)
    figure()
    tg = uitabgroup; % tabgroup
    j=1;%keep track of which output to plot
    for ii = 1:ntabs
      thistab = uitab(tg,'Title',['Page ' num2str(ii)]); % build iith tab
      axes('Parent',thistab); % somewhere to plot
      for k=1:nplotspertab
        if j<=size(Cin,1)
            subplot(nplotrows,nplotcols,k)
            plot(y{exp}.OutputData(:,j));
            hold on
            exp_data = getexp(all_data,exp);
            plot(exp_data.OutputData(:,j));
            title(flsignals_names(j))
            j=j+1;
        end
      end
    end
end


% Lext_fit = [Lext_fit{:}];
% Rext_fit = [Rext_fit{:}];
% icoil = [fitcoillist{:}];
% ext_fit = variables2struct(icoil, Lext_fit, Rext_fit);
% save('ext_fit', 'ext_fit')











