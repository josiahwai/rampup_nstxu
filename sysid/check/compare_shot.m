ccc
vacuum_system = load('NSTXU_vaccum_system.mat').NSTXU_vacuum_system;
ncx = vacuum_system.ncx;
Rxx = vacuum_system.Rxx;
Mxx = vacuum_system.Mxx;

load('nstxu_obj_config2016_6565.mat');
circ = nstxu2016_circ(tok_data_struct);

% INSERTING RESULTS FROM VESSEL FITTING HERE
% Mvc = load('sys_est_opt3.mat').sys_est.Structure.Parameters(2).Value;
% Mxx(circ.iicx, circ.iivx) = Mvc';
% Mxx(circ.iivx, circ.iicx) = Mvc;
% Rvv = load('sys_est_opt1.mat').sys_est.Structure.Parameters(4).Value;
% Rxx(circ.iivx,circ.iivx) = diag(Rvv/1000);
% 
% 
% Mcc = Mxx(circ.iicx, circ.iicx);
% Mvv = Mxx(circ.iivx, circ.iivx);
% % load('Mvc.mat')
% % load('Mvv.mat')
% Mss = [Mcc Mvc'; Mvc Mvv];
% A = -inv(Mss)*Rxx;
% e = esort(eig(A)); e(1:5)
% 
% A = -inv(Mxx)*Rxx;
% e = esort(eig(A)); e(1:5)

% rxx = diag(Rxx);
% % rvv = load('Rvv_fit.mat').Rvv_fit;
% rvv = diag(circ.Pvv' * diag(tok_data_struct.resv) * circ.Pvv);
% rxx(ncx+1:end) = rvv;
% Rxx = diag(rxx);

odefun = 'coil_dynamics';

%..............
% Load the data

include_coils = {'OH', 'PF1aU', 'PF1bU', 'PF1cU', 'PF2U', 'PF3U', 'PF4', ...
        'PF5', 'PF3L', 'PF2L', 'PF1cL', 'PF1bL', 'PF1aL'};
    
% shotlist = [202135,202145,202142,202125,202118,202283,202085,202073,202092,202114];
% starttimes = [0,0,0,0,0,0,0,0,0,0];
% endtimes = [0.6,1.5,1.5,0.5,0.5,1.5,0.5,0.5,0.5,0.5];
% fitcoillist = {[10],[9],[11],[5],[8],[1],[13],[2],[6],[7]};

shotlist = [204660];
starttimes = [0];
endtimes = [1];
fitcoillist = {[1     2     5     6     8     9    10    13]};
fitcoilnames = include_coils(fitcoillist{1});


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
    
    tstart = starttimes(i);
    tend = endtimes(i);
    tsample = tstart:Ts:tend;

    icts = timeseries(icsigs_filtered,ictime);
    vobjcts = timeseries(vobjcsigs,vobjctime);
    icts = resample(icts,tsample);
    vobjcts = resample(vobjcts,tsample);
        
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
    
    load('ext_fit.mat')
    for j = 1:length(fit_coils)
      k = find(fit_coils(j) == ext_fit.icoil);
      Rext_mOhm(j) = ext_fit.Rext_fit(k) * 1000;
      Lext_mH(j)   = ext_fit.Lext_fit(k) * 1000;
      
%       Rext_mOhm(j) = 0;  % 40
%       Lext_mH(j) = 0; % 4
    end
        
    parameters = {'Mxx', Mxx; 'Rxx', Rxx_use;...
        'fit_coils',fit_coils; 'ncx', ncx; 'Rext_mOhm',Rext_mOhm; ...
        'Lext_mH',Lext_mH; };

    fcn_type = 'd';
    sys = idgrey(odefun, parameters, fcn_type,[],Ts,'InputDelay',3);
            
    load('sim_inputs204660_smoothed.mat')    
    x0 = sim_inputs.traj.x(k,:)';
    x0 = sim_inputs.traj.x(1,:)';
    x0(end) = [];
    
    [y,t,x] = lsim(sys,shotdata{i}.InputData,tsample,x0);
    
    x0 = x(1,:)';
    [y,t,x] = lsim(sys,shotdata{i}.InputData,tsample,x0);
            
    % compensate for unipolar feature of PF1A and PF2    
    searchstrs = {'PF1aU','PF2U','PF2L','PF1aL'};
    for k = 1:length(searchstrs)      
      ii_unipolar(k) = find(strcmp(fitcoilnames, searchstrs{k}));
    end    
        
%     % METHOD 1:
%     y_comp = y;        
%     v = smoothdata(shotdata{i}.InputData,1);    
%     for k = ii_unipolar
%       ioff = find( v(:,k) > 0, 1, 'first') + 100;
%       yshift = shotdata{i}.OutputData(ioff,k) - y(ioff,k);
%       y_comp(:,k) = y_comp(:,k) + yshift;
%     end
%     figure
%     plot(t,v(:,ii_unipolar))    
    
    
%     % METHOD 2:
%     y_comp = y;
%     tt = floor( 0.9 * length(y));
%     yshift = shotdata{i}.OutputData(tt,ii_unipolar) - y(tt,ii_unipolar);
%     y_comp(:,ii_unipolar) = y_comp(:,ii_unipolar) + yshift;      
%     v = smoothdata(shotdata{i}.InputData,1);    
%     for k = ii_unipolar
%       ioff = find( v(:,k) > 0, 1, 'first');
%       y_comp(1:ioff,k) = 0;            
%     end    
%     figure
%     plot(t,v(:,ii_unipolar))    
        
    
    % NO Compensation
    y_comp = y;

    
    figure
    hold on
    plot(t,y_comp(:,1),'b')
    plot(t, shotdata{i}.OutputData(:,1), '--r')           
    plot(t,y_comp,'b')
    plot(t, shotdata{i}.OutputData, '--r') 
    legend('Simulation','Experiment','Fontsize',14)
    title('Coil Currents 204660', 'Fontsize', 14)
    
%     figure
%     hold on
%     plot(t,x,'--r')
%     plot(t,y,'b')

    load('coils_greybox.mat')
    figure
    hold on
    plot(t,x(:,circ.iivx),'b')
    plot(coils.t, coils.iv,'--r')
    xlim([0 0.9])

end











