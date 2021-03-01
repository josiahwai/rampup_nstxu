function [taupen,ymeas,t,dbgout] = calc_taupen(tok_data_struct,ipf,Rmeas,Zmeas,options);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  USAGE:  >> taupen = calc_taupen(tok_data_struct,ipf,Rmeas,Zmeas,options);
%
%  PURPOSE: Function to calculate pentration time for generic tokamak
%		by driving PF coils with step CURRENT
%		to produce VV currents and calculate field at meas point
%
%  INPUTS:
%	tok_data_struct = standard TokSys vacuum objects structure
%	ipf = index of PF coil to drive and calculate measurement time
%	Rmeas,Zmeas = R,Z values at which to measure penetration [m]
%	options.iplot = (optional) flag to enable plotting of time history for  
%		estimation of penetration time for PF ipf at (Rmeas,Zmeas).
%		Default = 1
%	options.Gpc = (optional) Green function object desired for penetration
%		time estimation. Must be a PF-coil-to-grid object. 
%		Default = mpc
%	options.Gpv = (optional) Green function object desired for penetration
%		time estimation. Must be a VV-element-to-grid object. 
%		Default = mpv
%	options.taupf = (optional) time constant for PF current rise [s]
%		Default = 1e-3. 
%	options.tmax = (optional) max time for evolution [s]. Default = 50e-3
%	options.nt = (optional) # of time points in evolution. 
%		Default = calculated on basis of taupf and tmax: 
%		nt = 10*(tmax/taupf)
%	options.resv_fix = (optional) modified resv vector to use for
%		calculation.  Default = resv from tok_data_struct.
%	
%  OUTPUTS:
%	taupen = estimated penetration time for ipf and Rmeas,Zmeas
%	outputs = structure with other useful quantities
%	(Plots showing fit, if iplot=1)
%
%  RESTRICTIONS:
%
%  METHOD:  
%	Run lsim for step CURRENT applied to PF coil ipf, generate field
%	history at (Rmeas,Zmeas) measurement point with field type specified
%	by options.Gpc, Gpv, fit result to exponential rise to extract approx
%	single time constant taupen. Fit only to early time corre. to one
%	e-folding (if were pure single-tau exponential) to extract effective
%	single time constant for dominant rise time... Remove all other PF coils from problem,
%	so they're not doing any shielding. 
%
%  USE EXAMPLE:
%	>> tok_data_struct = load_tok_objects('kstar','2010','6565');
%	>> resv_fix = tok_data_struct.resv;
%	>> resv_fix(71:end) = resv_fix(71:end)*1000;  %kill all but VV walls
%	>> options.resv_fix = resv_fix;
%	>> taupfs = zeros(14,1);   %just SC PF's
%	>> for ii=1:14
%	>>   [taupen,ymeas,t] = calc_taupen(tok_data_struct,ii,1.8,0,options);
%	>>   taupfs(ii)=taupen;
%	>> end
%

%  WRITTEN BY:  Dave Humphreys  ON	8/2/14
%
%  MODIFICATION HISTORY:
%
%    8/2/14  DAH  Derived from calc_penetration.m, script to calculate
%		penetration times somewhat different way
%
%  NOTES:
%     Calc'd VV penetration time for EAST: ~11 msec for vertical field from
%	PF coils 6,13 to Bprobe #1 (midplane inboard). Faster than present
%	reported penetration time of 18 ms from Du  (8/8/05)
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Extract options and tok_data_struct objects:
  if exist('options')==1, struct_to_ws(options); end
  struct_to_ws(tok_data_struct);

% Default inputs:
   if exist('iplot')~=1, iplot=1; end
   if exist('Gpc')~=1, Gpc=mpc; end
   if exist('Gpv')~=1, Gpv=mpv; end
   if exist('taupf')~=1, taupf=1e-3; end
   if exist('tmax')~=1, tmax = 100e-3; end
   if exist('nt')~=1, nt = 10*(tmax/taupf); end
   if exist('resv_fix')==1, resv=resv_fix; end

% Prelims and Constants:
   mu0 = 0.4*pi;
   taupen = 0;    %initialize for testing...

% Derived Values:
                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 rvv = diag(resv);

 amat = -inv(mvv)*rvv;
 nx = size(amat,1);  
 vals = eigsort(amat);
 tauwall = -1/vals(1);
 disp(['VV dominant time const tauwall = ',num2str(tauwall),' sec'])

 bmat = -inv(mvv)*mcv(ipf,:)';
 nu = size(bmat,2);
 
%Interpolate output location:
  tmp = reshape(Gpc(:,ipf),nz,nr);
  Gmpf = interp2(rg,zg,tmp,Rmeas,Zmeas);   %scalar coupling PF to meas
  Gmvv = zeros(1,nv);
  for ii=1:nv 
    tmp = reshape(Gpv(:,ii),nz,nr);
    Gmvv(1,ii) = interp2(rg,zg,tmp,Rmeas,Zmeas);   %scalar coupling VV to meas 
  end
  
 cmat = eye(nx,nx);   %output is all VV currents *only*
 dmat = zeros(nx,nu);

 sys = ss(amat,bmat,cmat,dmat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up for simulation with lsim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 t = linspace(0,tmax,nt)';
 pfcurt = (1 - exp(-t/taupf));    % peak current = 1 A
 dpfcurdt = deriv(pfcurt,t);       %dIpf/dt

 umat = dpfcurdt;   

 [y,t,x] = lsim(sys,umat,t);

 totcurvv = sum(y')';
 ymeas = (Gmvv*y')' + (Gmpf*pfcurt')';
 ymeaspk = Gmpf*1;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit for approx taupen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ymeasnorm = (ymeaspk - ymeas)/ymeaspk;    % approx 1*exp(-t/taupen)
   tautst = -1./deriv(log(ymeasnorm),t); 
   %BAD taupen: finds slow long-lasting mode only...
   %taupen = mean(tautst(fix(nt/5):nt));   %only avg last 80% of sim time
   idx = find_near(ymeasnorm,exp(-1));  %find when signal is one e-fold
   taupen = mean(tautst(1:idx));  %avg only period from 
  disp(['Approx. penetration time taupen = ',num2str(taupen),' sec'])
  ymeastst = ymeaspk*(1 - exp(-t/taupen));
  dbgout = struct('tautst',tautst,'ymeasnorm',ymeasnorm,'ymeaspk',ymeaspk, ...
               'ymeas',ymeas);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define CC geometry for plotting:
  if exist('ecdata')==1   %if there are Ecoils
      nec = size(ecnames,1);  % # of Ecoils ending up in cc set
      if ipf<=nec    %if user selected Ecoils
        Rpf=mean(ecdata(2,:)); Zpf=mean(ecdata(1,:)); 
      else    %if user selected an F-coil
        Rpf=fcdata(2,ipf-nec); Zpf=fcdata(1,ipf-nec); 
      end
  else    %if there are no Ecoils in tok_data_struct, ipf is an fcdata idx
      Rpf=fcdata(2,ipf); Zpf=fcdata(1,ipf); 
  end
figure(1),clf,hold off
plot_tok_geo(tok_data_struct)
hold on
plot(Rmeas,Zmeas,'gx','MarkerSize',14,'LineWidth',3)
plot(Rpf,Zpf,'co','MarkerSize',14,'LineWidth',3)

figure(2),clf,hold off
subplot(3,1,1)
  plot(t,totcurvv)
  grid on
  ylabel('Iv [A]')
title('Penetration Data')
subplot(3,1,2)
  plot(t,ymeas)
  grid on
  ylabel('Meas Signal')
subplot(3,1,3)
  plot(t,ymeasnorm)
  grid on
  ylabel('Norm Meas Signal')
xlabel('t [sec]')

 figure(3),clf,hold off
  plot(t,tautst,'r')
  hold on
  plot(t,taupen*ones(length(t),1),'g--')
  grid on
  xlabel('t [sec]')
  title(['Approx time constant = ',num2str(taupen)],'FontSize',15)

figure(4),clf,hold off
   plot(t,ymeas,'r')
   hold on
   plot(t,ymeastst,'g--')
   grid on
   xlabel('t [sec]')
   title(['Test fit of time constant = ',num2str(taupen)],'FontSize',15)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   disp('All done.')
               
