   function [out_data] = linsim_tok_zctrl(tok_system,tau_ps,T_ps,Vvec, ...
			Vlims,Kpvec,Kdvec,taupd,iplot,zdisp0,tmax,zmax, ...
			Gpsweep,Gdsweep,idxacsw)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  USAGE:       [out_data] = linsim_tok_zctrl(tok_system,tau_ps,T_ps,Vvec, ...
%			Vlims,Kpvec,Kdvec,taupd,iplot,zdisp0,tmax,zmax, ...
%			Gpsweep,Gdsweep,idxacsw)
%
%  PURPOSE: Function to simulate linear vertical control with standard
%		tau_ps, T_ps power supply definition. Called by script
%		sim_iter_zctrl.m. Simulates response to zdisp0 initial 
%		condition
%
%  INPUTS:  
%	tok_system = structure gen. by build_* script containing system data 
%	tau_ps = power supply time constant [s]
%	T_ps = power supply delay time (in addition to one-pole tau_ps) [s]
%	Vvec = vector specifying active coils (entries: 0=not active, 
%		+/-1 = active positive(negative) voltage/current
%	Vlims = vector giving voltage saturation levels [V]
%	Kpvec   = vector giving proportional gain [V/m]
%	Kdvec   = vector giving derivative gain [V/m/s]
%	taupd  = time constant for 1-pole filter on PD operation [s]
%	iplot = (opt) flag to select plotting (1=plot, 0=don't(default))
%	zdisp0 = (opt) specific value to calculate displacement trajectory for
%		  and plot (if selected) 
%	tmax = (opt) max time for simulation (def=1.0 sec)
%	zmax = (opt) max z for plotting (def = 0.9 m)
%	Gpsweep = vector with Gpmin, Gpmax, ngp
%	Gdsweep = vector with Gdmin, Gdmax, ngd
%	idxacsw = idx scalar for selected active coil to plot contours
%		(this must be idx of coil whose peak I, V are negative,
%		 since code finds min of time history for contouring)
%	
%  OUTPUTS: 
%	out_data = structure that includes histories
%	+ plots (if selected)
%
%  RESTRICTIONS:
%
%  METHOD:  

%  WRITTEN BY:  Dave Humphreys  ON	3/30/08
%
%  MODIFICATION HISTORY:
%	DAH  3/30/08  Deriving elements from calc_tok_dzmax.m
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hardwired inputs:
   npade = 2;        %order for Pade approx to delay
   tau_filt = 1.0e-3;  
   nt = 200;	     %# of time pts
   nt0 = 500;        %# of time pts for zdisp0 calc (if selected)

% Defaults:
   if exist('tmax')~=1
     tmax = 1;	     %max time for simulation (sec)
   end
   if exist('zmax')~=1
     zmax = 0.9;       %max z for plotting dzmax trajectories
   end

% Prelims and Constants:
   mu0 = 0.4*pi;
   s = tf('s');    %define Laplace frequency variable

% Derived Values:
                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build system for lsim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   amat = tok_system.amat;
   nx = size(amat,1);
   bmat = tok_system.bmat;
   nu = size(bmat,2);
   cmat1 = eye(nx);
   Pxx =  tok_system.Pxx; 
   if isfield(tok_system,'scldzdis')   %default scldzdis or get from struct
     scldzdis = tok_system.scldzdis;
   else
     scldzdis = 1.0;
   end
   dzdis = scldzdis*tok_system.rzrig_data.dzdis*Pxx;    %z-response object
   %cmat = dzdis;      %make output = actual z
   cmat = cmat1;
   ny = size(cmat,1);
   dmat = zeros(ny,nu);
   t = linspace(0,tmax,nt)';

% System with no PS:
   sys = ss(amat,bmat,cmat,dmat); 
   [vecs,vals] = eigsort(amat);
   disp(['OL dominant pole (amat)= ',num2str(vals(1))])

% System with cascade PS + filter model:
   PS = 1/(tau_ps*s + 1);  set(PS,'OutputDelay',T_ps); 
   PS = pade(PS,npade);
   Filt = 1/(tau_filt*s + 1);   
   PSF = series(PS,Filt);
   sys_wps = series(PSF, sys); 
   [a,b,c,d] = ssdata(sys_wps);
   tmp = eigsort(a);
   disp(['OL dominant pole (eig of A from sys_wps)= ',num2str(tmp(1))])
   tmp = pole(sys_wps);
   idx = find(real(tmp) == max(real(tmp)));
   disp(['OL dominant pole (pole( sys_wps))= ',num2str(tmp(idx(1)))])
   c1 = dzdis*c;     %just z output
   d1 = d(size(c1,1),:);
   sys_wps = ss(a,b,c1,d1);
   nx1 = size(a,1);

% Execute dZmax simulation (saturated voltages): 
   x0=vecs(:,1); x0=zdisp0*x0/(dzdis*x0); x0 = [x0; zeros(nx1-nx,1)];
   U = ones(nt,1) * (Vvec.*Vlims)'; [y,t,x] = lsim(sys_wps,U,t,x0);
   zt_dispmax = (dzdis*x(:,1:nx)')';     
   zt_dispmaxy = (y')';
   idxac = find(Vvec.*Vlims ~= 0);
   Vcoils = U(:,idxac);		 %voltage in all active coils
   Icoils = x(:,idxac);    %currents in all active coils

% Build CL system:
  PD = [1/(1+taupd*s); s/(1+taupd*s)];     %PD operation on zerror
  if isempty(Kdvec), PD=1/(1+taupd*s); end    %set Kdvec to [] for P ctrl
  Kmat = [Kpvec(:) Kdvec(:)]*PD;  %gain mat maps PD output to Vcom
  
  sys_cl =  feedback(series(Kmat,sys_wps),tf(1),1,1,-1);
  itst = isproper(sys_cl);
  if itst==1
     disp('sys_cl is proper...')
  else
     disp('sys_cl is NOT proper...')
  end

% Execute CL simulation with initial condition z0 = zdisp0:  
  [acl,bcl,ccl,dcl] = ssdata(sys_cl);
  nxcl = size(acl,1);
  x0=vecs(:,1); x0=zdisp0*x0/(dzdis*x0); x0 = [x0; zeros(nxcl-nx,1)];
  [y,t,x] = initial(sys_cl,x0,t);
  zt_init = (dzdis*x(:,1:nx)')';     
  zt_inity = (y')';
  idxac = find(Vvec.*Vlims ~= 0);
  Icoilscl = x(:,idxac);    %currents in all active coils
  dt = t(2)-t(1);
  Ieff = sqrt( sum(Icoilscl.*Icoilscl.*dt)/t(end) );
  disp(['Max Zt_init = ',num2str(max(zt_init))])
  disp(['Min Icoilscl = ',num2str(min(Icoilscl))])
  disp(['Max Icoilscl = ',num2str(max(Icoilscl))])
  disp(['Ieff = ',num2str(Ieff)])

%Calculate voltage from Controller system Kmat:
  %Vcoilscl = U(:,idxac);    %WRONG: STILL MUST CALC VOLTAGE CL
  [a,b,c,d]=ssdata(Kmat);    %extract ABCD for controller system
  U = -y;
  xk0 = zeros(size(a,1),1);
  [yk,tk,xk] = lsim(Kmat,U,t,xk0);
  Vcoilscl = yk(:,idxac);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sweep gains to search for performance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isweep = 0;   %initialize flag to say no sweep as default...
if exist('Gpsweep')==1 & exist('Gdsweep')==1  & ~isempty(Gpsweep) ...
	& ~isempty(Gdsweep)     %then do sweep specified
  isweep = 1;   %set flag to say sweep executed
  gpmin=Gpsweep(1); gpmax=Gpsweep(2); ngp=Gpsweep(3);
  gdmin=Gdsweep(1); gdmax=Gdsweep(2); ngd=Gdsweep(3);
  gps = linspace(gpmin,gpmax,ngp)';
  gds = linspace(gdmin,gdmax,ngd)';
  gammasw = zeros(ngp,ngd);
  Icsw = zeros(ngp,ngd);
  Ieffsw = zeros(ngp,ngd);
  Vcsw = zeros(ngp,ngd);
  disp('Beginning gain sweep...')
  for ii=1:ngp
    gp = gps(ii);
    for jj=1:ngd
     gd = gds(jj);
     Kmat = [gp*Kpvec(:) gd*Kdvec(:)]*PD;  %gain mat maps PD output to Vcom
     sys_cl =  feedback(series(Kmat,sys_wps),tf(1),1,1,-1);
     tmp = pole(sys_cl); 
     idx = find(real(tmp) == max(real(tmp)));
     gammasw(ii,jj) = tmp(idx(1));    %real part of dominant pole
     [y,t,x] = initial(sys_cl,x0,t);   %uses x0 from above
     %zt_init = (dzdis*x(:,1:nx)')';     
     Icoilsw = x(:,idxacsw);    %current in idxacsw specified active coil
     Icsw(ii,jj) = min(Icoilsw);
     dt = t(2)-t(1);
     Ieff = sqrt( sum(Icoilsw.*Icoilsw.*dt)/t(end) );
     Ieffsw(ii,jj) = Ieff;
     %Get Active Coil voltage:
     [a,b,c,d]=ssdata(Kmat);    %extract ABCD for controller system
     U = -y;
     xk0 = zeros(size(a,1),1);
     [yk,tk,xk] = lsim(Kmat,U,t,xk0);
     Vcoilsw = yk(:,idxacsw);
     Vcsw(ii,jj) = min(Vcoilsw);

    end
  end
  disp('Gain sweep complete.')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot things if iplot==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iplot==1
    figure(5),clf,hold off
    plot(t,zt_dispmax,'r--','LineWidth',3);  % the maximum controllable zdisp
    hold on
    plot(t,zt_dispmaxy,'g:','LineWidth',2); %zt meas by cmat
    axis0 = axis;
    axis([axis0(1:2) 0 zmax])
    grid on
    title('Vertical Displacement w/ PS','FontSize',15)
    ylabel('Z [m]','FontSize',14)   
    xlabel('t [s]','FontSize',14)
    text(0.95*axis0(2),0.85*zmax,['Max Z_0\newlinefor control:\newline', ... 
         num2str(zdisp0)],'hori','ri','fonts',15)

    figure(6),clf,hold off
    nrows=length(idxac); ncols=2;
    for ii=1:nrows
     subplot(nrows,ncols,2*ii-1)
	plot(t,Vcoils(:,ii),'r','LineWidth',3);
        grid on
	ylabel('Vcoil [V/t]')
     if 2*ii-1==1
	title('Voltage and Current History in Active Coils','FontSize',15)
     end
     if 2*ii-1==2*nrows-1
	xlabel('t [s]','FontSize',14)
     end
     subplot(nrows,ncols,2*ii)
	plot(t,Icoils(:,ii),'r','LineWidth',3);
        grid on
	ylabel('Icoil [A]')
     if 2*ii==2*nrows
	xlabel('t [s]','FontSize',14)
     end
    end   %end for ii=1:nrows

 %Plot CL results:
    figure(7),clf,hold off
    plot(t,zt_init,'r--','LineWidth',3);  % the maximum controllable zdisp
    hold on
    plot(t,zt_inity,'g:','LineWidth',2); %zt meas by cmat
    axis0 = axis;
    axis([axis0(1:2) 0 zmax])
    grid on
    title('Vertical Displacement in CL','FontSize',15)
    ylabel('Z [m]','FontSize',14)   
    xlabel('t [s]','FontSize',14)
    text(0.95*axis0(2),0.85*zmax,['Max Z_0\newlinefor control:\newline', ... 
         num2str(zdisp0)],'hori','ri','fonts',15)

    figure(8),clf,hold off
    nrows=length(idxac); ncols=2;
    for ii=1:nrows
     subplot(nrows,ncols,2*ii-1)
	plot(t,Vcoilscl(:,ii),'r','LineWidth',3);
        grid on
	ylabel('Vcoil [V/t]')
     if 2*ii-1==1
	title('Voltage and Current History in Active Coils','FontSize',15)
     end
     if 2*ii-1==2*nrows-1
	xlabel('t [s]','FontSize',14)
     end
     subplot(nrows,ncols,2*ii)
	plot(t,Icoilscl(:,ii),'r','LineWidth',3);
        grid on
	ylabel('Icoil [A]')
     if 2*ii==2*nrows
	xlabel('t [s]','FontSize',14)
     end
    end   %end for ii=1:nrows

 %Plot gain sweep results (if done):
 if isweep==1
   figure(9),clf,hold off
    [c,h]=contour(gds,gps,real(gammasw),100);
    clabel(c,h,'labelspacing',400);
    hold on
    [c,h]=contour(gds,gps,real(gammasw),[0 0],'r--');
    set(h(1),'Linewidth',3)
    xlabel('Gd','FontSize',14)
    ylabel('Gp','FontSize',14)
    title('Contours of CL Gamma + Stable FB Roots','FontSize',15)

   figure(10),clf,hold off
    [c,h]=contour(gds,gps,Icsw,100);
    clabel(c,h,'labelspacing',400);
    hold on
    [c,h]=contour(gds,gps,Icsw,[0 0],'r--');
    set(h(1),'Linewidth',3)
    xlabel('Gd','FontSize',14)
    ylabel('Gp','FontSize',14)
    title('Contours of Peak AC Current','FontSize',15)

   figure(11),clf,hold off
    [c,h]=contour(gds,gps,Ieffsw,100);
    clabel(c,h,'labelspacing',400);
    set(h(1),'Linewidth',3)
    xlabel('Gd','FontSize',14)
    ylabel('Gp','FontSize',14)
    title('Contours of Ieff','FontSize',15)

 %  figure(11),clf,hold off
 %   [c,h]=contour(gds,gps,Vcsw,100);
 %   clabel(c,h,'labelspacing',400);
 %   hold on
 %  [c,h]=contour(gds,gps,Vcsw,[0 0],'r--');
 %   set(h(1),'Linewidth',3)
 %   xlabel('Gd','FontSize',14)
 %   ylabel('Gp','FontSize',14)
 %   title('Contours of Peak AC Voltage','FontSize',15)

 end  %end if isweep

end   %end if iplot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct output structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define output structure:
   out_data.nt = nt;
   out_data.nt0 = nt0;
   out_data.tmax = tmax;
   out_data.Vcoils = Vcoils;
   out_data.Icoils = Icoils;
   out_data.Vcoilscl = Vcoilscl;
   out_data.Icoilscl = Icoilscl;
   out_data.zt_dispmax = zt_dispmax;
   if isweep
     out_data.Vcsw = Vcsw;
     out_data.Icsw = Icsw;
     out_data.gammasw = gammasw;
   end



