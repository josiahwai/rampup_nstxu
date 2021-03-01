 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  get_efit_data
%
%  PURPOSE:  Script to extract many slices of EFIT data from set of a0-eqdsk 
%  files assumed present in default directory and plot selected data.
%
%  INPUT:
%	eqdir   = directory for EFIT data files (default='./');
%		  eqdir string must end in "/" (eg '/users/humphrys/matlab/')
%       shot	= shot number
%    	iplot 	= 1 if want to plot many frames on one plot (default)
%	      	= 2 if want to plot Z,R,Ip(t) on different plots
%		= 0 to skip plotting
%    EFIT time slices are specified using EITHER of the following methods:
%       tmiefit	= 1st efit time slice (msec)
%       dt 	= delta time (msec)
%       nslices	= # of efit slices
%    OR:
%       tefit 	= vector of efit times (msec) to read from a0 files.  
%    If tefit is specified, input variables tmiefit, dt, nslices are ignored.
%       ipltref = present figure if want to start plotting after present fig
%		  (default=0) e.g. if ipltref=2, EFIT plots start at fig 3 
%		  Currently used only for iplot=1
%
%  OUTPUT:
%       Extracts data from a0 files and plots things.
 
%  RESTRICTIONS:
%     A0-files specified by inputs must exist in default directory.
%	In present version, efit times must all have 4 digits (ie t must be
%	> 999, and < 9999)...
%
%  METHOD:  
%     Uses read_afile to read A0-eqdsk file for series of EFIT slices.
%
%  WRITTEN BY:  Dave Humphreys ON 	6/25/95 Adapted from compare_efitz.m
%
%  MODIFICATION HISTORY:
%      Finally completed to work properly. DAH 2/96 
%      Modified to accept inputs from Matlab environment instead of definitions
%         at beginning of script.  MLW 5/15/97
%      Plotting logic modified.  DAH 5/16/97
%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial definitions:
fprintf('shot = %d\n',shot)
if ~exist('iplot'), iplot = 1; end
fprintf('iplot = %d\n',iplot)
if(exist('tefit'))
   fprintf('time slices specified by tefit \n')
else
   fprintf('tmiefit = %d\n',tmiefit)
   fprintf('dt = %d\n',dt)
   fprintf('nslices = %d\n',nslices)
   tefit = linspace(tmiefit,tmiefit+dt*(nslices-1),nslices)';  %time vec in msec
end

% Default for EFIT data directory:
  if ~exist('eqdir')
    eqdir = './';    %default to default directory
  end

%Derived Quantities:
  nt = length(tefit);
  shotno=int2str(shot);

% Data vectors:
  t=tefit*1e-3;    %time vector in secs
  t1=t(1); t2=max(t);
  zmagt=t; rmagt=t; zcurt=t; rcurt=t; routt=t; ipt=t; kappat=t;
  aoutt=t; zsep1t=t; rsep1t=t; zsep2t=t; rsep2t=t; zoutt=t;
  gaptopt=t; gapint=t; gapoutt=t;
  betapt=t; lit=t; qet=t;
  cmpr2t = [];

for ii=1:nt
  if shot<10000
     stmp='00'
  elseif shot<100000
     stmp='0'
  else 
     stmp=''; 
  end
  if tefit(ii) < 1e4
      padstr = '.0';
      if tefit(ii) < 1e3
         padstr = '.00';
         if tefit(ii) < 1e2
            padstr = '.000';
            if tefit(ii) < 1e1
               padstr = '.0000';
            end
         end
      end
   end
  filename = [eqdir,'a',stmp,shotno,padstr,int2str(tefit(ii))];
  read_afile;
  zmagt(ii) = zmagx*0.01;   %cm->m
  rmagt(ii) = rmagx*0.01;
  zcurt(ii) = zcurrt*0.01;
  rcurt(ii) = rcurrt*0.01;
  zoutt(ii) = zout*0.01;
  routt(ii) = rout*0.01;
  aoutt(ii) = aout*0.01;
  kappat(ii) = eout;
  zsep1t(ii) = zseps1*0.01;
  rsep1t(ii) = rseps1*0.01;
  zsep2t(ii) = zseps2*0.01;
  rsep2t(ii) = rseps2*0.01;
  gaptopt(ii) = otop*0.01;
  gapoutt(ii) = oright*0.01;
  gapint(ii) = oleft*0.01;
  ipt(ii) = pasmat*1e-6;    %A->MA
  betapt(ii) = betap;
  lit(ii) = ali;
  qet(ii) = qpsib;
  cmpr2t = [cmpr2t; cmpr2'];
  cpasmat(ii) = cpasma; 
end

%Derived and others:
  zlim = -1.3675;   %floor limiter vertical position [m]
  zrelt = (zcurt-zlim);  %vert posit rel to floor limiter
  dzdipt = (zrelt(2:nt)-zrelt(1:nt-1))./(ipt(2:nt)-ipt(1:nt-1));

%Plot things:
if iplot==1
if ~exist('ipltref'), ipltref=0; end
figure(1+ipltref),clf,hold off
nplots=4;
subplot(nplots,1,1)
  plot(t,rcurt)
  hold on
  plot(t,rmagt,'m--')
  plot(t,routt,'c-.')
  axis0=axis; axis([t1 t2 axis0(3:4)])
  grid
  ylabel('Rcur,mag,out [m]')
title([shotno,': EFIT Ovrvu 1 '])
subplot(nplots,1,2)
  plot(t,zcurt)
  hold on
  plot(t,zmagt,'m--')
  plot(t,zoutt,'c-.')
  axis0=axis; axis([t1 t2 axis0(3:4)])
  grid
  ylabel('Zcur,mag,out [m]')
subplot(nplots,1,3)
  plot(t,aoutt)
  axis0=axis; axis([t1 t2 axis0(3:4)])
  grid
  ylabel('a [m]')
subplot(nplots,1,4)
  plot(t,kappat)
  axis0=axis; axis([t1 t2 axis0(3:4)])
  grid
  ylabel('Kappa')
xlabel('t [sec]')

figure(2+ipltref),clf,hold off
nplots=4;
subplot(nplots,1,1)
  plot(t,rsep1t)
  axis0=axis; axis([t1 t2 axis0(3:4)])
  grid
  ylabel('Rx1 [m]')
title([shotno,': EFIT Ovrvu 2'])
subplot(nplots,1,2)
  plot(t,zsep1t)
  axis0=axis; axis([t1 t2 axis0(3:4)])
  grid
  ylabel('Zx1 [m]')
subplot(nplots,1,3)
  plot(t,rsep2t)
  axis0=axis; axis([t1 t2 axis0(3:4)])
  grid
  ylabel('Rx2 [m]')
subplot(nplots,1,4)
  plot(t,zsep2t)
  axis0=axis; axis([t1 t2 axis0(3:4)])
  grid
  ylabel('Zx2 [m]')
xlabel('t [sec]')

figure(3+ipltref),clf,hold off
nplots=4;
subplot(nplots,1,1)
  plot(t,zsep1t)
  axis0=axis; axis([t1 t2 axis0(3:4)])
  grid
  ylabel('Zx1 [m]')
title([shotno,': EFIT Ovrvu 3'])
subplot(nplots,1,2)
  plot(t,gapoutt)
  axis0=axis; axis([t1 t2 axis0(3:4)])
  grid
  ylabel('Gout [m]')
subplot(nplots,1,3)
  plot(t,gaptopt)
  axis0=axis; axis([t1 t2 axis0(3:4)])
  grid
  ylabel('Gtop [m]')
subplot(nplots,1,4)
  plot(t,gapint)
  axis0=axis; axis([t1 t2 axis0(3:4)])
  grid
  ylabel('Gin [m]')
xlabel('t [sec]')

figure(4+ipltref),clf,hold off
nplots=4;
subplot(nplots,1,1)
  plot(t,ipt)
  axis0=axis; axis([t1 t2 axis0(3:4)])
  grid
  ylabel('Ip [MA]')
title([shotno,': EFIT Ovrvu 4'])
subplot(nplots,1,2)
  plot(t,lit)
  axis0=axis; axis([t1 t2 axis0(3:4)])
  grid
  ylabel('li')
subplot(nplots,1,3)
  plot(t,betapt)
  axis0=axis; axis([t1 t2 axis0(3:4)])
  grid
  ylabel('Betap')
subplot(nplots,1,4)
  plot(t,qet)
  axis0=axis; axis([t1 t2 axis0(3:4)])
  grid
  ylabel('q_{edge}')
xlabel('t [sec]')

end  %end iplot==1


if iplot==2
  figure(1),clf,hold off
  plot(t,zrelt)
  hold on
  plot(t,zmagt-zlim,'m--')   %vert posit of mag axis rel to floor
  xlabel('t [sec]')
  ylabel('Zrel [m]')
  title([shotno,': Z Relative to Floor'])
  grid

figure(2),clf,hold off
  plot(t,rcurt)
  hold on
  plot(t,rmagt,'m--')
  xlabel('t [sec]')
  ylabel('Rcur [m]')
  title([shotno,': Rcur '])
  grid

figure(3),clf,hold off
  plot(t,ipt)
  xlabel('t [sec]')
  ylabel('Ip_meas [MA]') 
  title([shotno,': Plasma Current']) 
  grid

%figure(3),clf,hold off
%  plot(t(1:nt-1),dzdipt)
%  hold on
%  plot(t(1:nt-1),smoothdt(dzdipt,t,4e-3),'m--')
%  xlabel('t [sec]')
%  ylabel('dz/dIp [m/MA]')
%  title([shotno,': dZ/dIp During VDE'])  
%  grid

end  %end iplot==2

