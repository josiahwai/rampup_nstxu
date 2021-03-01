 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  USAGE:  >>efit_movie
%
%  PURPOSE: Script to display movie of EFIT flux patterns from set of EFITs.
%
%  INPUTS:
%	eqdir = string giving directory containing EFIT g-files 
%		(final character is '/', e.g. '/u/humphrys/efit/')
%	shot = integer shot number
%	t1ms = initial EFIT slice time (msec)
%	t2ms = final EFIT slice time (msec)
%	dtms = EFIT slices time interval (msec)
%
%  OUTPUTS:
%	EFIT movie displayed.
%
%  RESTRICTIONS:
%	EFIT g0-files identified by shot, t1ms, t2ms, dtms must exist in eqdir.
%	Input variables must be defined; no defaults provided EXCEPT for
%	eqdir, defaulted to './' (so uses default directory).
%
%  METHOD:  

%  WRITTEN BY:  Dave Humphreys  ON	1/1/97
%
%  MODIFICATION HISTORY:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default inputs:
  if ~exist('eqdir')
     eqdir = './';
  end

% Prelims:
   mu0 = 0.4*pi;
   clf,hold off  %clear current figure

% Derived Values:
  npts = (t2ms-t1ms)/dtms + 1;

% Acquire data and make movie:
  for ii=1:npts
    disp(['ii=',int2str(ii),' of ',int2str(npts)])
    itime = t1ms + (ii-1)*dtms;
    filename = [eqdir,'g0',int2str(shot),'.0',int2str(itime)];
    hold off
    plot_efit
    if ii==1, M=moviein(npts); end  %initialize M first time
    M(:,ii) = getframe;
  end

% Play movie:
  disp('Playing movie once...')
  movie(M)
  disp('Movie complete. To play again N times, do >>movie(M,N)')


            
