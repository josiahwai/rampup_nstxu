function []=hdcopy(figures,keep_ps,device)
 %
%  SYNTAX:  hdcopy(figures,keep_ps,device)
%
%  PURPOSE:  Create hardcopy plots of figures on screen.
%	(See hdcopysys.m for hardcopies of SIMULINK block diagrams.)
%
%  INPUT:
%   figures = list of figure numbers 
%	      (optional: if omitted, current figure is plotted)
%   keep_ps = set to 1 to retain postscript file (hdcopytemp.ps) in local directory, 
%		if a string, then plots are stored in file with that name,
%	      	else deleted (optional, default=0)
%   device  = string defining print device (optional, default='hplj21')
%		- note that device can be 'eps' or 'epsc', but set 
%		  keep_ps=1 to keep it.
%               - Set device=0 to prevent printing (e.g. when keep_ps=1).
%
%  OUTPUT:
%	hardcopy plots to printer 
%	hdcopytemp.ps = local postscript file (if keep_ps=1)
%		(or hdcopytemp.eps if device='eps' or 'epsc')
%
%  METHOD:  Creates postscript file using matlab print command, then sends
%  to printer using !lp.  
 
%  RESTRICTIONS:

%  WRITTEN BY:  Mike Walker 	ON 	?/93
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)hdcopy.m	1.4 10/01/13

if nargin < 1, figures=[]; end;
if nargin < 2, keep_ps=0; end;
if(isstr(keep_ps))
   psfile = keep_ps;
   keep_ps=1;
end
if(nargin < 3 & ~exist('psfile'))
   device = 'hplj21';
   keep_eps = 0;
else
   if strncmp(device,'eps',3) & keep_ps
      keep_eps = 1;
   else
      keep_eps = 0;
   end
end;

% remove existing postscript file
!rm hdcopytemp.ps

if (isempty(figures))	% print only the active figure window
   if keep_eps
      eval(['print -d' device ' hdcopytemp.eps']);
%      print -deps hdcopytemp.eps
   else
      print -dpsc hdcopytemp.ps
   end
else			% else print all specified figures
   figure(figures(1));
   if keep_eps
         eval(['print -d' device ' hdcopytemp' int2str(figures(1)) '.eps'])
   else
         print -dpsc hdcopytemp.ps
   end
   for k=2:length(figures)
      figure(figures(k));
      if keep_eps	% can't append eps files
         eval(['print -d' device ' hdcopytemp' int2str(figures(k)) '.eps'])
      else
         print -dpsc -append hdcopytemp.ps
      end
   end;
end;

if ~keep_eps
   if(device~=0 & ~strcmp(device,'eps') & ~strcmp(device(1:2),'ps'))
      eval(['!lp -d ' device ' hdcopytemp.ps'])
   end
end

if(keep_ps & exist('psfile'))
   eval(['!mv hdcopytemp.ps ' psfile])
elseif ~keep_ps
   !rm hdcopytemp.ps
end;
