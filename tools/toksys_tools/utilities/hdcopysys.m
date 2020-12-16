function []=hdcopysys(system,keep_ps,device)
 %
%  SYNTAX:  hdcopysys(system,keep_ps,device)
%
%  PURPOSE:  Create hardcopy plots of system block diagram.
%	(See hdcopy.m for ordinary figures.)
%
%  INPUT:
%   system  = text string specifying system name 
%   keep_ps = set to 1 to retain postscript file in local directory, 
%	      else deleted (optional, default=0)
%   device  = string defining print device (optional, default='hplj21')
%		- Set = 0 to prevent printing (e.g. when keep_ps=1).
%		- Device can be 'eps', in which case keep_ps must=1 to
%		  save the eps file.
%
%  OUTPUT:
%	hardcopy plot to printer 
%	hdcopytemp.ps = local postscript file (if keep_ps=1)
%
%  RESTRICTIONS: The system to be printed must be opened by SIMULINK within
%	your MATLAB session in order for this function to work.
%
%  METHOD:  Creates postscript file using matlab print command, then sends
%  to printer using !lp.  Modify printer device for your local printer.

%  WRITTEN BY:  Mike Walker 	ON 	10/2/96
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
   printf('must specify system name\n'); 
   return;
end;
if nargin < 2, keep_ps=0; end;
if nargin < 3
   device = 'hplj21';
   keep_eps=0;
else
   if strcmp(device,'eps') & keep_ps
      keep_eps = 1;
   else
      keep_eps = 0;
   end
end;

% remove existing postscript file
!rm hdcopytemp.ps

if keep_eps
   eval([' print -deps hdcopytemp.eps -s''' system ''''])
else
   [' print -dps hdcopytemp.ps -s''' system '''']
   eval([' print -dps hdcopytemp.ps -s''' system ''''])
   if(device~=0)
      eval(['!lp -d ' device ' hdcopytemp.ps'])
   end
end

if ~keep_ps
   !rm hdcopytemp.ps
end;
