function status = mdsclose(para)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  status = mdsclose(para)
%
%  PURPOSE:  function to mdsclose on remote server
%
%  INPUT:
%	para = ??	(optional)
% 		(if para ~= 0, perform a disconnect of remote server)
%
%  OUTPUT:
%	status
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Basil P. DUVAL	ON 	Oct 1998
%  MODIFIED BY: Mike Walker	ON	6/5/01		for use at DIII-D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%status = mdsremote(3);
%status = mdsipmex(3);
%status  = mdsipmex('MDSLIB->MDS$CLOSE()');
if(nargin < 1);para=0;end
if(para)
   status = mdsdisconnect;
else
   status  = mdsipmex('TreeClose()');
end
