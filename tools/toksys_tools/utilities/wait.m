function wait(message)
 %
%  SYNTAX:  wait(message)
%
%  PURPOSE:  Procedure to implement a pause with a message to let user know that
% 	carriage return is required in order to continue.  
% 	(This eliminates confusion over whether long wait for MATLAB is due to 
% 	a pause or to a very long computation.)
%
%   If the environemt variable 'SKIPWAIT' exists on the workspace, then the pause is
%   bypassed. This is useful for automated builds.
%
%  INPUT:
% 	message = optional string to print out at pause 
%
%  OUTPUT:
%	printed message to terminal

%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	12/7/94
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('\n');
if nargin<1
   fprintf(' Enter return to continue\n');
else
%   fprintf([message '\n']);
   disp(message)
   fprintf(' Enter return to continue\n');
end;

skipwait = getenv('SKIPWAIT');
if isempty(skipwait)
    pause
end
fprintf('continuing...\n')
