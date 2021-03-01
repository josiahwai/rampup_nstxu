function status = mdsdisconnect()
 %
%  SYNTAX:  status = mdsdisconnect
%
%  PURPOSE:  function to close communication with remote mds server
%
%  INPUT: none
%
%  OUTPUT:
%	status
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Basil P. DUVAL	ON	Oct 1998
%  MODIFIED BY: Mike Walker     ON      6/5/01          for use at DIII-D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @(#)mdsdisconnect.m	1.2 03/03/10

%status = mdsremote(4);
status = mdsipmex(4);
