function [shot, time] = gfile_shot_time(filename,typeflag)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  [shot, time] = gfile_shot_time(filename,typeflag)
%
%  PURPOSE:  Return shot and time values from gfile or afile name.
%
%  INPUT:
%	filename = name of gfile or afile
%	typeflag = 0 to return as integer values (time in ms),
%			1 to return as string values (default)
%
%  OUTPUT:
%	shot = shot number stripped from file name
%	time = time (in ms) stripped from file name 
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	5/30/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nargin < 2)
   typeflag=1;
end

idx2 = find(filename=='.');
idx1 = find(filename=='/'); 
if isempty(idx1)
   idx1 = 0;
else
   idx1 = idx1(end);
end

% parse out of filename

shot = filename(idx1+2:idx2(end)-1);
time = filename(idx2(end)+1:end);

% strip off leading zeros

shot = str2num(shot);
time = str2num(time);

% return type requested

if(typeflag)
   shot = int2str(shot);
   time = int2str(time);
end

