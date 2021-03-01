
   function run_efit(shot,tims,dtms,nslices,snap,snapidx,efitver,efitdir)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  run_efit(shot,tims,dtms,nslices,snap,snapidx,efitver,efitdir)
%
%  PURPOSE: 
%	Run EFIT from inside Matlab: uses basic input format just as in
%	standard invokation from command line (so will perform many slices
%	as usual with starting time, delta time, nslices).
%
%  INPUTS:
%	shot = shot number
%	tims = starting time (in ms) at which to calculate EFIT
%	dtms = delta time (in ms) between each slice
%	nslices = # of slices to calculate
%	snap = (optional) string specifying:
%		if snapidx=2: data file to read input data from
%		if snapidx=3: use snap file in default directory (or in /link)
%		if snapidx=7: snap file extension (eg 'jta_t') 
%	snapidx = (optional) EFIT input mode (2=input file, 3=snap file)
%		Default=3 (used if nargin<6). Note if snapidx=2, the shot,
%		tims,dtms,nslices inputs are ignored, and snap specifies the 
%		input file to use.
%	efitver = (optional) string specifying EFIT version (eg 'efitd65yd' or 'efitd6565d')
%	efitdir = (optional) string specifying directory in which EFIT version
%		is to be found (default = /link/efit/ or /d/linux/efit/)
%
%  OUTPUTS:
%	None
%
%  RESTRICTIONS:
%
%  METHOD:  
%	Uses Matlab unix() command to execute command string built in 
%	this function for desired inputs. 

%  WRITTEN BY:  Dave Humphreys 	ON	11/1/00
%
%  VERSION @(#)run_efit.m	1.1 08/27/10
%
%  MODIFICATION HISTORY: ASW 4/30/04 runs on either thor or uscws9 transparently
%                        ASW 8/27/10 changed from int2str to num2str 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define defaults:
   if nargin<5
      snap='';       %no default: use snap file in default dir...
      snapidx = 3;   %select snap option "3" in efit call (use file in def dir)
   end   
   if nargin==5
      snapidx = 7;   %select snap option "7" in efit call (use snap file ext)
   end   
   if nargin<7, efitver='efitd65yd'; end   %default version (recommend efitd6565d for 65x65 grid)
   if nargin<8
      efitdir='/link/efit/'; % default efit directory for hp
      % Check to see if we are on Thor (a linux machine)
      thordisk = getenv('THORDISK');
      if length(thordisk)==5
        if thordisk=='/home'
	   efitdir='/d/linux/efit/';
	end
      end
   end

% Prelims:
   mu0 = 0.4*pi;

% Derived Values:

% Build command string:
  if snapidx==3   %use snap file input, use snap file "efit_snap.dat"
     efitcom = ['echo ',' "','3\n', ...
        num2str(shot),',',num2str(tims),',',num2str(dtms),',', ...
          num2str(nslices),'" ','| ',efitdir,efitver];
  end

  if snapidx==7  %use snap file input, use snap file with extension snap
     efitcom = ['echo ',' "','7\n', ...
	snap,'\n', ...
        num2str(shot),',',num2str(tims),',',num2str(dtms),',', ...
          num2str(nslices),'" ','| ',efitdir,efitver];
  end

  if snapidx==2  %use input file input, input file specified by snap variable
     efitcom = ['echo ',' "','2\n', ...
          '1\n',snap,'\n','" ','| ',efitdir,efitver];
  end

% Execute EFIT:
  unix(efitcom);


