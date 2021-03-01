function [cc_filename,ier] = write_cc_file(eq);
%  USAGE:  [cc_filename,ier] = write_cc_file(equilibria);
%
%  PURPOSE: return absolute path to newly generated cc_file
%
%  INPUTS:  filename:name of geqdsk file
%
%  OUTPUTS: cc_filename: string containing name of cc_file
%
%  RESTRICTIONS: 
%
%  METHOD:   
%
%     cc_file=  coil current file name [cc=load(cc_file)] or coil current vector
%               use when coil currents not available from efit (or to override)
%
%               ohmic coil current is first, then poloidal shaping currents
%		(units = MA-turns)
%
%  WRITTEN BY:  Matthew J. Lanctot
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %inits
  cc_filename='';
  ier = 0;

  cc_filename = ['cc' num2str(eq.shotnum) '.' sprintf('%05g',eq.time*1e3)];
  
  if isfield(eq,'brsp')
    brsp = eq.brsp;	%have to fit PF coil, even with small weight
  else 
    ier=1;
  end
  
  if isfield(eq,'ecurrt')
    ecurrt = eq.ecurrt; %have to fit OH coil, even with small weight
  else 
    ier=2;
  end 
 
  switch ier
    case 1
      warning('gfile does not contain brsp data. were PF coils fit?') 
      return   
    case 2
      warning('gfile does not contain ecurrt data. was OH coil fit?')
      return
    otherwise %write cc file
      nf = length(unique( eq.fcid));
      data = [ecurrt brsp(1:nf)']';
      dlmwrite(cc_filename,data,'\t',0,0,14);
  end
  
  
