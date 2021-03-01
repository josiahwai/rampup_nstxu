function [equilibria, neq, eqraw, cc_filename] = write_eqdsk(shotnum,tefit,efit_source,tokamak,write_cc)
%  USAGE: [equilibria, neq, eqraw] = write_eqdsk(shotnum,tefit,efit_source,tokamak,write_cc)
%
%  PURPOSE: Retrieves eqdsk data from MDS and saves g,a files to dir
%
%  INPUTS:
%
%  	 shotnum: Shotnumber for equilibrium eg. 146970
%        tefit: time of efit, can be two element array specifying time range as in read_eq
%  	 efit_source: EFIT source eg. 'EFIT01'
%  	 tokamak: device name eg. 'd3d'
%
%  OUTPUTS:  eq = equilibrium data on toksys form, that can be used with
%              cc_efit_to_tok and build_tokamak_system
%            neq = number of equilibria returned
%            eqraw = 'raw' eq data on original form before conversion to toksys convention
%
%  RESTRICTIONS: 
%
%  METHOD:  
%
%
%  WRITTEN BY:  Matthew J. Lanctot 
%
%  MODIFICATION HISTORY:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  


[equilibria, neq, eqraw] = read_eq(shotnum,tefit,efit_source,tokamak);

cc_filename='';
if write_cc
  if neq==1
    [cc_filename,ier] = write_cc_file(equilibria);
  else
    for k=1:neq
      [cc_filename{k},ier] = write_cc_file(equilibria.gdata(k));
    end
  end
end

if isfield(equilibria,'gdata')
  % to simply gfile names, do not support sub-ms time resolution if times are unique
  uttmp = unique(equilibria.time-mod(equilibria.time,1e-3));
  if length(uttmp)==length(equilibria.time)
    equilibria.time = equilibria.time-mod(equilibria.time,1e-3);
    fix_gtime = 1;
  else
    fix_gtime = 0;
  end
 
  %write files
  for k=1:length(equilibria.time)
    if fix_gtime==1
      equilibria.gdata(k).time=equilibria.gdata(k).time-mod(equilibria.gdata(k).time,1e-3);
    end
    write_gfile(equilibria.gdata(k));
    % afile from read_eq doesn't have rzero
    %if isfield(equilibria,'adata'), write_afile(equilibria.adata(k));end
  end
else
  write_gfile(equilibria);
end
