  function  eqmod= eq_mod(eq)
%
%  SYNTAX:
%         eqmod= eq_mod(eq);
%
%  PURPOSE:  Puts time index as last index of each array. 
%
%  INPUT: <default>
%    eq         = structure from mds_eq
%
%  OUTPUT:
%    eqmod      = structure containing EFIT eqdsk variables with time last entry
%                 Example: eq.BCENTR(255,1) => eqmod.BCENTER(1,255)
%
%    see also mds_eq
%
%  WRITTEN BY:  Jim Leuer    ON      6/11/03
% ==========================================================================
% 
% Defaults
  if nargin==0
     disp('% eq_mod needs at least a "eq" argument')
     help eq_mod
     return
  end

% -----------------------------------------------
  eqmod= eq;

 if isfield(eq,'GTIME') & isfield(eq,'gnames')

  nt= length(eq.GTIME);
    
% loop over all variables
  for ii= 1:length(eq.gnames);
%    ii=0; ii=ii+1;
     nam= deblank(eq.gnames(ii,:));
     
     str= ['ndim= ndims(eq.' nam ');'];
     eval(str);
     if ndim<=2
        str= ['siz= size(eq.' nam ');'];
        eval(str);
	if siz(1)==nt
           str= ['eqmod.' nam '= eq.' nam ''' ;'];
	   eval(str);
	end
      end
   end

 end
 
% AFILE
    
if isfield(eq,'ATIME') & isfield(eq,'anames')

  nt= length(eq.ATIME);
    
% loop over all variables
  for ii= 1:length(eq.anames);
%    ii=0; ii=ii+1;
     nam= deblank(eq.anames(ii,:));
     
     str= ['ndim= ndims(eq.' nam ');'];
     eval(str);
     if ndim<=2
        str= ['siz= size(eq.' nam ');'];
        eval(str);
	if siz(1)==nt
           str= ['eqmod.' nam '= eq.' nam ''' ;'];
	   eval(str);
	end
      end
   end

 end

% MFILE:


if (isfield(eq,'TIME') | isfield(eq,'MTIME')) & isfield(eq,'mnames')

   if isfield(eq,'MTIME') 
     nt= length(eq.MTIME);
   else
     nt= length(eq.TIME);
   end
    
% loop over all variables
  for ii= 1:length(eq.mnames);
%    ii=0; ii=ii+1;
     nam= deblank(eq.mnames(ii,:));
     
     str= ['ndim= ndims(eq.' nam ');'];
     eval(str);
     if ndim<=2
        str= ['siz= size(eq.' nam ');'];
        eval(str);
	if siz(1)==nt
           str= ['eqmod.' nam '= eq.' nam ''' ;'];
	   eval(str);
	end
      end
   end

 end

 return
  
% Testing
%  eq= mds_eq(114504, 'EFIT01', 0);
%  eqmod= eq_mod(eq)

