   function [hr,hz,fl] = fine(ri,ro,zl,zu,cur,rp,zp,idebug)

% fine.m: Axisymmetric Magnetics from uniform current rectangle to point  
% SYNTAX: [hr,hz,fl] =fine(ri,ro,zl,zu,cur,rp,zp,<idebug>)
%
%  Magnetics function to determine the fields from an 
%  axisymmetric rectangular uniform current region 
%  (ri,ro,zl,zu,cur) to a point (rp,zp).
%
%  SYNTAX: [hr,hz,fl] =fine(ri,ro,zl,zu,cur,rp,zp,<idebug>)
%
%  INPUTS:
%           ri,ro,zl,zu = inner,outer,lower,upper dimensions of region [m]
%           cur         = current [Amp]
%           rp,zp       = receiving point [m]
%           idebug      = debug parameter, 0=no, >0 = debug print out
%
%           NOTE: FINE accepts vectors as long as all vectors are same size
%                 (ie ri(10),ro(100,zl(10),zu(10),cur(10),rp(10),zp(10)
%
%  OUTPUTS:
%           hr, hz      = fields at points rp,zp [Amp/m]
%           fl          = flux at point [Weber/mu_o]
%
%  NOTE:           multiply results by mu_o= 4*pi*1e-7 to get Tesla and Weber
%
%  CAUTION: Accuracy of routine becomes problematic for highly elongated
%           rectangles at small radius - SOLUTION: break up into group of
%           smaller rectangles
%
%  NOTE:    same routine can be executed as a MEX file (for faster execution)
%           SYNTAX: [hr,hz,fl] =fine1(ri,ro,zl,zu,cur,rp,zp,<idebug>)
%            

%  Uses MEX File fine1c.mexlgx
%  Gateway function: fine1c.c
%  C function: fine.c
%
% GA Jim Leuer 4-95, revised by Brian Sammuli 4/30/08

% ---------------------------------------------------------------------
% CHECK FOR PROPER NUMBER OF ARGUMENTS
%
      if nargin <= 6
        nargin
        disp('%ERR fine: Number of input arguments must be at least 7')
        disp('%SYNTAX: [hr,hz,fl] =fine(ri,ro,zl,zu,cur,rp,zp,<idebug>)')
        return
      end

      if nargin == 7 idebug= 0; end

      if nargout <= 2
        nargout
        disp('%ERR fine: Number of output arguments must be 3')
        disp('%SYNTAX: [hr,hz,fl] =fine(ri,ro,zl,zu,cur,rp,zp,<idebug>)')
        return
      end
%
% ---------------------------------------------------------------------

      [hr,hz,fl] = fine1(ri,ro,zl,zu,cur,rp,zp,idebug);

   return
%
% ---------------------------------------------------------------------
