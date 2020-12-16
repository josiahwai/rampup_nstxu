  function [cphi,ipfil,izg,irg]= ...
           makelpcur(zgg,rgg,ip,zp,rp,ap,kap,limod,nfil,iaddctr)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX: (function)
%[cphi,ipfil,izg,irg]=makelpcur(zgg,rgg,ip,zp,rp,ap,kap,limod,nfil,iaddctr)
%
%  PURPOSE:  
%   Function to make elliptical cross-section, flat to parabolic plasma 
%   current distr. defined on standard EPGenv grid. 
%
%  INPUTS:
%    zgg = array of vertical position of grid points on nzxnr grid [m]
%    rgg = array of maj. radial position of grid points on nzxnr grid [m]
%    ip = total plasma current = sum(sum(cphi))
%           [conventionally MA, but can be any units]
%    zp = current centroid vertical position [m]
%    rp = current centroid major radial position [m]
%    ap = minor radius at plasma midplane [m]
%    kap = elongation of elliptical current distribution
%    limod = (optional) specifies the distribution type & peakedness of profile:
%             if limod>0 (continuous distribution):
%		limod=0.5 -> flat profile, limod>=1.0 -> parabolic profile.
%               Default = 1.0 (parabolic distribution)
%             if limod<=0 (discrete filamentary distribution, MFIT-like):
%               limod=0 -> filament at nearest grid point to rp,zp
%               limod=-0.5 -> filaments at 0.5ap, kap ellipse.
%               limod=-1 -> filaments at edge of ap,kap ellipse.
%    nfil = (optional; required if limod<=0) number of filaments. Should be
%            even number >= 4 to have 2 filaments on midplane, equal number 
%            above and below midplane. However, algorithm will work with any
%            number.
%    iaddctr = (optional; use only if limod<0) if =1, adds filament at grid
%            point nearest (zp,rp) location.
%
%  OUTPUTS:
%    cphi = plasma current distribution on nzxnr grid [MA, but see ip above]
%           Note that cphi consists of current filaments, NOT CURRENT DENSITY!!
%            Thus ip = total(cphi);
%    ipfil = vector of currents in each filament
%    irg = vector of radial indices corresponding to locations of current fils
%          in order of ipfil (indices are indices of nrx1 rg vector)
%    izg = vector of vertical indices corresponding to locations of current fils
%          in order of ipfil (indices are indices of nzx1 zg vector)
%
%  RESTRICTIONS:
%     Must have loaded DIII-D geometry environment (>> load_d3denv) or
%       otherwise defined rgg, zgg. nfil should be even number >=4 to have
%       up-down symmetric distribution with 2 filaments on plasma midplane.
%       However, algorithm will work with any number. Note that the final 
%       number of filaments may not equal nfil if the original distribution
%       of nfil filaments is too dense for the grid. Note also that
%       limod=0 will produce a one-filament distribution regardless of nfil
%       (since the nfil points will all be located at rp,zp, resulting in a
%       single grid point being used for the filament.
%
%  METHOD:  
%     Taken from IDL makelpcur.pro. 
%     Filament ring algorithm (limod<=0) developed after Matlab version.
%       For this case, the filaments are placed in the grid points nearest to
%       the originally calculated locations on the ellipse. The current in each
%       filament is simply ip/nfil, so that the resulting current centroid 
%       is not necessarily at (rp,zp). Should fix this at some point... DAH
%     Another option is to make limod (when negative) tell it what fraction of
%      current goes in the central filamnet relative to the edge ones.  DAH

%  WRITTEN BY:  Dave Humphreys 	ON	8/96
%
%  MODIFICATION HISTORY:
%     Adding routine to allow production of ring of filaments (MFIT-like)
%         representing the plasma current.  DAH  10/97
%     Fixing logic for limod>0 so defines cphi correctly for limod ~= 1
%	  (parabolic default)  DAH  1/07
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prelims:
   mu0 = 0.4*pi;

% Derived Values:


% Define Current Distribution (multifil now  5/31/91):
    if ~exist('zgg'), disp('D3D Env not loaded: RETURNING.'), return, end
    ss=size(zgg); nz=ss(1); nr=ss(2);
    bp = kap*ap;
    x2 = (rgg-rp).^2;   z2 = (zgg-zp).^2;
    parab = ones(size(zgg)) - x2./ap.^2 - z2./bp.^2;
    plasma=zeros(nz,nr); 
    if nargin>7
     if limod>0   %parabolic to flat CONTINUOUS distribution
       f = 2*limod-1;  %fraction corr to limod: limod=0.5->1 corr to f=0->1
       %Old definition (< 1/2/07) of cphi for limod>0:
        %cphi=(1-f)*ones(size(zgg)) + f*(ones(size(zgg))-x2./ap.^2-z2./bp.^2);
       %New definition of cphi for limod>0:
       cphi = parab.^f; %f=0 gives flat cphi, f=1 gives parabolic...
       itmp = find(parab>=0);   %defines edge of elliptical region
       plasma(itmp) = cphi(itmp);
       cphi = plasma;
     else  %discrete filamentary distribution
       if nargin<8, disp('No rho argument: RETURNING!'), return, end
       a=abs(limod)*ap; b=kap*a;       
       thetas = linspace(0,2*pi,nfil+1)';
       thetas=thetas(1:nfil);
       if limod==0
         rs=0; 
       else
         rs = sqrt(a^2*b^2./(a^2*sin(thetas).^2 + b^2*cos(thetas).^2));
       end
       xs = rp + rs.*cos(thetas);
       ys = zp + rs.*sin(thetas);
       if nargin>9
         if iaddctr==1, xs=[xs; rp]; ys=[ys; zp]; end
       end
       cphi=zeros(nz*nr,1); rggv=rgg(:); zggv=zgg(:);
       for ii=1:length(xs)
         [tmp,itmp]=min((rggv-xs(ii)).^2 + (zggv-ys(ii)).^2);
         cphi(itmp(1)) = 1;
       end
       cphi = reshape(cphi,nz,nr);
     end
    else
     cphi = parab;
     itmp = find(parab>=0);
     plasma(itmp) = cphi(itmp);
     cphi = plasma;
    end
    cphi = ip*cphi/sum(cphi(:)); %Normlz:Note cphi is currfils, not dens

% Determine rg, zg indices of filament locations (ir, iz):
   irgg=ones(nz,1)*(1:nr); irgg=irgg(:);
   izgg=(1:nz)'*ones(1,nr); izgg=izgg(:);
   cphi1 = cphi(:);
   itmp = find(cphi1~=0);          
   irg = irgg(itmp);
   izg = izgg(itmp);
   ipfil = cphi1(itmp);


