  function [b]= crect2bprobe(rc,zc,dr,dz,cur,rp,zp,theta,leng)

% PURPOSE:
% Generated B-field at magnetic probe located at rp,zp with angle theta and
% length l from circular coil locatd at rc,zc with rectanular width and height
% dr and dz and carrying current cur.
% Note: All Units are MKS (meters, Amps, Tesla)
%
% SYNTAX:
%        [b]= crect2bprobe(rc,zc,dr,dz,i,r,z,theta,leng)
%
% INPUT:
%       rc,zc=        coordinates of center of rectangular coil
%       dr,dz=        width and height of coil
%       cur=            current in coil (+=counter clockwise)
%       rp,zp=        coordinates of the center of the probe
%       angle=        Angle of the coil with the horozontal (deg)
%       leng=         length of the probe in the theta direction
%
% OUTPUT:
%       [b] ==        B field at probes
%
% CAUTION: Present version does not use leng. It computes field only at center
%
% NOTE: uses magnetics routing fine.m which is single precision
%

% Jim Leuer, General Atomics, 2-98

% ---------------
% initialization stuff


  if nargin <= 7
    prog_name= 'crect2bprobe';
    disp(['%ERROR ',prog_name,': Must have at least 7 arguments']);
    eval(['help ',prog_name,]) % print out help from prog_name
    return
  elseif nargin <=8
    leng=zeros(size(rp));
  end

% ---------------------------------------------------------------------------
% compute magnetics using fine assuming length of probe isnt important

  nc= length(rc);
  np= length(rp);
  tr= pi/180*theta;
  ct= cos(tr);
  cs= sin(tr);
  amu= 4e-7*pi;

  for ic= 1:nc

    ri= rc(ic) - 0.5*dr(ic);
    ro= ri     + dr(ic);
    zl= zc(ic) - 0.5*dz(ic);
    zu= zl     + dz(ic);
    ii= cur(ic);
    
    ri= ri*ones(np,1);
    ro= ro*ones(np,1);
    zl= zl*ones(np,1);
    zu= zu*ones(np,1);
    ii= ii*ones(np,1);    
    [hr,hz,hfl]= fine(ri,ro,zl,zu,ii,rp,zp);
    b(:,ic)= amu*(hr.*ct + hz.*cs);
  end % for ic

  return

% --------------------
% check out routine below

  nargin= 8
  rc=[1;2]; zc=[0;1]; dr=[0;1]; dz=[0;1]; cur=[1;1];
  rp=[0;1;1]; zp=[0;0;-1]; theta=[90;0;-90]; leng=[0;0];
%  crect2bprobe
  [b]= crect2bprobe(rc,zc,dr,dz,cur,rp,zp,theta,leng)
                                            
