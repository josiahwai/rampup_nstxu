  function [br,bz]= cir_filament(rc,zc,ic,rp,zp)
%
% Computes the field: [br,bz] at field point rp,zp
% from a circular filament located at rc,zc with current ic
%
% usage:   [br,bz]= cir_filament(rc,zc,ic,rp,zp)
% where:
%		rc=  radius of filament [m]
%		zc=  z position of filament [m]
%		ic=  current in filament [A]
%		rp=  radial position of field point [m]
%		zp=  z position of field point [m]
%		br=  radial field at field point [T]
%		bz=  axial field at field point [T]
%

% Jim Leuer 6-95 base on elliptic integral formulation

%  idx=  find( rp < eps );
%  rp(idx)= eps*ones(length(idx),1); % puts small number in to prevent devide by zero in br

  dz=    zp-zc;
  dz2=   dz.^2;
  rc2=   rc.^2;
  rp2=   rp.^2;
  den=   (rc+rp).^2 + dz2;
  den2=  sqrt(den);
  k2=    rc.*rp./den*4;
  [k,e]= ellipke(k2);
  bot=   (rp-rc).^2 + dz2;
  br=   -ic.*dz./(rp.*den2).*( k - e.*(rp2+rc2+dz2)./bot )*2.0e-7; % ?-? 030403
  bz=    ic./den2          .*( k - e.*(rp2-rc2+dz2)./bot )*2.0e-7;

  idx=  find(isnan(br)==1);      % find not_a_number caused by rp=0 (devide by zero)
  br(idx)= zeros(length(idx),1); % put in zero for all NaN

  idx=  find(isinf(br)==1);      % find infinite caused by rp=0 (devide by zero)
  br(idx)= zeros(length(idx),1); % put in zero for all inf

  return

% ==========================================================
% test below = Sign Error found 04/03/02 see cirfil.mth for comparison
% error found in sign of br on 030403 when compared with excel basic
% now corrected - look for other problems in other mag/ codes??
  rc= [1 1 4 4];
  zc= [.2 2 .2 2]; 
  ic= [1 1 1 1];

  rp= [3 3 3 3]
  zp= [1 1 1 1]

 [br,bz]= cir_filament(rc,zc,ic,rp,zp);

  del= .001;  
  ri= rc*(1-del);
  ro= rc*(1+del);
  zl= zc - del*rc;
  zu= zc + del*rc;
  cur= 4*pi*1e-7*ones(size(ri));

  [br1,bz1,fl1]= fine(ri,ro,zl,zu,cur,rp,zp); % note cur contains muo

  [br',br1'] % this has sign problem??
  [bz',bz1']
