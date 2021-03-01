%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gs_field_topology
%
%  PURPOSE: Find important nulls in the flux and trace the boundary
%
%  INPUTS: psizr, flux on grid
%
%  OUTPUTS:  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander  ON	2/20/14
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Magnetic axis
gs_update_maxis
if ~maxis
  gs_find_maxis
end

gs_nulls

gs_touches

n = nulls.count;
m = touches.count;
mm = touches.bdefpossible;
mm(m+1:end) = logical(0);
p = n+sum(mm);
bdefs.count = p;
bdefs.r(1:p) = [nulls.r(1:n) touches.r(mm)];
bdefs.z(1:p) = [nulls.z(1:n) touches.z(mm)];
bdefs.psi(1:p) = [nulls.psi(1:n) touches.psi(mm)];
bdefs.w(1:p,:) = [nulls.w(1:n,:); touches.w(mm,:)];
bdefs.ii(1:p,:) = [nulls.ii(1:n,:); touches.ii(mm,:)];
bdefs.drdpsi(1:p,:) = [nulls.drdpsi(1:n,:); touches.drdpsi(mm,:)];
bdefs.dzdpsi(1:p,:) = [nulls.dzdpsi(1:n,:); touches.dzdpsi(mm,:)];
bdefs.limiter(1:p) = [zeros(1,n) touches.limiter(mm)];
bdefs.placursign(1:p) = [zeros(1,n) touches.placursign(mm)];
[~, kk(1:p)] = sort(-sign(cpasma)*bdefs.psi(1:p));

trace_status = -1;
ibdef = 0;
while trace_status ~= 0 & ibdef < bdefs.count
  ibdef = ibdef+1;
  % SHOULD PUT TEST HERE OF THE BDEF POINT
  rbdef = bdefs.r(ibdef);
  zbdef = bdefs.z(ibdef);
  wb = bdefs.w(ibdef,:);
  iib = bdefs.ii(ibdef,:)';
  drbdefdpsi = bdefs.drdpsi(ibdef,:);
  dzbdefdpsi = bdefs.dzdpsi(ibdef,:);
  psibry = bdefs.psi(ibdef);
  
  xwannabe.count = 0;
  for i = 1:bdefs.count
    if i ~= ibdef & bdefs.limiter(i) == 0
      xwannabe.count = xwannabe.count+1;
      n = xwannabe.count;
      xwannabe.r(n) = bdefs.r(i);
      xwannabe.z(n) = bdefs.z(i);
      xwannabe.ii(n,:) = bdefs.ii(i,:);
      xwannabe.drdpsi(n,:) = bdefs.drdpsi(i,:);
      xwannabe.dzdpsi(n,:) = bdefs.dzdpsi(i,:);
    end
  end
  
  lim = bdefs.limiter(ibdef);
  gs_trace_boundary
  if trace_status == 0
    % Check that all of the boundary with the possible exception
    % of a touch point is inside the limiter
    if ~all(isinpoly(rbbbs(1:nbbbs),zbbbs(1:nbbbs)) | gbbbs(1:nbbbs) == 0);
      trace_status = 3;
    end
  end

  if ibdef == bdefs.count & trace_status ~= 0
    % This is a serious error caused by failure of gs_find_nulls
    % or by large "correction" of psizr that puts phony nulls in the plasma
    psicorr1 = 0; % psicorr1 is known to be too small
    psimag = wa*psizr(iia);
    psicorr2 = 0.999*(psimag-psibry); % psicorr2 is known to be large enough
    ix = -1;
    rx = nan;
    wx = nan(1,16);
    drxdpsi = nan(1,16);
    dzxdpsi = nan(1,16);
    psix = nan;
    while abs(psicorr2-psicorr1) > 1e-9
      psicorr = (psicorr1+psicorr2)/2;
      psibry = bdefs.psi(ibdef)+psicorr;
      gs_trace_boundaryx
      if trace_status == 0
        psicorr2 = psicorr;
      else
        psicorr1 = psicorr;
      end
      if plotit > 2
        gs_plot_progress
        %title(hs(4),'Special Tracing','color','red','fonts',24)
        drawnow
      end
    end
    psibry = bdefs.psi(ibdef)+psicorr2;
    gs_trace_boundaryx
    if plotit > 2
      gs_plot_progress
      %title(hs(4),'Special Tracing','color','red','fonts',24)
      drawnow
    end
    rbdef = rx;
    if exist('zx','var')
      zbdef = zx;
    else
      zbdef = 0;
    end
    wb = wx;
    if ~exist('iix','var') % Will fix better later
      nbbbs = 0;
    end
    if nbbbs
      iib = iix;
      drbdefdpsi = drxdpsi;
      dzbdefdpsi = dzxdpsi;
      psibry = psix;
    end
  end
end

% Remove bbbs points that are too close to each other
keep_bbbs = logical(ones(nbbbs_max,1));
r = rbbbs(1);
z = zbbbs(1);
k = 1; % k will be the last approved bbbs point
must_remove_bbbs_point = false;
for i = 2:nbbbs
  if ((r-rbbbs(i))/dr)^2+((z-zbbbs(i))/dz)^2 < 1e-8 % If true then remove
    if gbbbs(i) == 0 % Let's always keep the boundary-defining point
      keep_bbbs(k) = false;
      r = rbbbs(i);
      z = zbbbs(i);
      k = i;
    else
      keep_bbbs(i) = false;
    end
    must_remove_bbbs_point = true;
  else % update r, z
    r = rbbbs(i);
    z = zbbbs(i);
    k = i;
  end
end
if must_remove_bbbs_point
  old_nbbbs = nbbbs;
  nbbbs = sum(keep_bbbs(1:nbbbs));
  gbbbs(1:nbbbs)      = gbbbs(keep_bbbs(1:old_nbbbs));
  ibbbs(1:nbbbs)      = ibbbs(keep_bbbs(1:old_nbbbs));
  rbbbs(1:nbbbs)      = rbbbs(keep_bbbs(1:old_nbbbs));
  zbbbs(1:nbbbs)      = zbbbs(keep_bbbs(1:old_nbbbs));
  wbbbs(1:nbbbs,:)    = wbbbs(keep_bbbs(1:old_nbbbs),:);
  thbbbs(1:nbbbs)     = thbbbs(keep_bbbs(1:old_nbbbs));
  rhobbbs(1:nbbbs)    = rhobbbs(keep_bbbs(1:old_nbbbs));
  dpsibbbsdr(1:nbbbs) = dpsibbbsdr(keep_bbbs(1:old_nbbbs));
  dpsibbbsdz(1:nbbbs) = dpsibbbsdz(keep_bbbs(1:old_nbbbs));
end

kr0 = min(nr-3,max(1,floor((rbdef-rg(1))/dr)));
kz1 = min(nz-2,max(2,ceil((zbdef-zg(1))/dz)));
k = kr0*nz+kz1;
iib = k+neighbors'; % iia indexes 16 points around magnetic axis
pp = psizr(iib);
tr = (rbdef-rgg(k))/dr;
tz = (zbdef-zgg(k))/dz;
wb = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
wbr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
wbz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
wbrr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]/dr^2*mx,1,16);
wbzz = reshape(([0 0 2 6*tz]/dz^2*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
wbrz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);

