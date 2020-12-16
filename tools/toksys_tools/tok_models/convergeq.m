function [ceq, eqx] = convergeq(eq,tok_data_struct,options,idoplot)
%
%  USAGE:   [ceq, eqx] = convergeq(eq,tok_data_struct,options,idoplot)
%
%  PURPOSE: Converge equilibrium, i.e. minimize psizr_app+psizr_pla - psizr
%           where psizr_app is calculated from mpc, cc and psizr_pla from mpp, jphi
%           and jphi in turn is calculated using psizr, pprime, ffprim
%
%           The converged equilibrium is on the grid in tok_data_struct
%           which can be different from the grid in the original efit.
%
%           Bus constraint such as the VFI on DIII-D can be imposed on cc
%
%  INPUTS:  eq: An equilibrium containing at least psizr, pprime, ffprim, cc
%           tok_data_struct: standard toksys object for tokamak description
%           options.converrmax = upper limit on resulting flux error
%             max(abs(psizr_app(:)+psizr_pla(:)-psizr(:))), default 1e-12
%           options.iconstraints: (default is 3)
%             0. conserve pres           and fpol as functions of normalized flux
%             1: conserve pres           and total plasma current
%             2: conserve thermal energy and fpol
%             3: conserve thermal energy and total plasma current
%           options.maxiter: maximum number of iterations, default = 25;
%           options.icx: Specify a subset of cc that can change, default icx=1:nc
%           options.cccirc: standard toksys assignment of circuit numbers to coils,
%           options.bus_code: indices are 1 for coils on bus, 0 for others
%              bus_code*cc = 0 (bus_code is row vector and cc is column vector)
%              cc will be adjusted to meet the constraint: bus_code*cc = 0
%              For DIII-D: PP_objs =  get_PP_objs(shot)
%                          options.bus_code = [0 0 PP_objs.bus_code]              
%           idoplot: plot each iteration for idoplot seconds, default 0 = no plot
%
%  OUTPUTS: ceq: structure containing converged equilibrium on tok_data_struct grid
%           eqx: structure containing extra information (see eqx.descriptions)
%	
%  RESTRICTIONS: The following quantities are updated for converged equilibrium:
%    		 rbbbs, zbbbs, jphi, psizr, pprime, ffprim, pres, fpol
%                cc, psimag, psibry, cpasma, rmaxis, zmaxis, rg, zg, dr, dz
%                nw, nh, psirz, ssibry, ssimag, qpsi
%
%  METHOD:  The Newton-Rhapson method is used to converge to the solution.

%           In each iteration the equilibrium is corrected by response 
%           to flux error and adjustment of cc.
%           The adjustments of cc strive to preserve the boundary without
%           excessive changes of cc. This is achieved by minimizing norm of: 
%           [flux errors at boundary points; change of cc]
%           When nulls are close to plasma or are defining boundary then 
%           flux differences between these and any touch point are controlled
%           with a factor 10 extra weight
	
%  VERSION @(#)convergeq.m	1.17 12/05/13
%
%  WRITTEN BY:  Anders Welander  ON	11/19/10
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  mu0 = .4e-6*pi;
  twopi = 2*pi;
  
  if nargin == 0
    help convergeq
    return
  end
  
  if ~exist('options','var')
    options.converrmax = 1e-12;
  end
  if ~isfield(options,'converrmax')
    options.converrmax = 1e-12;
  end
  if ~isfield(options,'maxiter')
    options.maxiter = 25;
  end
  if ~isfield(options,'nxiter') % Number of extra iterations after converrmax achieved
    options.nxiter = 0;
  end
  if ~isfield(options,'iconstraints')
    options.iconstraints = 3;
  end
  if ~isfield(options,'iverbose')
    options.iverbose = 0;
  end
  if ~exist('idoplot','var')
    idoplot = 0;
  end
  
  % Unpack the equilibrium, eq
  struct_to_ws(eq);
  % Return directly if this is garbage
  if max(abs(jphi(:))) == 0 | ~isempty(find(isnan(jphi(:))))
    ceq = eq;
    eqx.convFailure = 1;
    return
  end
  na = nbbbs-1;
  rbbbs = rbbbs(1:nbbbs);
  zbbbs = zbbbs(1:nbbbs);
  if isfield(options,'Iptarget')
    Iptarget = options.Iptarget;
  else
    Iptarget = cpasma;
  end  
  psigrid = linspace(psimag,psibry,nw)'; % Wb
  dpsigrid = psigrid(2)-psigrid(1);
  rgg = ones(nh,1)*rg(:)';
  iplasma = find(jphi);
  c0 = spline(1:nw,pprime,1+1/64:1/32:nw);
  pres(nw,1) = 0;
  for k = nw-1:-1:1
    pres(k) = pres(k+1) - sum(c0(k*32-31:k*32))*dpsigrid/64/pi;
  end
  P = spline(psigrid,pres,psizr);
  Wth0 = 3*pi*rgg(iplasma)'*P(iplasma)*dr*dz;
  if isfield(options,'Wtarget')
    Wtarget = options.Wtarget;
  else
    Wtarget = Wth0;
  end
 
  % Unpack tok_data_struct
  struct_to_ws(tok_data_struct);
  dr = mean(diff(rg)); dz = mean(diff(zg)); ngg = nr*nz;
  
  % Compare grids and rescale if needed
  griderror = nw ~= nr | nh ~= nz | (max(abs([eq.rg(1)-rg(1) eq.rg(end)-rg(end) eq.zg(1)-zg(1) eq.zg(end)-zg(end)])) > 0.001);
  if griderror % Convert eq to grid in tok_data_struct
    psizr = interp2(eq.rg,eq.zg,psizr,rgg,zgg,'spline');
    psibarzr = (psizr-psimag)/(psibry-psimag);
    pprime = spline(linspace(0,1,nw),pprime,linspace(0,1,nr))';
    ffprim = spline(linspace(0,1,nw),ffprim,linspace(0,1,nr))';
    pres =   spline(linspace(0,1,nw),pres,linspace(0,1,nr))';
    qpsi =   spline(linspace(0,1,nw),qpsi,linspace(0,1,nr))';
    delstarpsizr = (-[zeros(nz,1) psizr(:,3:end)-psizr(:,1:end-2) zeros(nz,1)]/2/dr./rgg+...
      [zeros(nz,1) psizr(:,3:end)+psizr(:,1:end-2)-2*psizr(:,2:end-1) zeros(nz,1)]/dr^2+...
      [zeros(1,nr);psizr(3:end,:)+psizr(1:end-2,:)-2*psizr(2:end-1,:);zeros(1,nr)]/dz^2)/2/pi;
    jphi = -delstarpsizr/mu0./rgg/1e6;
    ivac = find(zgg<min(zbbbs) | zgg>max(zbbbs) | psibarzr>1);
    jphi(ivac) = 0;
    iplasma = find(jphi);
  end

  % Remember original pprime, ffprim, psizr
  pprime0 = pprime;
  ffprim0 = ffprim;
  psizr0 = psizr;
  
  % Create structure: converged equilibrium
  ceq = eq; % ceq will become Converged equilibrium if convergence succeeded
  
  % Store higher resolution limiter in rl, zl
  if size(limdata,2)<size(limdata,1), limdata=limdata'; end
  rlimdata = limdata(2,:); zlimdata = limdata(1,:);
  if rlimdata(1)~=rlimdata(end) | zlimdata(1)~=zlimdata(end)
    rlimdata = rlimdata([1:end 1]); zlimdata = zlimdata([1:end 1]);
  end
  nlim = length(rlimdata);
  rl = interp1(1:nlim,rlimdata,1:.01:nlim);
  zl = interp1(1:nlim,zlimdata,1:.01:nlim);
  
  % Get I.vc0t to calculate applied flux from vessel currents
  I = cc_efit_to_tok(tok_data_struct,eq);
  
  % Switch to efit convention for calculating flux from cc
  efmpc = zeros(ngg,length(cc));
  efgbr2c = zeros(ngg,length(cc));
  efgbz2c = zeros(ngg,length(cc));
  qx = eq;
  for j = 1:length(cc);
    qx.cc = cc*0;
    qx.cc(j) = 1;
    Ix = cc_efit_to_tok(tok_data_struct,qx);
    efmpc(:,j) = mpc*Ix.cc0t;
    efgbr2c(:,j) = gbr2c*Ix.cc0t;
    efgbz2c(:,j) = gbz2c*Ix.cc0t;
  end
    
  % The value bus_code_lumped*cc will be controlled to zero.
  if isfield(options,'bus_code')
    bus_code = options.bus_code;
  else
    bus_code = zeros(1,length(cc));
  end
  if sum(abs(bus_code))>0
    bus_code_lumped = mean(ccnturn(find(bus_code)))*bus_code./ccnturn';
  else
    bus_code_lumped = zeros(1,length(cc));
  end
  
  % Only the coils indexed by icx will be adjusted to converge equilibrium
  if isfield(options,'icx')
    icx = options.icx;
  else
    icx = 1:length(cc);
  end
  
  % Transform cc, efmpc, icx, bus_code for circuits if cccirc exists
  if isfield(options,'cccirc')
    cccirc = options.cccirc;
    ncc = max(abs(cccirc));
    efmpc2 = zeros(ngg,ncc);
    efgbr2c2 = zeros(ngg,ncc);
    efgbz2c2 = zeros(ngg,ncc);
    cc2 = zeros(ncc,1);
    bus_code2 = zeros(1,ncc);
    bus_code_lumped2 = zeros(1,ncc);
    icxflags = zeros(1,length(cc)); icxflags(icx) = 1;
    icxflags2 = zeros(1,ncc);
    for j = 1:ncc % j = number of circuit
      k = find(abs(cccirc)==j); % k = numbers of the coils in ciruit j      
      n = length(k);
      for l = k
        s = sign(cccirc(l));
        cc2(j) = cc2(j) + s*cc(l)/n;
        efmpc2(:,j) = efmpc2(:,j) + s*efmpc(:,l);
        efgbr2c2(:,j) = efgbr2c2(:,j) + s*efgbr2c(:,l);
        efgbz2c2(:,j) = efgbz2c2(:,j) + s*efgbz2c(:,l);
	bus_code2(j) = bus_code2(j) + s*bus_code(l);
	bus_code_lumped2(j) = bus_code_lumped2(j) + s*bus_code_lumped(l);
	icxflags2(j) = icxflags2(j) + s*icxflags(l)/n;
      end
    end
    cc = cc2;
    efmpc = efmpc2;
    efgbr2c = efgbr2c2;
    efgbz2c = efgbz2c2;
    bus_code = bus_code2;
    bus_code_lumped = bus_code_lumped2;
    icx = find(icxflags2);
  end
  
  dcc = cc*0; % dcc will hold boundary-stabilizing adjustment of cc
  ncx = length(icx);

  % Grid points that can be considered in the plasma have a maximum initial rho compared to bbbs
  rhobbbs = sqrt((rbbbs(2:end)-rmaxis).^2+(zbbbs(2:end)-zmaxis).^2);
  thbbbs = angle((rbbbs(2:end)-rmaxis)+sqrt(-1)*(zbbbs(2:end)-zmaxis)); thbbbs = thbbbs(:);
  rhogg = sqrt((rgg-rmaxis).^2+(zgg-zmaxis).^2);
  thgg = angle((rgg-rmaxis)+sqrt(-1)*(zgg-zmaxis));
  rhomaxgg = interp1([thbbbs-2*pi; thbbbs; thbbbs+2*pi],[rhobbbs; rhobbbs; rhobbbs],thgg);
  
  % Find weights in grid points to calculate value at a point using cubic Hermite spline (ref. wikipedia)
  mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2;
  neighbors = reshape([-1-nz -1 -1+nz -1+2*nz;-nz 0 nz 2*nz;1-nz 1 1+nz 1+2*nz;2-nz 2 2+nz 2+2*nz],1,16);
  
  % Find magnetic axis
  for j=1:9
    kr0 = min(nr-3,max(1,floor((rmaxis-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
    kz1 = min(nz-2,max(2,ceil((zmaxis-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
    k = kr0*nz+kz1;
    iia = k+neighbors; % iia indexes 16 points around magnetic axis
    pp = psizr(iia);
    tr = (rmaxis-rgg(k))/dr; tz = (zmaxis-zgg(k))/dz;
    wa = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
    war = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
    waz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
    warr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]/dr^2*mx,1,16);
    wazz = reshape(([0 0 2 6*tz]/dz^2*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
    warz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
    ashift_Ba = -inv([pp*warr' pp*warz';pp*warz' pp*wazz']);
    c = ashift_Ba*[pp*war'; pp*waz']; rmaxis = rmaxis+c(1); zmaxis = zmaxis+c(2);
  end
  psimag = wa*psizr(iia)';

  % Find TOUCH POINT candidate
  psibarzr = (psizr-psimag)/(psibry-psimag);
  psibarlim = interp2(rg,zg,psibarzr,rl,zl,'spline');
  for j=1:length(rl) % Distance check to rbbbs,zbbbs needed
    d2minlimbbbs(j) = min((rbbbs-rl(j)).^2+(zbbbs-zl(j)).^2);
  end
  k = find(zl<max(zbbbs) & zl>min(zbbbs) & d2minlimbbbs<dr^2);
  if isempty(k)
    [dum, k1] = min(d2minlimbbbs);
    psibartouch = psibarlim(k1);
  else
    [psibartouch, l] = min(psibarlim(k)); k1 = k(l);
  end
  rlim = rl(k1); zlim = zl(k1);
  k2 = k1+1; if k2 > length(rl), k2 = 2; end
  ulim = [rl(k2)-rl(k1) , zl(k1)-zl(k2)];
  ulim = ulim/norm(ulim);
  kr0 = min(nr-3,max(1,floor((rlim-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
  kz1 = min(nz-2,max(2,ceil((zlim-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
  k = kr0*nz+kz1;
  iilim = k+neighbors;
  tr = (rlim-rgg(k))/dr; tz = (zlim-zgg(k))/dz;
  wl = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
  wld = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]*mx*ulim(1)/dr+...
                 ([0 1 2*tz 3*tz^2]*mx)'*ulim(2)/dz*[1 tr tr^2 tr^3]*mx,1,16);
		 
  % Check for UPPER NULL
  psizr_r = [zeros(nz,1) psizr(:,3:end)-psizr(:,1:end-2) zeros(nz,1)]/2/dr;
  psizr_z = [zeros(1,nr);psizr(3:end,:)-psizr(1:end-2,:);zeros(1,nr)]/2/dz;
  j = find(zbbbs>zmaxis);
  [dum, k] = min(interp2(rg,zg,psizr_r.^2+psizr_z.^2,rbbbs(j),zbbbs(j),'spline'));
  rx1 = rbbbs(j(k)); zx1 = zbbbs(j(k));
  kr0 = min(nr-3,max(1,floor((rx1-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
  kz1 = min(nz-2,max(2,ceil((zx1-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
  k = kr0*nz+kz1;
  c_previous = zg(end)-zg(1); cx1 = c_previous*.99; j = 9;
  % Try zooming in on x-point with Newton Rhapson.
  while j>0 & k>1+nz & k<=ngg-2-2*nz & cx1<c_previous
    % Find indices and weights for grid points around the x-point
    iix1 = k+neighbors; % iix indexes 16 points around x point
    pp = psizr(iix1); j = j-1; c_previous = cx1;
    tr = (rx1-rgg(k))/dr; tz = (zx1-zgg(k))/dz;
    wx1 = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
    wx1r = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
    wx1z = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
    wx1rr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]/dr^2*mx,1,16);
    wx1zz = reshape(([0 0 2 6*tz]/dz^2*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
    wx1rz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
    x1shift_Bx = -inv([pp*wx1rr' pp*wx1rz';pp*wx1rz' pp*wx1zz']);
    cx1 = x1shift_Bx*[pp*wx1r'; pp*wx1z']; rx1 = rx1+cx1(1); zx1 = zx1+cx1(2);
    kr0 = min(nr-3,max(1,floor((rx1-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
    kz1 = min(nz-2,max(2,ceil((zx1-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
    k = kr0*nz+kz1;
  end
  psibarx1 = psibarzr(iix1)*wx1'; psix1 = psizr(iix1)*wx1';
  if min(sqrt((rbbbs-rx1).^2+(zbbbs-zx1).^2)) > 2*sqrt(dr^2+dz^2) | norm(cx1)>max(dr,dz)/10
    psibarx1 = 1e9;
    ztop = max(zbbbs)+dz;
  else
    ztop = zx1; % Max z for where plasma can be
  end

  % Check for LOWER NULL
  j = find(zbbbs<zmaxis);
  [dum, k] = min(interp2(rg,zg,psizr_r.^2+psizr_z.^2,rbbbs(j),zbbbs(j),'spline'));
  rx2 = rbbbs(j(k)); zx2 = zbbbs(j(k));
  kr0 = min(nr-3,max(1,floor((rx2-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
  kz1 = min(nz-2,max(2,ceil((zx2-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
  k = kr0*nz+kz1;
  c_previous = zg(end)-zg(1);; cx2 = c_previous*.99; j = 9;
  while j>0 & k>1+nz & k<=ngg-2-2*nz & cx2<c_previous
    iix2 = k+neighbors; % iix2 indexes 16 points around 2:nd x point
    pp = psizr(iix2); j = j-1; c_previous = cx2;
    tr = (rx2-rgg(k))/dr; tz = (zx2-zgg(k))/dz;
    wx2 = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
    wx2r = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
    wx2z = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
    wx2rr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]/dr^2*mx,1,16);
    wx2zz = reshape(([0 0 2 6*tz]/dz^2*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
    wx2rz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
    x2shift_Bx = -inv([pp*wx2rr' pp*wx2rz';pp*wx2rz' pp*wx2zz']);
    cx2 = x2shift_Bx*[pp*wx2r'; pp*wx2z']; rx2 = rx2+cx2(1); zx2 = zx2+cx2(2);
    kr0 = min(nr-3,max(1,floor((rx2-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
    kz1 = min(nz-2,max(2,ceil((zx2-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
    k = kr0*nz+kz1;
  end
  psibarx2 = psibarzr(iix2)*wx2'; psix2 = psizr(iix2)*wx2';
  if min(sqrt((rbbbs-rx2).^2+(zbbbs-zx2).^2)) > 2*sqrt(dr^2+dz^2) | norm(cx2)>max(dr,dz)/10
    psibarx2 = 1e9;
    zbot = min(zbbbs)-dz;
  else
    zbot = zx2; % Min z for where plasma can be
  end
  
  % Which point defines boundary? 
  if psibartouch < min(psibarx1,psibarx2)
    ilimited = 1;
    r0 = rlim;
    z0 = zlim;
  else
    ilimited = 0;
    if psibarx1 < psibarx2
      r0 = rx1;
      z0 = zx1;
    else
      r0 = rx2;
      z0 = zx2;
    end
  end
  kr0 = min(nr-3,max(1,floor((r0-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
  kz1 = min(nz-2,max(2,ceil((z0-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
  k = kr0*nz+kz1;
  ii0 = k+neighbors; % ii0 indexes 16 points around limiting point
  tr = (r0-rgg(k))/dr; tz = (z0-zgg(k))/dz;
  w0 = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
  
  % Index of boundary defining point in ibdef
  [dum, ibdef] = min((rbbbs-r0).^2+(zbbbs-z0).^2);
  
  % If there are nulls close to plasma control their locations
  if min(sqrt((rbbbs-rx1).^2+(zbbbs-zx1).^2)) < 3*dz & psibarx1 < 1.2
    icontrolNULL1 = 1;
  else
    icontrolNULL1 = 0;
  end
  if min(sqrt((rbbbs-rx2).^2+(zbbbs-zx2).^2)) < 3*dz & psibarx2 < 1.2
    icontrolNULL2 = 1;
  else
    icontrolNULL2 = 0;
  end

  % Calculate response matrix Minv
  psigrid = linspace(psimag,psibry,nr)'; % Wb
  Pprime = spline(psigrid,pprime,psizr);
  FFprim = spline(psigrid,ffprim,psizr);
  iplasma = find(jphi);
  ivac = find(jphi==0);
  P = spline(psigrid,pres,psizr);
  pb = (spline(psigrid,pprime,psigrid+1e-6)-spline(psigrid,pprime,psigrid-1e-6))*pi*1e6;
  Pbis = spline(psigrid,pb,psizr);
  gb = (spline(psigrid,ffprim,psigrid+1e-6)-spline(psigrid,ffprim,psigrid-1e-6))*pi*1e6;
  Gbis = spline(psigrid,gb,psizr);
  djdp = (rgg.*Pbis+Gbis/mu0./rgg)/twopi; djdp(ivac)=0;
  % Some variables to facilitate the making of M
  x = psibarzr(iplasma); a = 1-x; R = rgg(iplasma); % For more readability below
  jpla = 1e6*jphi(iplasma)/(psimag-psibry); djdpsi = djdp(iplasma);
  coupling = iplasma';
  ilt = find(jphi(1+nz:ngg)~=0 & jphi(1:ngg-nz)==0)+nz;
  irt = find(jphi(1:ngg-nz)~=0 & jphi(1+nz:ngg)==0);
  idn = find(jphi(1+1:ngg)~=0 & jphi(1:ngg-1)==0)+1;
  iup = find(jphi(1:ngg-1)~=0 & jphi(1+1:ngg)==0);
  psizr_r = [zeros(nz,1) psizr(:,3:end)-psizr(:,1:end-2) zeros(nz,1)]/2/dr;
  psizr_z = [zeros(1,nr);psizr(3:end,:)-psizr(1:end-2,:);zeros(1,nr)]/2/dz;
  diltdpsi = +(rgg(ilt)*pprime(end)+ffprim(end)/mu0./rgg(ilt))./psizr_r(ilt)*dz;
  dirtdpsi = -(rgg(irt)*pprime(end)+ffprim(end)/mu0./rgg(irt))./psizr_r(irt)*dz;
  didndpsi = +(rgg(idn)*pprime(end)+ffprim(end)/mu0./rgg(idn))./psizr_z(idn)*dr;
  diupdpsi = -(rgg(iup)*pprime(end)+ffprim(end)/mu0./rgg(iup))./psizr_z(iup)*dr;
  iused = union(iplasma,ii0);
  for k = 1:na
    kr0 = min(nr-3,max(1,floor((rbbbs(k)-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
    kz1 = min(nz-2,max(2,ceil((zbbbs(k)-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
    bi(k) = kr0*nz+kz1;
    tr = (rbbbs(k)-rgg(bi(k)))/dr; tz = (zbbbs(k)-zgg(bi(k)))/dz;
    bw(k,:) = reshape(mx'*[1 tz tz^2 tz^3]'*[1 tr tr^2 tr^3]*mx,1,16);
    iused = union(iused,bi(k)+neighbors);
  end
  nused = length(iused);
  inotused = setdiff(1:ngg,iused);
  iplasmac = find(ismember(iused,iplasma));
  nplasmac = length(iplasmac);
  iiac = find(ismember(iused,iia));
  ii0c = find(ismember(iused,ii0));
  iltc = find(ismember(iused,ilt));
  irtc = find(ismember(iused,irt));
  idnc = find(ismember(iused,idn));
  iupc = find(ismember(iused,iup));
  % Equations for flux change on grid.
  M = eye(nused); % M = dpsizr_app/dpsizr
  for j = 1:nused
    iru = ceil(iused(j)/nz);
    for k=1:length(iplasma)
      irp = ceil(iplasma(k)/nz);
      idz = 1+abs(iused(j)-iplasma(k)+(irp-iru)*nz);
      coupling(k) = mpp(idz+(iru-1)*nz,irp);
    end
    for k=1:length(ilt)
      irp = ceil(ilt(k)/nz);
      idz = 1+abs(iused(j)-ilt(k)+(irp-iru)*nz);
      mutlt(k) = mpp(idz+(iru-1)*nz,irp);
    end
    for k=1:length(irt)
      irp = ceil(irt(k)/nz);
      idz = 1+abs(iused(j)-irt(k)+(irp-iru)*nz);
      mutrt(k) = mpp(idz+(iru-1)*nz,irp);
    end
    for k=1:length(idn)
      irp = ceil(idn(k)/nz);
      idz = 1+abs(iused(j)-idn(k)+(irp-iru)*nz);
      mutdn(k) = mpp(idz+(iru-1)*nz,irp);
    end
    for k=1:length(iup)
      irp = ceil(iup(k)/nz);
      idz = 1+abs(iused(j)-iup(k)+(irp-iru)*nz);
      mutup(k) = mpp(idz+(iru-1)*nz,irp);
    end
    % Change of integrand in area integral
    M(j,iplasmac) = M(j,iplasmac)-dr*dz*coupling.*djdpsi';
    M(j,iiac) = M(j,iiac)+dr*dz*coupling*(a.*djdpsi+jpla)*wa;
    M(j,ii0c) = M(j,ii0c)+dr*dz*coupling*(x.*djdpsi-jpla)*w0;
    M(j,nused+1) = -dr*dz*coupling*(R.*Pprime(iplasma));
    M(j,nused+2) = -dr*dz*coupling*(FFprim(iplasma)/mu0./R);
    % Boundary displacement contribution (change of area in area integral)
if 0
    M(j,iltc) = M(j,iltc)-mutlt.*diltdpsi;
    M(j,irtc) = M(j,irtc)-mutrt.*dirtdpsi;
    M(j,idnc) = M(j,idnc)-mutdn.*didndpsi;
    M(j,iupc) = M(j,iupc)-mutup.*diupdpsi;
    M(j,ii0c) = M(j,ii0c)+(mutlt*diltdpsi'+mutrt*dirtdpsi'+mutdn*didndpsi'+mutup*diupdpsi')*w0;
end
  end
  % Perturbation of thermal energy
  M(nused+1,iplasmac) = dr*dz*3/2*R'.*Pprime(iplasma)';
  M(nused+1,iiac) = M(nused+1,iiac)-dr*dz*3/2*sum(a.*R.*Pprime(iplasma))*wa;
  M(nused+1,ii0c) = M(nused+1,ii0c)-dr*dz*3/2*sum(x.*R.*Pprime(iplasma))*w0;
  M(nused+1,nused+1) = Wth0;
  % Perturbation of total plasma current
  M(nused+2,iplasmac) = dr*dz*djdpsi';
  M(nused+2,iiac) = M(nused+2,iiac) - dr*dz*sum(a.*djdpsi+jpla)*wa;
  M(nused+2,ii0c) = M(nused+2,ii0c) - dr*dz*sum(x.*djdpsi-jpla)*w0;
  M(nused+2,nused+1) = dr*dz*sum(R.*Pprime(iplasma));
  M(nused+2,nused+2) = dr*dz*sum(FFprim(iplasma)/mu0./R);

  % Indices of equations and variables to use
  if     options.iconstraints == 0
    ii = 1:nused;
  elseif options.iconstraints == 1
    ii = [1:nused nused+2];
  elseif options.iconstraints == 2
    ii = 1:nused+1;
  elseif options.iconstraints == 3
    ii = 1:nused+2;
  end
  C = zeros(nused+2,1);

  % Solve equation for perturbed equilibrium
  Minv = inv(M(ii,ii));
  dpsizr = zeros(nz,nr);

  % Create dbdcc, where b is vector of boundary points and cc is the efit coil current vector
  dbdcc = zeros(na,length(cc)); db = zeros(na,1);
  % Control flux on b to psibry and, if diverted, field at x-point to 0
  for j = 1:length(cc)
    C(1:nused) = efmpc(iused,j);
    pert = Minv*C(ii);
    dpsizr(iused) = pert(1:nused);
    for k = 1:na
      dbdcc(k,j) = bw(k,:)*dpsizr(bi(k)+neighbors')-w0*dpsizr(ii0');
    end
    if ilimited
      dbdcc(na+1,j) = wld*dpsizr(ii0');
      k = k+1;
    end
    if icontrolNULL1
      dbdcc(k+1,j) = wx1r*dpsizr(iix1')*10;
      dbdcc(k+2,j) = wx1z*dpsizr(iix1')*10;
      k = k+2;
    end
    if icontrolNULL2
      dbdcc(k+1,j) = wx2r*dpsizr(iix2')*10;
      dbdcc(k+2,j) = wx2z*dpsizr(iix2')*10;
      k = k+2;
    end
    if icontrolNULL1 & icontrolNULL2
      dbdcc(k+1,j) = (wx2*dpsizr(iix2')-wx1*dpsizr(iix1'))*10;
      k = k+1;
    end
    if ilimited & icontrolNULL1
      dbdcc(k+1,j) = (wx1*dpsizr(iix1')-w0*dpsizr(ii0'))*10;
      k  = k+1;
    end
    if ilimited & icontrolNULL2
      dbdcc(k+1,j) = (wx2*dpsizr(iix2')-w0*dpsizr(ii0'))*10;
      k  = k+1;
    end    
  end

  % Begin iterations
  iterations = 0;
  convFailure = 0;
  converrors = [];
  xiter = -1;
  while xiter < options.nxiter & iterations < options.maxiter & convFailure == 0
    
    iterations = iterations+1;
    
    % Update rmaxis, zmaxis, psimag  
    for j=1:9
      kr0 = min(nr-3,max(1,floor((rmaxis-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
      kz1 = min(nz-2,max(2,ceil((zmaxis-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
      k = kr0*nz+kz1;
      iia = k+neighbors; % iia indexes 16 points around magnetic axis
      pp = psizr(iia);
      tr = (rmaxis-rgg(k))/dr; tz = (zmaxis-zgg(k))/dz;
      wa = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      war = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
      waz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      warr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]/dr^2*mx,1,16);
      wazz = reshape(([0 0 2 6*tz]/dz^2*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      warz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
      ashift_Ba = -inv([pp*warr' pp*warz';pp*warz' pp*wazz']);
      c = ashift_Ba*[pp*war'; pp*waz']; rmaxis = rmaxis+c(1); zmaxis = zmaxis+c(2);
    end
    psimag = wa*psizr(iia)';
    
    % Update r0, z0, psibry  
    if ilimited
      kr0 = min(nr-3,max(1,floor((r0-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
      kz1 = min(nz-2,max(2,ceil((z0-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
      k = kr0*nz+kz1;
      ii0 = k+neighbors;
      tr = (rlim-rgg(k))/dr; tz = (zlim-zgg(k))/dz;
      wl = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      wld = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]*mx*ulim(1)/dr+...
                     ([0 1 2*tz 3*tz^2]*mx)'*ulim(2)/dz*[1 tr tr^2 tr^3]*mx,1,16);
      wlb = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]*mx*ulim(1)^2/dr^2+...
                     ([0 0 2 6*tz]*mx)'*ulim(2)^2/dz^2*[1 tr tr^2 tr^3]*mx+...
		   2*([0 1 2*tz 3*tz^2]*mx)'*ulim(2)/dz*[0 1 2*tr 3*tr^2]*mx*ulim(1)/dr,1,16);
      dpsitouchdulim = psizr(ii0)*wld';
      d2psitouchdulim2 = psizr(ii0)*wlb';
      r0 = r0-ulim(1)*dpsitouchdulim/d2psitouchdulim2;
      z0 = z0-ulim(2)*dpsitouchdulim/d2psitouchdulim2;
    else    
      for j=1:9
	kr0 = min(nr-3,max(1,floor((r0-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
	kz1 = min(nz-2,max(2,ceil((z0-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
	k = kr0*nz+kz1;
	ii0 = k+neighbors;
	pp = psizr(ii0);
	tr = (r0-rgg(k))/dr; tz = (z0-zgg(k))/dz;
	w0 = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
	w0r = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
	w0z = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
	w0rr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]/dr^2*mx,1,16);
	w0zz = reshape(([0 0 2 6*tz]/dz^2*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
	w0rz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
	xshift_Bx = -inv([pp*w0rr' pp*w0rz';pp*w0rz' pp*w0zz']);
	c = xshift_Bx*[pp*w0r'; pp*w0z']; r0 = r0+c(1); z0 = z0+c(2);
      end
    end
    psibry = w0*psizr(ii0)';
   
    % Update, iplasma, psigrid, pres, ffprim
    iplasma0 = iplasma;
    iplasma = find(zgg < ztop & zgg > zbot & (psizr-psimag)/(psibry-psimag) <= 1 & rhogg-rhomaxgg < dr);
    if length(setdiff(iplasma,iplasma0))>0 | length(setdiff(iplasma0,iplasma))>0
      % Code to switch to Minv without boundary displacement effect
    end
    %iplasma = iplasma0;
    ivac = setdiff(1:ngg,iplasma);
    psigrid = linspace(psimag,psibry,nr)'; % Wb
    dpsigrid = psigrid(2)-psigrid(1);
    P = spline(psigrid,pres,psizr);
    Wth = 3*pi*rgg(iplasma)'*P(iplasma)*dr*dz;    
    if     options.iconstraints == 0 % Keeping pres and fpol
      pprime = pprime0*(eq.psimag-eq.psibry)/(psimag-psibry);
      ffprim = ffprim0*(eq.psimag-eq.psibry)/(psimag-psibry);
    elseif options.iconstraints == 1 % Keeping pres
      pprime = pprime0*(eq.psimag-eq.psibry)/(psimag-psibry);
      c0 = spline(1:nr,ffprim,1+1/64:1/32:nr);
      half_fpol_squared(nr,1) = 0.5*eq.fpol(end)^2;
      for k = nr-1:-1:1
	half_fpol_squared(k) = half_fpol_squared(k+1) - sum(c0(k*32-31:k*32))*dpsigrid/64/pi;
      end
      fpol = sign(eq.fpol(end))*sqrt(half_fpol_squared*2);
    elseif options.iconstraints == 2 % Keeping fpol
      ffprim = ffprim0*(eq.psimag-eq.psibry)/(psimag-psibry);      
      c0 = spline(1:nr,pprime,1+1/64:1/32:nr);
      pres(nr,1) = 0;
      for k = nr-1:-1:1
	pres(k) = pres(k+1) - sum(c0(k*32-31:k*32))*dpsigrid/64/pi;
      end
    elseif options.iconstraints == 3
      c0 = spline(1:nr,pprime,1+1/64:1/32:nr);
      pres(nr,1) = 0;
      for k = nr-1:-1:1
	pres(k) = pres(k+1) - sum(c0(k*32-31:k*32))*dpsigrid/64/pi;
      end
      c0 = spline(1:nr,ffprim,1+1/64:1/32:nr);
      half_fpol_squared(nr,1) = 0.5*eq.fpol(end)^2;
      for k = nr-1:-1:1
	half_fpol_squared(k) = half_fpol_squared(k+1) - sum(c0(k*32-31:k*32))*dpsigrid/64/pi;
      end
      fpol = sign(eq.fpol(end))*sqrt(half_fpol_squared*2);
    end
    
    Pprime = spline(psigrid,pprime,psizr);
    FFprim = spline(psigrid,ffprim,psizr);
    jphi = (rgg.*Pprime+FFprim/mu0./rgg)/1e6; jphi(ivac) = 0;
    psizr_pla = gscalc(jphi,rg,zg,mpp);
    psizr_app = reshape(efmpc*cc+mpv*I.vc0t,nz,nr);
    
    % At this point we have (starting with a psizr) derived:
    % r0,z0, ra, za, psibry, psimag, fpol, ffprim, pres, pprime, jphi, psizr_pla
    % We have also calculated Wth, Ip
    
    % Calculate errors
    psizr_err = psizr_app+psizr_pla - psizr; 
    dW = Wtarget-Wth;
    dI = Iptarget - 1e6*sum(jphi(:))*dr*dz;
    
    % Calculate the dpsizr that removes the errors in psizr-psizr_app-psizr_pla and possibly dW, dI    
    siused = size(iused);
    C = [reshape(psizr_err(iused),siused(1),siused(2)); dW; dI];
    pert1 = Minv*C(ii); % pert1 corrects convergence and W, Ip errors
    
    psizr(iused) = psizr(iused)+reshape(pert1(1:nused),siused(1),siused(2)); % psizr is corrected for the errors but not yet pres, fpol, etc
    
    % Control the boundary with isoflux control at original rbbbs, zbbbs points
    % Convergence is much faster if dcc only corrects for last dpsizr
    dpsizr(iused) = pert1(1:nused)';
    for k = 1:na
      db(k) = bw(k,:)*dpsizr(bi(k)+neighbors')-w0*dpsizr(ii0');
    end
    if ilimited
      db(na+1) = wld*dpsizr(ii0');
      k = k+1;
    end
    if icontrolNULL1
      db(k+1) = wx1r*dpsizr(iix1');
      db(k+2) = wx1z*dpsizr(iix1');
      k = k+2;
    end
    if icontrolNULL2
      db(k+1) = wx2r*dpsizr(iix2');
      db(k+2) = wx2z*dpsizr(iix2');
      k = k+2;
    end
    if icontrolNULL1 & icontrolNULL2
      db(k+1) = (wx2*dpsizr(iix2')-wx1*dpsizr(iix1'))*10;
      k  = k+1;
    end
    if ilimited & icontrolNULL1
      db(k+1) = (wx1*dpsizr(iix1')-w0*dpsizr(ii0'))*10;
      k  = k+1;
    end
    if ilimited & icontrolNULL2
      db(k+1) = (wx2*dpsizr(iix2')-w0*dpsizr(ii0'))*10;
      k  = k+1;
    end
    
    % Modify coil currents to correct flux errors at boundary points
    % Do it in such a way that dcc is changed as little as possible
    % This is a compromise between preserving boundary and cc

    dcc(icx) = -pinv([dbdcc(:,icx); eye(ncx); bus_code(icx)])*[db; zeros(ncx,1); bus_code_lumped*cc];
    % This could use some improvement regarding the bus constraint
    
    % Calculate the perturbation that controls boundary and cc
    C = [efmpc(iused,:)*dcc; 0; 0];
    pert2 = Minv*C(ii);
    psizr(iused) = psizr(iused)+reshape(pert2(1:nused),siused(1),siused(2));
    psizr(inotused) = psizr_app(inotused)+psizr_pla(inotused);
    
    if     options.iconstraints == 0
      cp = 0;
      cf = 0;
    elseif options.iconstraints == 1
      cp = 0;
      cf = pert1(end)+pert2(end);
    elseif options.iconstraints == 2
      cp = pert1(end)+pert2(end);
      cf = 0;
    elseif options.iconstraints == 3
      cp = pert1(end-1)+pert2(end-1);
      cf = pert1(end)+pert2(end);
    end
        
    % Now we have a dpsizr that corrects convergence error and W, Ip error
    % while preserving rbbbs,zbbbs, cc in a good way
    ;
    dpsizr(iused) = reshape(pert1(1:nused),siused(1),siused(2)) + ...
                    reshape(pert2(1:nused),siused(1),siused(2));
    
    converrors(end+1) = max(abs(dpsizr(iused)));
    % If we have converged then do options.nxiter extra iterations
    if min(converrors) < options.converrmax | xiter > -1
      xiter = xiter+1;
    end

    if iterations > 2 & (converrors(end)/converrors(end-1) > 2 | converrors(end) > 10*min(converrors))
      convFailure = 1;
    end

    % If this is the min(converrors) that means it requires less correction than any previous psizr
    % We therefore archive this as the solution
    if min(converrors) == converrors(end)
      cpasma = 1e6*dr*dz*sum(jphi(:));
      
      % Pack converged equilibrium
      ceq.pcurrt = 1e6*dr*dz*jphi;
      ceq.jphi = jphi;
      ceq.psizr = psizr_app+psizr_pla; % psizr only contains data for points iused
      ceq.pprime = pprime;
      ceq.ffprim = ffprim;
      ceq.pres = pres;
      ceq.fpol = fpol;
      ceq.bzero = fpol(end)/ceq.rzero;
      ceq.cc = cc+dcc;
      ceq.psimag = psimag;
      ceq.psibry = psibry;
      ceq.psibnd = psibry;
      ceq.cpasma = cpasma;
      ceq.rmaxis = rmaxis;
      ceq.zmaxis = zmaxis;
      ceq.rg = rg;
      ceq.zg = zg;
      ceq.dr = dr;
      ceq.dz = dz;
      ceq.qpsi = qpsi;
      ceq.nw = nr;
      ceq.nh = nz;
      ceq.psirz = -ceq.psizr'/2/pi;
      ceq.ssibry = -psibry/2/pi;
      ceq.ssimag = -psimag/2/pi;
      
      % Create structure: extra information about equilibrium  
      psizr_r = [zeros(nz,1) ceq.psizr(:,3:end)-ceq.psizr(:,1:end-2) zeros(nz,1)]/2/dr;
      psizr_z = [zeros(1,nr);ceq.psizr(3:end,:)-ceq.psizr(1:end-2,:);zeros(1,nr)]/2/dz;
      Bp2zr = (psizr_r.^2+psizr_z.^2)/twopi^2./rgg.^2;
      eqx.V = sum(2*pi*rgg(iplasma)*dr*dz);
      DX = diff(rbbbs); DY = diff(zbbbs);
      cpieces = sqrt(DX.^2+DY.^2);
      Cl = sum(cpieces); % Contour length
      bp2flx = (mu0*cpasma/Cl)^2;
      Bp2V = sum(dr*dz*twopi*rgg(iplasma).*Bp2zr(iplasma));
      eqx.dcc = dcc;
      eqx.drbbbs = rbbbs - eq.rbbbs(1:nbbbs);
      eqx.dzbbbs = zbbbs - eq.zbbbs(1:nbbbs);
      eqx.betap = 4/3*mu0*Wth/eqx.V/bp2flx; % 2*mu0 * Volume-averaged pressure / bp2flx
      eqx.li = Bp2V/eqx.V/bp2flx; % Volume-averaged Bp^2 / bp2flx
      eqx.Wth = Wth;
      eqx.converrors = converrors;
      eqx.Iff = sum(FFprim(iplasma)/mu0./rgg(iplasma))*dr*dz;
      eqx.Ipp = sum(rgg(iplasma).*Pprime(iplasma))*dr*dz;
      eqx.convFailure = convFailure;
      eqx.ilimited = ilimited;
      eqx.rx1 = rx1;
      eqx.zx1 = zx1;
      eqx.rx2 = rx2;
      eqx.zx2 = zx2;
      eqx.r0 = r0;
      eqx.z0 = z0;
      eqx.dpsizr = dpsizr;

    end    
    
    dpsimag = dpsizr(iia)*wa'; dpsibry = dpsizr(ii0)*w0';
    dpres = cp*pres;
    dpprime = cp*pprime-pprime*(dpsimag-dpsibry)/(psimag-psibry);
    dffprim = cf*ffprim-ffprim*(dpsimag-dpsibry)/(psimag-psibry);
    
    % psizr was updated earlier. Now update the rest
    pres = pres+dpres;
    pprime = pprime+dpprime;
    ffprim = ffprim+dffprim;
    cc = cc+dcc;
    
    % Plot errors
    if idoplot
      clf
      h=subplot(2,2,2); set(h,'pos',[0.60 0.58 0.37 0.34]);
      plot(dpsizr(iused))
      title(['\delta\Psi, no: ' num2str(iterations)],'FontSize',15)
      set(gca,'FontSize',15)
      set(gca,'fontw','b')
      a = axis; a(2)=length(dpsizr(iused)); axis(a);
      h=subplot(2,2,4); set(h,'pos',[0.60 0.09 0.37 0.34]);
      plot(dcc,'r')
      title(['dcc'],'FontSize',15)
      set(gca,'FontSize',15)
      set(gca,'fontw','b')
      h=subplot(2,2,1); set(h,'pos',[0.1 0.08 0.4 0.84]);
      hold on
      plot(rbbbs,zbbbs,'g',rl,zl,'m','linew',3)
      contour(rg,zg,psizr,psigrid)
      set(gca,'FontSize',15)
      set(gca,'fontw','b')
      axis('image')
      axis([rg(1) rg(end) zg(1) zg(end)])
      inew = setdiff(iplasma,iplasma0);
      for k = 1:length(inew)
        rectangle('Position',[rgg(inew(k))-dr/2,zgg(inew(k))-dz/2,dr,dz],'facecolor','red')
      end
      iold = setdiff(iplasma0,iplasma);
      for k = 1:length(iold)
        rectangle('Position',[rgg(iold(k))-dr/2,zgg(iold(k))-dz/2,dr,dz],'facecolor','blue')
      end
      if idoplot < 0
        dum = input('Enter new idoplot value or just press return for next plot. ');
	if ~isempty(dum), idoplot = dum; end
      else
        pause(idoplot)
      end
    end
  end % End of iterations

  
  % Trace flux contours, get rbbbs, zbbbs, qpsi
  % Restore best converged equilibrium (not always the last one) to update rbbbs, zbbbs
  struct_to_ws(ceq);
  % If cccirc was specified then put currents from circuits into original coil vector
  if isfield(options,'cccirc')
    cc = eq.cc;
    for j = 1:length(cccirc)
      cc(j) = sign(cccirc(j))*ceq.cc(abs(cccirc(j)));
    end
    ceq.cc = cc;
    eqx.dcc = ceq.cc-eq.cc;
  end
  r0 = eqx.r0; z0 = eqx.z0;
  psigrid = linspace(psimag,psibry,nw)'; % Wb
  rbbbs = rbbbs(1:nbbbs);
  zbbbs = zbbbs(1:nbbbs);
  rs = zeros(na,nw); zs = zeros(na,nw);
  for ia = 1:na
    th = angle(rbbbs(ia)-rmaxis+sqrt(-1)*(zbbbs(ia)-zmaxis));
    ths(ia) = th;
    % Find initial rho
    rho = sqrt((rbbbs(ia)-rmaxis)^2+(zbbbs(ia)-zmaxis)^2)-dr;
    R = rmaxis+rho*cos(th); Z = zmaxis+rho*sin(th);
    kr0 = min(nr-3,max(1,floor((R-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
    kz1 = min(nz-2,max(2,ceil((Z-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
    k = kr0*nz+kz1;
    iip = k+neighbors;
    tr = (R-rgg(k))/dr; tz = (Z-zgg(k))/dz;
    wp = reshape(mx'*[1 tz tz^2 tz^3]'*[1 tr tr^2 tr^3]*mx,1,16);
    psiprev = psizr(iip)*wp';
    for j = nw:-1:1
      % Fine-tune the rho
      drho = -dr/5; l = 0;
      while l < 15
	rho = rho+drho;
	R = rmaxis+rho*cos(th); Z = zmaxis+rho*sin(th);
        kr0 = min(nr-3,max(1,floor((R-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
        kz1 = min(nz-2,max(2,ceil((Z-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
        k = kr0*nz+kz1;
        iip = k+neighbors;
	tr = (R-rgg(k))/dr; tz = (Z-zgg(k))/dz;
	wp = reshape(mx'*[1 tz tz^2 tz^3]'*[1 tr tr^2 tr^3]*mx,1,16);
	psi = psizr(iip)*wp';
	if abs(psi-psigrid(j)) > abs(psiprev-psigrid(j)) % Switch direction
	  drho = -drho/2; l = l+1;
	end
	if j > nw/2 & psi/(psibry-psimag) > psigrid(j)/(psibry-psimag)
	  drho = -abs(drho);
	end
	psiprev = psi;
      end % end of l loop
      wpr = reshape(mx'*[1 tz tz^2 tz^3]'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
      wpz = reshape(mx'*[0 1 2*tz 3*tz^2]'/dz*[1 tr tr^2 tr^3]*mx,1,16);
      rs(ia,j) = R;
      zs(ia,j) = Z;
      br(ia,j) = -psizr(iip)*wpz'/R;
      bz(ia,j) = +psizr(iip)*wpr'/R;
    end
  end % End of ia loop

  rs(ibdef,end) = r0;
  zs(ibdef,end) = z0;
  if ~ilimited
    br(ibdef,end) = 0;
    bz(ibdef,end) = 0;
  end
  rbbbs(1:na) = rs(1:na,end);
  zbbbs(1:na) = zs(1:na,end);
  rbbbs(end) = rbbbs(1);
  zbbbs(end) = zbbbs(1);
  ceq.rbbbs = rbbbs;
  ceq.zbbbs = zbbbs;
  % Now calculate qpsi
  qpsi = zeros(nw,1);
  for j = 2:nw
    for ia1 = 1:na
      ia2 = ia1+1-na*(ia1==na);
      bp1 = sqrt(br(ia1,j)^2+bz(ia1,j)^2);
      bp2 = sqrt(br(ia2,j)^2+bz(ia2,j)^2);
      d = sqrt((rs(ia1,j)-rs(ia2,j))^2+(zs(ia1,j)-zs(ia2,j))^2);
      if bp1 == 0, bp1=5e-2; end
      if bp2 == 0, bp2=5e-2; end
      if bp1 ~= 0 & bp2 ~= 0
        qpsi(j) = qpsi(j) - d*fpol(j)*(1/bp1/rs(ia1,j)^2+1/bp2/rs(ia2,j)^2)/2;
      elseif bp1 == 0 % Model Bp*R/Bphi = t = k*y^e close to x-point, where y is distance to x-point
        ia3 = ia2+1-na*(ia2==na);
	d2 = sqrt((rs(ia3,j)-rs(ia2,j))^2+(zs(ia3,j)-zs(ia2,j))^2);
        bp3 = sqrt(br(ia3,j)^2+bz(ia3,j)^2);
	t2 = bp2*rs(ia2,j)^2/fpol(j);
	t3 = bp3*rs(ia3,j)^2/fpol(j);
	e = min(0.9,(log(t3)-log(t2))/(log(d+d2)-log(d)));
	k = t2/d^e;
	qpsi(j) = qpsi(j) - d^(1-e)/k/(1-e);
      elseif bp2 == 0
        ia0 = ia1-1+na*(ia1==1);
	d2 = sqrt((rs(ia0,j)-rs(ia1,j))^2+(zs(ia0,j)-zs(ia1,j))^2);
        bp0 = sqrt(br(ia0,j)^2+bz(ia0,j)^2);
	t1 = bp1*rs(ia1,j)^2/fpol(j);
	t0 = bp0*rs(ia0,j)^2/fpol(j);
	e = min(0.9,(log(t0)-log(t1))/(log(d+d2)-log(d)));
	k = t1/d^e;
	qpsi(j) = qpsi(j) - d^(1-e)/k/(1-e);
      end
    end
  end
  p = polyfit(2:5,qpsi(2:5)',2); qpsi(1) = polyval(p,1);
  ceq.qpsi = qpsi;
  
  eqx.drbbbs = rbbbs - eq.rbbbs(1:nbbbs);
  eqx.dzbbbs = zbbbs - eq.zbbbs(1:nbbbs);
  eqx.rs = rs;
  eqx.zs = zs;
  eqx.br = br;
  eqx.bz = bz;
  eqx.rsurf = (max(rbbbs(1:nbbbs))+min(rbbbs(1:nbbbs)))/2;
  eqx.bt = abs(ceq.fpol(end)/eqx.rsurf);
  eqx.aminor = (max(rbbbs(1:nbbbs))-min(rbbbs(1:nbbbs)))/2;
  eqx.ip_normalized = abs(cpasma/1e6/eqx.aminor/eqx.bt);
  eqx.betat = 100*2/3*eqx.Wth/eqx.V/eqx.bt^2*2*mu0;
  eqx.betan = eqx.betat/eqx.ip_normalized;
  eqx.converrors = converrors;
  eqx.convFailure = convFailure;
  eqx.rcur = jphi(:)'*rgg(:)/sum(jphi(:));
  eqx.zcur = jphi(:)'*zgg(:)/sum(jphi(:));

  eqx.descriptions.V = 'plasma volume';
  eqx.descriptions.dcc = 'the change of cc';
  eqx.descriptions.drbbbs  = 'rbbbs - rbbbs0: radial change of boundary';
  eqx.descriptions.dzbbbs  = 'zbbbs - zbbbs0: vertical change of boundary';
  eqx.descriptions.betap  = '4/3*mu0*Wth/V/bp2flx: betap per efit definition';
  eqx.descriptions.li  = 'Bp2V/V(end)/bp2flx: li per efit definition';
  eqx.descriptions.Wth = 'Total thermal energy';
  eqx.descriptions.converrors  = 'convergence errors for all iterations';
  eqx.descriptions.Iff  = 'plasma current due to ffprim';
  eqx.descriptions.Ipp  = 'plasma current due to pprime';
  eqx.descriptions.convFailure = '1 if convergence failed, 0 otherwise';
  eqx.descriptions.ilimited = '1 if plasma touches limiter, 0 otherwise';
  eqx.descriptions.rx1 = 'R of upper null';
  eqx.descriptions.zx1 = 'Z of upper null';
  eqx.descriptions.rx2 = 'R of lower null';
  eqx.descriptions.zx2 = 'Z of lower null';
  eqx.descriptions.r0 = 'R of point that decides boundary';
  eqx.descriptions.z0 = 'Z of point that decides boundary';
  eqx.descriptions.dpsizr = 'Remaining error: psizr-psizr_pla-psizr_app';
  eqx.descriptions.rs = 'Radius of nbbbs-1 points on nw flux surfaces from axis to boundary';
  eqx.descriptions.zs = 'Z position of nbbbs-1 points on nw flux surfaces from axis to boundary';
  eqx.descriptions.br = 'radial magnetic field at points rs,zs';
  eqx.descriptions.bz = 'vertical magnetic field at points rs,zs';
  eqx.descriptions.rsurf = 'average of max and min radius of boundary';
  eqx.descriptions.aminor = 'half the difference of max and min radius of boundary';
  eqx.descriptions.bt = 'abs(fpol(end)/rsurf)';
  eqx.descriptions.ip_normalized = 'abs(cpasma/1e6/aminor/bt) = (plasma current in MA)/aminor/bt';
  eqx.descriptions.betat = '100*2/3*Wth/V/bt^2*2*mu0 = 2*mu0*(volume-averaged pressure)/bt^2 * 100';
  eqx.descriptions.betan = 'betat/ip_normalized';
  eqx.descriptions.rcur = 'R of current centroid';
  eqx.descriptions.zcur = 'Z of current centroid';

  if idoplot
    clf
    h = subplot(2,2,2); set(h,'pos',[0.60 0.58 0.37 0.34]);
    plot(eqx.dpsizr(iused))
    title(['\delta\Psi, no: ' num2str(iterations)],'FontSize',15)
    set(gca,'FontSize',15)
    set(gca,'fontw','b')
    a = axis; a(2)=length(dpsizr(iused)); axis(a);
    h = subplot(2,2,4); set(h,'pos',[0.60 0.09 0.37 0.34]);
    plot(eqx.dcc,'r')
    title(['dcc'],'FontSize',15)
    set(gca,'FontSize',15)
    set(gca,'fontw','b')
    h = subplot(2,2,1); set(h,'pos',[0.1 0.08 0.4 0.84]);
    hold on
    plot(rbbbs,zbbbs,'g',rl,zl,'m','linew',3)
    contour(rg,zg,psizr,psigrid)
    title(['Final result with new boundary'],'FontSize',15)
    set(gca,'FontSize',15)
    set(gca,'fontw','b')
    axis('image')
    axis([rg(1) rg(end) zg(1) zg(end)])
    if idoplot < 0
      wait
    else
      pause(idoplot)
    end
    hold off
  end
