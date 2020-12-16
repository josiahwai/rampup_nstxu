%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gs_find_design_dx
%
%  PURPOSE: Find the change of x needed to approach targets
%
%  INPUTS:  targets, weights, limits, locks - gsdesign variables
%
%  OUTPUTS: dx, the change of x (=[ic;iv;sp;sf;er])
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander  ON	3/5/14
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ev = []; % Will be the error vector containing what we have minus what we want
devdx = []; % How the error vector responds to x = [ic; iv; sp; sf; er]
strev = ''; % Will hold description of each target

iev = 0; % Index of item in error vector

outofrange = false;

% Separatrix
if isfield(targets,'rsep') & isfield(targets,'zsep')
  for i = 1:min(length(targets.rsep),length(targets.zsep))
    if isfield(weights,'sep') & length(weights.sep) >= i
      weight = weights.sep(i);
    else
      weight = 1;
    end
    r = targets.rsep(i);
    z = targets.zsep(i);
    if ~isnan(r) & ~isnan(z) & ~isnan(weight)
      ir = (r-rg(1))/dr+1;
      iz = (z-zg(1))/dz+1;
      kr = floor(ir);
      kz = floor(iz);
      tr = ir-kr;
      tz = iz-kz;
      if kr < 2 | kr > nr-2 | kz < 2 | kz > nz-2
	outofrange = true;
	if iteration_counter == 0
          disp(['warning in gsdesign: target rsep = ' num2str(r), ...
            ', zsep = ' num2str(z) ' is out of range, and will be ignored'])
	end
      else
	iev = iev+1;
	wsep = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
	iisep = kz+nz*(kr-1)+neighbors';
	psisep = wsep*psizr(iisep);
	dpsisepdx = wsep*dpsizrdx(iisep,:);
	ev(iev,1) = weight*(psisep - psibry);
	devdx(iev,:) = weight*(dpsisepdx - dpsibrydx);
	str = ['sep' num2str(i)];
	strev(iev,1:length(str)) = str;
      end
    end  
  end
end

% x-points
if isfield(targets,'rx') & isfield(targets,'zx')
  for i = 1:min(length(targets.rx),length(targets.zx))
    if isfield(weights,'x') & numel(weights.x) >= i
      weight = weights.x(i);
    else
      weight = 1;
    end
    r = targets.rx(i);
    z = targets.zx(i);
    if ~isnan(r) & ~isnan(z) & ~isnan(weight)
      ir = (r-rg(1))/dr+1;
      iz = (z-zg(1))/dz+1;
      kr = floor(ir);
      kz = floor(iz);
      tr = ir-kr;
      tz = iz-kz;
      if kr < 2 | kr > nr-2 | kz < 2 | kz > nz-2
	outofrange = true;
	if iteration_counter == 0
          disp(['warning in gsdesign: target rx = ' num2str(r), ...
            ', zx = ' num2str(z) ' is out of range, and will be ignored'])
	end
      else
	iix = kz+nz*(kr-1)+neighbors';
	wr0 = [1 tr tr^2 tr^3]*mx;
	wz0 = mx'*[1 tz tz^2 tz^3]';
	wr1 = [0 1 2*tr 3*tr^2]*mx/dr;
	wz1 = mx'*[0 1 2*tz 3*tz^2]'/dz;
	wxr = reshape(wz0*wr1, 1, 16);
	wxz = reshape(wz1*wr0, 1, 16);
	psixr = wxr*psizr(iix);
	psixz = wxz*psizr(iix);
	dpsixrdx = wxr*dpsizrdx(iix,:);
	dpsixzdx = wxz*dpsizrdx(iix,:);
	iev = iev+1;
	ev(iev,1) = weight*psixr;
	devdx(iev,:) = weight*dpsixrdx;
	str = ['\partial\psi_x_' num2str(i) '/\partialR'];
	strev(iev,1:length(str)) = str;
	iev = iev+1;
	ev(iev,1) = weight*psixz;
	devdx(iev,:) = weight*dpsixzdx;
	str = ['\partial\psi_x_' num2str(i) '/\partialZ'];
	strev(iev,1:length(str)) = str;
      end
    end
  end
end

% Snowflakes
if isfield(targets,'rsnf') & isfield(targets,'zsnf')
  for i = 1:min(length(targets.rsnf),length(targets.zsnf))
    if isfield(weights,'snf') & i <= length(weights.snf)
      weight = weights.snf(i);
    else
      weight = 1;
    end
    r = targets.rsnf(i);
    z = targets.zsnf(i);
    if ~isnan(r) & ~isnan(z) & ~isnan(weight)
      ir = (r-rg(1))/dr+1;
      iz = (z-zg(1))/dz+1;
      kr = floor(ir);
      kz = floor(iz);
      tr = ir-kr;
      tz = iz-kz;
      if kr < 2 | kr > nr-2 | kz < 2 | kz > nz-2
	outofrange = true;
	if iteration_counter == 0
          disp(['warning in gsdesign: target rsnf = ' num2str(r), ...
            ', zsnf = ' num2str(z) ' is out of range, and will be ignored'])
	end
      else
	iisnf = kz+nz*(kr-1)+neighbors';
	wr0 = [1 tr tr^2 tr^3]*mx;
	wz0 = mx'*[1 tz tz^2 tz^3]';
	wr1 = [0 1 2*tr 3*tr^2]*mx/dr;
	wz1 = mx'*[0 1 2*tz 3*tz^2]'/dz;
	wr2 = [0 0 2 6*tr]*mx/dr^2;
	wz2 = mx'*[0 0 2 6*tz]'/dz^2;
	wsnf = reshape(wz0*wr0, 1, 16);
	wsnfr = reshape(wz0*wr1, 1, 16);
	wsnfz = reshape(wz1*wr0, 1, 16);
	wsnfrr = reshape(wz0*wr2, 1, 16);
	wsnfrz = reshape(wz1*wr1, 1, 16);
	wsnfzz = reshape(wz2*wr0, 1, 16);
	% tsnf is the angle (theta) toward outer x-point
	if isfield(targets,'tsnf') & i <= length(targets.tsnf)
	  t = targets.tsnf(i); 
	else % Pick an angle aimed straight away from the plasma
          t = 40*sign(zmaxis-z);
	end
	% nsnf is number of rays going out from the snowflake
	if isfield(targets,'nsnf') & i <= length(targets.nsnf)
	  nsnf = targets.nsnf(i);
	else
          nsnf = 6;
	end
	% dt is angle in degrees between snowflake rays
	dt = 360/nsnf;
	% rhosnf is a snowflake "radius"
	if isfield(targets,'rhosnf') & i <= length(targets.rhosnf)
	  rhosnf = targets.rhosnf(i);
	else
          rhosnf = (dr+dz)/2;
	end
	iev = iev+1;
	ev(iev,1) = weight*10*wsnfr*psizr(iisnf);
	devdx(iev,:) = weight*10*wsnfr*dpsizrdx(iisnf,:);
	strev(iev,1:10) = 'snowflakes';
	iev = iev+1;
	ev(iev,1) = weight*10*wsnfz*psizr(iisnf);
	devdx(iev,:) = weight*10*wsnfz*dpsizrdx(iisnf,:);
	strev(iev,1:10) = 'snowflakes';
  if 01 % Control outer x-point
	r = targets.rsnf(i)+rhosnf*cos(t*pi/180)/2;
	z = targets.zsnf(i)+rhosnf*sin(t*pi/180)/2;
	ir = (r-rg(1))/dr+1;
	iz = (z-zg(1))/dz+1;
	kr = max(2,min(nr-2,floor(ir)));
	kz = max(2,min(nz-2,floor(iz)));
	tr = ir-kr;
	tz = iz-kz;
	iip = kz+nz*(kr-1)+neighbors';
	wr0 = [1 tr tr^2 tr^3]*mx;
	wz0 = mx'*[1 tz tz^2 tz^3]';
	wr1 = [0 1 2*tr 3*tr^2]*mx/dr;
	wz1 = mx'*[0 1 2*tz 3*tz^2]'/dz;
	wp = reshape(wz0*wr0, 1, 16);
	wpr = reshape(wz0*wr1, 1, 16);
	wpz = reshape(wz1*wr0, 1, 16);
	iev = iev+1;
	ev(iev,1) = weight/10*wpr*psizr(iip);
	devdx(iev,:) = weight/10*wpr*dpsizrdx(iip,:);
	strev(iev,1:10) = 'snowflakes';
	iev = iev+1;
	ev(iev,1) = weight/10*wpz*psizr(iip);
	devdx(iev,:) = weight/10*wpz*dpsizrdx(iip,:);
	strev(iev,1:10) = 'snowflakes';
  end
	for it = 1:nsnf % Control flux to form the snowflake
	  r = targets.rsnf(i)+rhosnf*cos((t+(it-0.5)*dt)*pi/180);
	  z = targets.zsnf(i)+rhosnf*sin((t+(it-0.5)*dt)*pi/180);
	  ir = (r-rg(1))/dr+1;
	  iz = (z-zg(1))/dz+1;
	  kr = max(2,min(nr-2,floor(ir)));
	  kz = max(2,min(nz-2,floor(iz)));
	  tr = ir-kr;
	  tz = iz-kz;
          iip = kz+nz*(kr-1)+neighbors';
          wp = reshape(mx'*[1 tz tz^2 tz^3]'*[1 tr tr^2 tr^3]*mx, 1, 16);
	  iev = iev+1;
	  ev(iev,1) = weight*(wp*psizr(iip)-wsnf*psizr(iisnf));
	  devdx(iev,:) = weight*(wp*dpsizrdx(iip,:)-wsnf*dpsizrdx(iisnf,:));
	  strev(iev,1:10) = 'snowflakes';
	end
      end
    end
  end
end

bdef.target.exists = false;
if isfield(targets,'rbdef') & isfield(targets,'zbdef') & ...
   ~isempty(targets.rbdef)  & ~isempty(targets.zbdef)
  if min(length(targets.rbdef),length(targets.zbdef)) > 1 & iteration_counter == 1
    disp('warning in gsdesign: Only first target for rbdef, zbdef used')
  end
  if isfield(weights,'bdef') & ~isempty(weights.bdef)
    weight = weights.bdef(1);
  else
    weight = 1;
  end
  r = targets.rbdef(1);
  z = targets.zbdef(1);
  if ~isnan(r) & ~isnan(z) & ~isnan(weight)
    ir = (r-rg(1))/dr+1;
    iz = (z-zg(1))/dz+1;
    kr = floor(ir);
    kz = floor(iz);
    tr = ir-kr;
    tz = iz-kz;
    xA = r+diff(zl')*[0 1];
    yA = z-diff(rl')*[0 1];
    xB = [rl(1:nl-1)' rl(2:nl)'];
    yB = [zl(1:nl-1)' zl(2:nl)'];
    [X, Y, flag, tA, tB] = line_intersection(xA, yA, xB, yB);
    dA = tA.*dl(1:nl-1);
    dlmin = min(abs(dA(flag > 1)));
    if kr < 2 | kr > nr-2 | kz < 2 | kz > nz-2
      outofrange = true;
      if iteration_counter == 0
        disp(['warning in gsdesign: target point rbdef = ' num2str(r), ...
          ', zbdef = ' num2str(z) ' is out of range and will be ignored'])
      end
    elseif dlmin > min(dr,dz)/100 & ~isinpoly(r,z)
      if iteration_counter == 0
        disp(['warning in gsdesign: target rbdef = ' num2str(r), ...
          ', zbdef = ' num2str(z) ' is outside limiter, and will be ignored'])    
      end
    else
      iev = iev+1;
      wsep = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      iisep = kz+nz*(kr-1)+neighbors';
      psisep = wsep*psizr(iisep);
      dpsisepdx = wsep*dpsizrdx(iisep,:);
      ev(iev,1) = weight*(psisep - psibry);
      devdx(iev,:) = weight*(dpsisepdx - dpsibrydx);
      str = '\psi_b_d_e_f';
      strev(iev,1:length(str)) = str;
      % Bad idea to include x points as part of assignment! Hence if 0
      if 0 & dlmin > min(dr,dz)/100 % Need a null at this point
 	wr0 = [1 tr tr^2 tr^3]*mx;
	wz0 = mx'*[1 tz tz^2 tz^3]';
	wr1 = [0 1 2*tr 3*tr^2]*mx/dr;
	wz1 = mx'*[0 1 2*tz 3*tz^2]'/dz;
	wxr = reshape(wz0*wr1, 1, 16);
	wxz = reshape(wz1*wr0, 1, 16);
	psixr = wxr*psizr(iisep);
	psixz = wxz*psizr(iisep);
	dpsixrdx = wxr*dpsizrdx(iisep,:);
	dpsixzdx = wxz*dpsizrdx(iisep,:);
	iev = iev+1;
	ev(iev,1) = weight*psixr;
	devdx(iev,:) = weight*dpsixrdx;
	str = ['\partial\psi_b_d_e_f_' num2str(i) '/\partialR'];
	strev(iev,1:length(str)) = str;
	iev = iev+1;
	ev(iev,1) = weight*psixz;
	devdx(iev,:) = weight*dpsixzdx;
	str = ['\partial\psi_b_d_e_f_' num2str(i) '/\partialZ'];
	strev(iev,1:length(str)) = str;       
      end
      % Archive data needed to make flux limits for other parts of boundary
      bdef.target.exists = true; % exists and inside limiter or very close
      bdef.target.dlmin = dlmin; % minimum distance to limiter
      bdef.target.r = r;
      bdef.target.z = z;
      bdef.target.w = wsep;
      bdef.target.ii = iisep;
      bdef.target.psi = psisep;
      bdef.target.dpsidx = dpsisepdx;            
    end
  end
end

if isfield(targets,'zcur')
  if isfield(weights,'zcur')
    weight = weights.zcur;
  else
    weight = 1;
  end
  if ~isnan(targets.zcur) & ~isnan(weight)
    iev = iev+1;
    ev(iev,1) = weight*(zcur - targets.zcur);
    devdx(iev,:) = weight*dzcurdx;
    strev(iev,1:7) = 'Z_c_u_r';
  end
end

if isfield(targets,'cpasma')
  if isfield(weights,'cpasma')
    weight = weights.cpasma;
  else
    weight = 1;
  end
  if ~isnan(targets.cpasma) & ~isnan(weight)
    iev = iev+1;
    ev(iev,1) = weight*(cpasma - targets.cpasma);
    devdx(iev,:) = weight*dcpasmadx;
    strev(iev,1:3) = 'I_p';
  end
end

if isfield(targets,'li')
  if isfield(weights,'li')
    weight = weights.li;
  else
    weight = 1;
  end
  if ~isnan(targets.li) & ~isnan(weight)
    iev = iev+1;
    ev(iev,1) = weight*(li - targets.li);
    devdx(iev,:) = weight*dlidx;
    strev(iev,1:3) = 'l_i';
  end
end

if isfield(targets,'betap')
  if isfield(weights,'betap')
    weight = weights.betap;
  else
    weight = 1;
  end
  if ~isnan(targets.betap) & ~isnan(weight)
    iev = iev+1;
    ev(iev,1) = weight*(betap - targets.betap);
    devdx(iev,:) = weight*dbetapdx;
    strev(iev,1:7) = '\beta_p';
  end
end

if isfield(targets,'betan')
  if isfield(weights,'betan')
    weight = weights.betan;
  else
    weight = 1;
  end
  if ~isnan(targets.betan) & ~isnan(weight)
    iev = iev+1;
    ev(iev,1) = weight*(betan - targets.betan);
    devdx(iev,:) = weight*dbetandx;
    strev(iev,1:7) = '\beta_n';
  end
end

if isfield(targets,'psibry')
  if isfield(weights,'psibry')
    weight = weights.psibry;
  else
    weight = 1;
  end
  if ~isnan(targets.psibry) & ~isnan(weight)
    iev = iev+1;
    ev(iev,1) = weight*(psibry - targets.psibry);
    devdx(iev,:) = weight*dpsibrydx;
    strev(iev,1:6) = '\psi_b';
  end
end

if isfield(targets,'psimag')
  if isfield(weights,'psimag')
    weight = weights.psimag;
  else
    weight = 1;
  end
  if ~isnan(targets.psimag) & ~isnan(weight)
    iev = iev+1;
    ev(iev,1) = weight*(psimag - targets.psimag);
    devdx(iev,:) = weight*dpsimagdx;
    strev(iev,1:6) = '\psi_a';
  end
end

if isfield(targets,'psipla')
  if isfield(weights,'psipla')
    weight = weights.psipla;
  else
    weight = 1;
  end
  if ~isnan(targets.psipla) & ~isnan(weight)
    iev = iev+1;
    ev(iev,1) = weight*(psipla - targets.psipla);
    devdx(iev,:) = weight*dpsipladx;
    strev(iev,1:6) = '\psi_p';
  end
end

if isfield(targets,'fl') & isfield(index_in_y,'fl')
  for i = 1:min(length(targets.fl),length(index_in_y.fl))
    if isfield(weights,'fl') & length(weights.fl) >= i
      weight = weights.fl(i);
    else
      weight = 1;
    end
    if ~isnan(targets.fl(i)) & ~isnan(weight)
      iev = iev+1;
      ev(iev,1) = weight*(y(index_in_y.fl(i)) - targets.fl(i));
      devdx(iev,:) = weight*dydx(index_in_y.fl(i),:);
      str = ['fl' num2str(i)];
      strev(iev,1:length(str)) = str;
    end
  end
end

if isfield(targets,'bp') & isfield(index_in_y,'bp')
  for i = 1:min(length(targets.bp),length(index_in_y.bp))
    if isfield(weights,'bp') & length(weights.bp) >= i
      weight = weights.bp(i);
    else
      weight = 1;
    end
    if ~isnan(targets.bp(i)) & ~isnan(weight)
      iev = iev+1;
      ev(iev,1) = weight*(y(index_in_y.bp(i)) - targets.bp(i));
      devdx(iev,:) = weight*dydx(index_in_y.bp(i),:);
      str = ['bp' num2str(i)];
      strev(iev,1:length(str)) = str;
    end
  end
end

if isfield(targets,'rog') & isfield(index_in_y,'rog')
  for i = 1:min(length(targets.rog),length(index_in_y.rog))
    if isfield(weights,'rog') & length(weights.rog) >= i
      weight = weights.rog(i);
    else
      weight = 1;
    end
    if ~isnan(targets.rog(i)) & ~isnan(weight)
      iev = iev+1;
      ev(iev,1) = weight*(y(index_in_y.rog(i)) - targets.rog(i));
      devdx(iev,:) = weight*dydx(index_in_y.rog(i),:);
      str = ['rog' num2str(i)];
      strev(iev,1:length(str)) = str;
    end
  end
end

if isfield(targets,'ic')
  for i = 1:min(length(targets.ic), nc)
    if isfield(weights,'ic') & length(weights.ic) >= i
      weight = weights.ic(i);
    else
      weight = 1;
    end
    if ~isnan(targets.ic(i)) & ~isnan(weight)
      iev = iev+1;
      ev(iev,1) = weight*(ic(i) - targets.ic(i));
      devdx(iev,indic(i)) = weight;
      str = ['ic' num2str(i)];
      strev(iev,1:length(str)) = str;
    end
  end
end

if isfield(targets,'iv')
  for i = 1:min(length(targets.iv), nv)
    if isfield(weights,'iv') & length(weights.iv) >= i
      weight = weights.iv(i);
    else
      weight = 1;
    end
    if ~isnan(targets.iv(i)) & ~isnan(weight)
      iev = iev+1;
      ev(iev,1) = weight*(iv(i) - targets.iv(i));
      devdx(iev,indiv(i)) = weight;
      str = ['iv' num2str(i)];
      strev(iev,1:length(str)) = str;
    end
  end
end

if isfield(targets,'sp')
  for i = 1:min(length(targets.sp), nkn+2)
    if isfield(weights,'sp') & length(weights.sp) >= i
      weight = weights.sp(i);
    else
      weight = 1;
    end
    if ~isnan(targets.sp(i)) & ~isnan(weight)
      iev = iev+1;
      ev(iev,1) = weight*(sp(i) - targets.sp(i));
      devdx(iev,indsp(i)) = weight;
      str = ['sp' num2str(i)];
      strev(iev,1:length(str)) = str;
    end
  end
end

if isfield(targets,'sf')
  for i = 1:min(length(targets.sf), nkn+2)
    if isfield(weights,'sf') & length(weights.sf) >= i
      weight = weights.sf(i);
    else
      weight = 1;
    end
    if ~isnan(targets.sf(i)) & ~isnan(weight)
      iev = iev+1;
      ev(iev,1) = weight*(sf(i) - targets.sf(i));
      devdx(iev,indsf(i)) = weight;
      str = ['sf' num2str(i)];
      strev(iev,1:length(str)) = str;
    end
  end
end

if isfield(targets,'psic')
  if ~isfield(targets,'dpsicdic')
    targets.dpsicdic = 0*targets.psic;
  end
  if ~isfield(targets,'ic0')
    targets.ic0 = 0*targets.psic;
  end
  psic = mcc*ic + mcv*iv + mpc'*pcurrt(:);
  dpsicdx = [mcc mcv zeros(nc,nx-nc-nv)] + mpc'*dpcurrtdx;
  for i = 1:min([length(targets.psic), ...
                 length(targets.dpsicdic), ...
                 length(targets.ic0), nc])
    if isfield(weights,'psic') & length(weights.psic) >= i
      weight = weights.psic(i);
    else
      weight = 1;
    end
    target = targets.psic(i) + targets.dpsicdic(i)*(ic(i)-targets.ic0(i));
    if ~isnan(target) & ~isnan(weight)
      iev = iev+1;
      ev(iev,1) = weight*(psic(i) - target);
      devdx(iev,:) = weight*dpsicdx(i,:);
      devdx(iev,i) = devdx(iev,i)-weight*targets.dpsicdic(i);
      str = ['\psic' num2str(i)];
      strev(iev,1:length(str)) = str;
    end
  end
end

if isfield(targets,'psiv')
  if ~isfield(targets,'dpsivdiv')
    targets.dpsivdiv = 0*targets.psiv;
  end
  targets.dpsivdiv = targets.dpsivdiv(:);
  if numel(targets.dpsivdiv) < nv
    targets.dpsivdiv(nv) = 0;
  end
  if numel(targets.dpsivdiv) > nv
    targets.dpsivdiv = targets.dpsivdiv(1:nv);
  end
  if ~isfield(targets,'iv0')
    targets.iv0 = 0*targets.psiv;
  end
  psiv = mcv'*ic + mvv*iv + mpv'*pcurrt(:);
  dpsivdx = [mcv' mvv-diag(targets.dpsivdiv) zeros(nv,nx-nc-nv)] + mpv'*dpcurrtdx;
  for i = 1:min([length(targets.psiv), ...
                 length(targets.dpsivdiv), ...
                 length(targets.iv0), nv])
    if isfield(weights,'psiv') & length(weights.psiv) >= i
      weight = weights.psiv(i);
    else
      weight = 1;
    end
    target = targets.psiv(i) + targets.dpsivdiv(i)*(iv(i)-targets.iv0(i));
    if ~isnan(target) & ~isnan(weight)
      iev = iev+1;
      ev(iev,1) = weight*(psiv(i) - target);
      devdx(iev,:) = weight*dpsivdx(i,:);
      devdx(iev,nc+i) = devdx(iev,nc+i)-weight*targets.dpsivdiv(i);
      str = ['\psiv' num2str(i)];
      strev(iev,1:length(str)) = str;
    end
  end
end

if isfield(targets,'frc')
  for i = 1:min(length(targets.frc), nc)
    if isfield(weights,'frc') & length(weights.frc) >= i
      weight = weights.frc(i);
    else
      weight = 1;
    end
    if ~isnan(targets.frc(i)) & ~isnan(weight)
      iev = iev+1;
      ev(iev,1) = weight*(frc(i) - targets.frc(i));
      devdx(iev,:) = weight*dfrcdx(i,:);
      str = ['frc' num2str(i)];
      strev(iev,1:length(str)) = str;
    end
  end
end

if isfield(targets,'fzc')
  for i = 1:min(length(targets.fzc), nc)
    if isfield(weights,'fzc') & length(weights.fzc) >= i
      weight = weights.fzc(i);
    else
      weight = 1;
    end
    if ~isnan(targets.fzc(i)) & ~isnan(weight)
      iev = iev+1;
      ev(iev,1) = weight*(fzc(i) - targets.fzc(i));
      devdx(iev,:) = weight*dfzcdx(i,:);
      str = ['fzc' num2str(i)];
      strev(iev,1:length(str)) = str;
    end
  end
end

if isfield(targets,'fluxexp') & ~isnan(targets.fluxexp)
  if isfield(targets,'rfluxexp') & isfield(targets,'zfluxexp')
    rfluxexp = targets.rfluxexp;
    zfluxexp = targets.zfluxexp;
    if isfield(targets,'dfluxexp')
      dfluxexp = targets.dfluxexp;
    else
      dfluxexp = 0;
    end
    if ~isnan(rfluxexp) & ~isnan(zfluxexp) & ~isnan(dfluxexp)
      eqfluxexp.rg = rg;
      eqfluxexp.zg = zg;
      eqfluxexp.rbbbs = rbbbs;
      eqfluxexp.zbbbs = zbbbs;
      eqfluxexp.nbbbs = nbbbs;
      eqfluxexp.psizr = psizr;
      eqfluxexp.psibry = psibry;
      [fluxexp, rfluxexp2, zfluxexp2] = calc_fluxexp(...
	eqfluxexp, rfluxexp, zfluxexp, dfluxexp);
      dfluxexpdpsizr = zeros(1,ngg);
      kz = floor((zbbbso-zg(1))/dz+1);
      kr = floor((rbbbso-rg(1))/dr+1);
      dfluxexpdpsizr(kz+nz*(kr-1)+neighbors) = nan;
      kz = floor((zfluxexp-zg(1))/dz+1);
      kr = floor((rfluxexp-rg(1))/dr+1);
      dfluxexpdpsizr(kz+nz*(kr-1)+neighbors) = nan;
      if dfluxexp > 0
      kz = floor((zfluxexp2-zg(1))/dz+1);
      kr = floor((rfluxexp2-rg(1))/dr+1);
      dfluxexpdpsizr(kz+nz*(kr-1)+neighbors) = nan;
      end
      dpsiba = psibry-psimag;
      for i = 1:ngg
	if isnan(dfluxexpdpsizr(i))
          eqfluxexp.psizr(i) = psizr(i)+dpsiba/1e6;
          dum = calc_fluxexp(eqfluxexp, rfluxexp, zfluxexp, dfluxexp);
          dfluxexpdpsizr(i) = (dum-fluxexp)/dpsiba*1e6;
	  eqfluxexp.psizr(i) = psizr(i);              
	end
      end
      dfluxexpdx = dfluxexpdpsizr*dpsizrdx;
      if isfield(weights,'fluxexp')
	weight = weights.fluxexp;
      else
	weight = 1;
      end
      iev = iev+1;
      ev(iev,1) = weight*(fluxexp - targets.fluxexp);
      devdx(iev,:) = weight*dfluxexpdx;
      str = ['\psiexp'];
      strev(iev,1:length(str)) = str;
    end
  else
    if iteration_counter == 0
      disp(['warning in gsdesign: targets.rfluxexp, targets.zfluxexp ', ...
        ' for calculating fluxexp are missing'])
    end
  end
end

if isfield(targets,'fluxerror')
  if isfield(weights,'fluxerror')
    weight = weights.fluxerror;
  else
    weight = 1;
  end
  if ~isnan(targets.fluxerror) & ~isnan(weight)
    iev = iev+1;
    ev(iev,1) = weight*(fluxerror - targets.fluxerror);
    devdx(iev,end) = -weight;
    strev(iev,1:10) = '\psi_e_r_r';
  end
end

strev(strev==0) = 32;

ilocked = logical(zeros(nxcirc,1));

% dxcirc will hold changes to be made
dxcirc = zeros(nxcirc,1);

% Changes set by locks on currents in circuits
for i = 1:nci
  if ~isnan(locks.ci(i))
    dxcirc(i) = locks.ci(i) - ci(i);
    ilocked(i) = true;  
  end
end

% Changes set by locks on currents in vessel elements
for i = 1:nv
  if ~isnan(locks.iv(i))
    dxcirc(i+nci) = locks.iv(i) - iv(i);
    ilocked(i+nci) = true;  
  end
end

% Change set by flux error correction
dxcirc(nxcirc) = -1;
ilocked(nxcirc) = 1;

if isempty(ev)
  iev = iev+1;
  ev(iev,1) = fluxerror;
  devdx(iev,nx) = -1;
  strev(iev,1:10) = '\psi_e_r_r';  
end

iopen = ~ilocked;
dxlocked = dxdxcirc(:,ilocked)*dxcirc(ilocked);
evc = ev + devdx*dxlocked;

n = sum(iopen);                   % Number of primal variables
devdpv = devdx*dxdxcirc(:,iopen); % d(error vector)/d(primal variable)
s = 1./sqrt(sum(devdpv.^2));
H = 2*devdpv'*devdpv;             % Hessian for devdpv
if ~any(H(:))
  H = eye(n);
end
h = (2*evc'*devdpv)';             % Hessian r.h.s.
Ae = zeros(0,n);                  % Equality constraints matrix
be = [];                          % Equality constraints r.h.s.
Ai = zeros(0,n);                  % Inequility constraints matrix
bi = [];                          % Inequility constraints r.h.s.

if isfield(locks,'rsep') & isfield(locks,'zsep')
  n = min(size(locks.rsep,1),size(locks.zsep,1));
  for i = 1:n
    r = locks.rsep(i);
    z = locks.zsep(i);
    if ~isnan(r) & ~isnan(z) & r > 0
      ir = (r-rg(1))/dr+1;
      iz = (z-zg(1))/dz+1;
      kr = floor(ir);
      kz = floor(iz);
      tr = ir-kr;
      tz = iz-kz;
      if kr < 2 | kr > nr-2 | kz < 2 | kz > nz-2
	outofrange = true;
	if iteration_counter == 0
          disp(['warning in gsdesign: LOCKED rsep = ' num2str(r), ...
            ', zsep = ' num2str(z) ' is out of range, and will be ignored'])
	end
      else
        wsep = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
        iisep = kz+nz*(kr-1)+neighbors';
        psisep = wsep*psizr(iisep);
        dpsisepdx = wsep*dpsizrdx(iisep,:);
        Ae(end+1,:) = (dpsisepdx-dpsibrydx)*dxdxcirc(:,iopen);
        if fluxerror < 0.01
          be(end+1,1) = psibry-dpsibrydx*dxlocked-psisep+dpsisepdx*dxlocked;  
        else
	  be(end+1,1) = psibry-psisep;
	end
      end
    end
  end
end
if isfield(locks,'rx') & isfield(locks,'zx')
  for i = 1:min(length(locks.rx),length(locks.zx))
    r = locks.rx(i);
    z = locks.zx(i);
    if ~isnan(r) & ~isnan(z) & r > 0
      ir = (r-rg(1))/dr+1;
      iz = (z-zg(1))/dz+1;
      kr = floor(ir);
      kz = floor(iz);
      tr = ir-kr;
      tz = iz-kz;
      if kr < 2 | kr > nr-2 | kz < 2 | kz > nz-2
	outofrange = true;
	if iteration_counter == 0
          disp(['warning in gsdesign: LOCKED rx = ' num2str(r), ...
            ', zx = ' num2str(z) ' is out of range, and will be ignored'])
	end
      else
	iix = kz+nz*(kr-1)+neighbors';
 	wr0 = [1 tr tr^2 tr^3]*mx;
	wz0 = mx'*[1 tz tz^2 tz^3]';
	wr1 = [0 1 2*tr 3*tr^2]*mx/dr;
	wz1 = mx'*[0 1 2*tz 3*tz^2]'/dz;
	wxr = reshape(wz0*wr1, 1, 16);
	wxz = reshape(wz1*wr0, 1, 16);
	dpsixrdx = wxr*dpsizrdx(iix,:);
	dpsixzdx = wxz*dpsizrdx(iix,:);
	psixr = wxr*psizr(iix) + dpsixrdx*dxlocked;
	psixz = wxz*psizr(iix) + dpsixzdx*dxlocked;
        Ae(end+1,:) = dpsixrdx*dxdxcirc(:,iopen);
	be(end+1,1) = -psixr;
        Ae(end+1,:) = dpsixzdx*dxdxcirc(:,iopen);
	be(end+1,1) = -psixz;
      end
    end
  end
end
if isfield(locks,'rsnf') & isfield(locks,'zsnf')
  for i = 1:min(length(locks.rsnf),length(locks.zsnf))
    r = locks.rsnf(i);
    z = locks.zsnf(i);
    if ~isnan(r) & ~isnan(z) & r > 0
      ir = (r-rg(1))/dr+1;
      iz = (z-zg(1))/dz+1;
      kr = floor(ir);
      kz = floor(iz);
      tr = ir-kr;
      tz = iz-kz;
      if kr < 2 | kr > nr-2 | kz < 2 | kz > nz-2
	outofrange = true;
	if iteration_counter == 0
          disp(['warning in gsdesign: LOCKED rsnf = ' num2str(r), ...
            ', zsnf = ' num2str(z) ' is out of range, and will be ignored'])
	end
      else
	iisnf = kz+nz*(kr-1)+neighbors';
	wr0 = [1 tr tr^2 tr^3]*mx;
	wz0 = mx'*[1 tz tz^2 tz^3]';
	wr1 = [0 1 2*tr 3*tr^2]*mx/dr;
	wz1 = mx'*[0 1 2*tz 3*tz^2]'/dz;
	wr2 = [0 0 2 6*tr]*mx/dr^2;
	wz2 = mx'*[0 0 2 6*tz]'/dz^2;
	wsnf = reshape(wz0*wr0, 1, 16);
	wsnfr = reshape(wz0*wr1, 1, 16);
	wsnfz = reshape(wz1*wr0, 1, 16);
	wsnfrr = reshape(wz0*wr2, 1, 16);
	wsnfrz = reshape(wz1*wr1, 1, 16);
	wsnfzz = reshape(wz2*wr0, 1, 16);
	dpsisnfrdx = wsnfr*dpsizrdx(iisnf,:);
	dpsisnfzdx = wsnfz*dpsizrdx(iisnf,:);
	psisnfr = wsnfr*psizr(iisnf) + dpsisnfrdx*dxlocked;
	psisnfz = wsnfz*psizr(iisnf) + dpsisnfzdx*dxlocked;
        Ae(end+1,:) = dpsisnfrdx*dxdxcirc(:,iopen);
	be(end+1,1) = -psisnfr;
        Ae(end+1,:) = dpsisnfzdx*dxdxcirc(:,iopen);
	be(end+1,1) = -psisnfz;
	% tsnf is the angle (theta) toward outer x-point
	if isfield(locks,'tsnf') & i <= length(locks.tsnf)
	  t = locks.tsnf(i); 
	else % Pick an angle aimed straight away from the plasma
          t = 40*sign(zmaxis-z);
	end
	% nsnf is number of rays going out from the snowflake
	if isfield(locks,'nsnf') & i <= length(locks.nsnf)
	  nsnf = locks.nsnf(i);
	else
          nsnf = 6;
	end
	% dt is angle in degrees between snowflake rays
	dt = 360/nsnf;
	% rhosnf is a snowflake "radius"
	if isfield(locks,'rhosnf') & i <= length(locks.rhosnf)
	  rhosnf = locks.rhosnf(i);
	else
          rhosnf = (dr+dz);
	end
	for it = 1:nsnf/2 % Control flux to form the snowflake
	  r = locks.rsnf(i)+rhosnf*cos((t+(it-0.5)*dt)*pi/180);
	  z = locks.zsnf(i)+rhosnf*sin((t+(it-0.5)*dt)*pi/180);
	  ir = (r-rg(1))/dr+1;
	  iz = (z-zg(1))/dz+1;
	  kr = max(2,min(nr-2,floor(ir)));
	  kz = max(2,min(nz-2,floor(iz)));
	  tr = ir-kr;
	  tz = iz-kz;
          iip = kz+nz*(kr-1)+neighbors';
          wp = reshape(mx'*[1 tz tz^2 tz^3]'*[1 tr tr^2 tr^3]*mx, 1, 16);
          psip = wp*psizr(iip);
          dpsipdx = wp*dpsizrdx(iip,:);
          Ae(end+1,:) = (dpsipdx-dpsibrydx)*dxdxcirc(:,iopen);
          if fluxerror < 0.01
            be(end+1,1) = psibry-dpsibrydx*dxlocked-psip+dpsipdx*dxlocked;  
          else
	    be(end+1,1) = psibry-psip;
	  end
	end
      end
    end
  end
end
if isfield(locks,'rbdef') & isfield(locks,'zbdef') & ...
  ~isempty(locks.rbdef) & ~isempty(locks.zbdef)
  if min(length(locks.rbdef),length(locks.zbdef)) > 1 & iteration_counter == 1
    disp('warning in gsdesign: Only first value for locked rbdef, zbdef used')
  end
  r = locks.rbdef(1);
  z = locks.zbdef(1);
  if ~isnan(r) & ~isnan(z)
    ir = (r-rg(1))/dr+1;
    iz = (z-zg(1))/dz+1;
    kr = floor(ir);
    kz = floor(iz);
    tr = ir-kr;
    tz = iz-kz;
    xA = r+diff(zl')*[0 1];
    yA = z-diff(rl')*[0 1];
    xB = [rl(1:nl-1)' rl(2:nl)'];
    yB = [zl(1:nl-1)' zl(2:nl)'];
    [X, Y, flag, tA, tB] = line_intersection(xA, yA, xB, yB);
    dA = tA.*dl(1:nl-1);
    dlmin = min(abs(dA(flag > 1)));
    if kr < 2 | kr > nr-2 | kz < 2 | kz > nz-2
      outofrange = true;
      if iteration_counter == 0
        disp(['warning in gsdesign: LOCKED rbdef = ' num2str(r), ...
          ', zbdef = ' num2str(z) ' is out of range, and will be ignored'])
      end
    elseif dlmin > min(dr,dz)/100 & ~isinpoly(r,z)
      if iteration_counter == 0
        disp(['warning in gsdesign: LOCKED rbdef = ' num2str(r), ...
          ', zbdef = ' num2str(z) ' is outside limiter, and will be ignored'])    
      end
    else
      wsep = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      iisep = kz+nz*(kr-1)+neighbors';
      psisep = wsep*psizr(iisep);
      dpsisepdx = wsep*dpsizrdx(iisep,:);
      Ae(end+1,:) = (dpsisepdx-dpsibrydx)*dxdxcirc(:,iopen);
      if fluxerror < 0.01
        be(end+1,1) = psibry+dpsibrydx*dxlocked-psisep-dpsisepdx*dxlocked;
      else
        be(end+1,1) = psibry-psisep;
      end
      % Bad idea to include x points as part of assignment! Hence if 0
      if 0 & dlmin > min(dr,dz)/100 % Need a null at this point
 	wr0 = [1 tr tr^2 tr^3]*mx;
	wz0 = mx'*[1 tz tz^2 tz^3]';
	wr1 = [0 1 2*tr 3*tr^2]*mx/dr;
	wz1 = mx'*[0 1 2*tz 3*tz^2]'/dz;
	wxr = reshape(wz0*wr1, 1, 16);
	wxz = reshape(wz1*wr0, 1, 16);
	dpsixrdx = wxr*dpsizrdx(iisep,:);
	dpsixzdx = wxz*dpsizrdx(iisep,:);
	psixr = wxr*psizr(iisep) + dpsixrdx*dxlocked;
	psixz = wxz*psizr(iisep) + dpsixzdx*dxlocked;
        Ae(end+1,:) = dpsixrdx*dxdxcirc(:,iopen);
	be(end+1,1) = -psixr;
        Ae(end+1,:) = dpsixzdx*dxdxcirc(:,iopen);
	be(end+1,1) = -psixz;
      end
      bdef.target.exists = true; % exists and inside limiter or very close
      bdef.target.dlmin = dlmin; % minimum distance to limiter
      bdef.target.r = r;
      bdef.target.z = z;
      bdef.target.w = wsep;
      bdef.target.ii = iisep;
      bdef.target.psi = psisep;
      bdef.target.dpsidx = dpsisepdx;            
    end
  end
end
if isfield(locks,'cpasma')
  dcpasmadpv = dcpasmadx*dxdxcirc(:,iopen);
  Ae(end+1,:) = dcpasmadpv;
  be(end+1,1) = locks.cpasma-cpasma-dcpasmadx*dxlocked;  
end
if isfield(locks,'li')
  dlidpv = dlidx*dxdxcirc(:,iopen);
  Ae(end+1,:) = dlidpv;
  be(end+1,1) = locks.li-li-dlidx*dxlocked;  
end
if isfield(locks,'betap')
  dbetapdpv = dbetapdx*dxdxcirc(:,iopen);
  Ae(end+1,:) = dbetapdpv;
  be(end+1,1) = locks.betap-betap-dbetapdx*dxlocked;  
end
if isfield(locks,'betan')
  dbetandpv = dbetandx*dxdxcirc(:,iopen);
  Ae(end+1,:) = dbetandpv;
  be(end+1,1) = locks.betan-betan-dbetandx*dxlocked;  
end
if isfield(locks,'psibry')
  dpsibrydpv = dpsibrydx*dxdxcirc(:,iopen);
  Ae(end+1,:) = dpsibrydpv;
  be(end+1,1) = locks.psibry-psibry-dpsibrydx*dxlocked;  
end
if isfield(locks,'psimag')
  dpsimagdpv = dpsimagdx*dxdxcirc(:,iopen);
  Ae(end+1,:) = dpsimagdpv;
  be(end+1,1) = locks.psimag-psimag-dpsimagdx*dxlocked;  
end
if isfield(locks,'psipla')
  dpsipladpv = dpsipladx*dxdxcirc(:,iopen);
  Ae(end+1,:) = dpsipladpv;
  be(end+1,1) = locks.psipla-psipla-dpsipladx*dxlocked;  
end
if isfield(locks,'fl')
  if ~exist('fl','var')
    if iteration_counter == 0
      disp('Warning gsdesign: locks.fl not applied since no fl found')
    end
  else
    n = min(length(locks.fl),length(fl));
    for i = 1:n
      dfldpv = dfldx(i,:)*dxdxcirc(:,iopen);
      Ae(end+1,:) = dfldpv;
      be(end+1,1) = locks.fl(i)-fl(i)-dfldx(i,:)*dxlocked;
    end
  end
end
if isfield(locks,'bp')
  if ~exist('bp','var')
    if iteration_counter == 0
      disp('Warning gsdesign: locks.bp not applied since no bp found')
    end
  else
    n = min(length(locks.bp),length(bp));
    for i = 1:n
      dbpdpv = dbpdx(i,:)*dxdxcirc(:,iopen);
      Ae(end+1,:) = dbpdpv;
      be(end+1,1) = locks.bp(i)-bp(i)-dbpdx(i,:)*dxlocked;
    end
  end
end
if isfield(locks,'rog')
  if ~exist('rog','var')
    if iteration_counter == 0
      disp('Warning gsdesign: locks.rog not applied since no rog found')
    end
  else
    n = min(length(locks.rog),length(rog));
    for i = 1:n
      drogdpv = drogdx(i,:)*dxdxcirc(:,iopen);
      Ae(end+1,:) = drogdpv;
      be(end+1,1) = locks.rog(i)-rog(i)-drogdx(i,:)*dxlocked;
    end
  end
end
if isfield(locks,'psic')
  n = min(length(locks.psic),nc);
  for i = 1:n
    k = indic(i);
    Ae(end+1,:) = dysdx(k,:)*dxdxcirc(:,iopen);
    be(end+1,1) = locks.psic(i)-ys(k)-dysdx(k,:)*dxlocked;
  end
end
if isfield(locks,'psiv')
  n = min(length(locks.psiv),nv);
  for i = 1:n
    k = indiv(i);
    Ae(end+1,:) = dysdx(k,:)*dxdxcirc(:,iopen);
    be(end+1,1) = locks.psiv(i)-ys(k)-dysdx(k,:)*dxlocked;
  end
end
if isfield(locks,'frc')
  n = min(length(locks.frc),nc);
  for i = 1:n
    Ae(end+1,:) = dfrcdx(i,:)*dxdxcirc(:,iopen);
    be(end+1,1) = locks.frc(i)-frc(i)-dfrcdx(i,:)*dxlocked;
  end
end
if isfield(locks,'fzc')
  n = min(length(locks.fzc),nc);
  for i = 1:n
    Ae(end+1,:) = dfzcdx(i,:)*dxdxcirc(:,iopen);
    be(end+1,1) = locks.fzc(i)-fzc(i)-dfzcdx(i,:)*dxlocked;
  end
end
ie = ~isnan(be);

if isfield(limits,'ci')
  dcidxcirc = eye(nci,nxcirc);
  dcidpv = dcidxcirc(:,iopen);
  for i = 1:size(limits.ci,1)
    Ai(end+1,:) = -dcidpv(i,:);
    bi(end+1,1) = -(limits.ci(i,1)-ci(i));
    Ai(end+1,:) = +dcidpv(i,:);
    bi(end+1,1) = +(limits.ci(i,2)-ci(i));
  end
end
if isfield(limits,'rsep') & isfield(limits,'zsep')
  n = min(size(limits.rsep,1),size(limits.zsep,1));
  for i = 1:n
    rs = linspace(limits.rsep(i,1),limits.rsep(i,2));
    zs = linspace(limits.zsep(i,1),limits.zsep(i,2));
    [is, ws] = gs_interp2(rg,zg,psizr,rs,zs,'WEIGHTS');
    ps = sum(ws'.*psizr(is)');
    [~, i1] = min(ps); % Upper limit for flux at i1 is psibry
    Ai(end+1,:) = +(ws(i1,:)*dpsizrdx(is(i1,:),:)-dpsibrydx)*dxdxcirc(:,iopen);
    bi(end+1,1) = +(psibry-ps(i1));
    [~, i2] = max(ps); % Lower limit for flux at i2 is psibry
    Ai(end+1,:) = -(ws(i2,:)*dpsizrdx(is(i2,:),:)-dpsibrydx)*dxdxcirc(:,iopen);
    bi(end+1,1) = -(psibry-ps(i2));
  end
end
if isfield(limits,'rx') & isfield(limits,'zx')
  gs_nulls_xlim
  m = nulls_xlim.count;
  if m > 0
    rxs = nulls_xlim.r(1:m);
    zxs = nulls_xlim.z(1:m);
    if fluxerror < 0.1
      % Correct positions due to dxlocked (which includes fluxerror correction)
      dpsizr(:) = dpsizrdx*dxlocked;
      for k = 1:m
	ii = nulls_xlim.ii(k,:)';
	rxs(k) = rxs(k) + nulls_xlim.drdpsi(k,:)*dpsizr(ii);
	zxs(k) = zxs(k) + nulls_xlim.dzdpsi(k,:)*dpsizr(ii);
      end
    end
    n = size(limits.rx,1);
    for i = 1:n
      ilimx = ~isnan(limits.rx(i,:));
      rs = limits.rx(i,ilimx);
      zs = limits.zx(i,ilimx);
      ks = inpolygon(rxs,zxs,rs,zs);
      rsm = mean(rs(1:end-1));
      zsm = mean(zs(1:end-1));
      [~,km] = sort((rxs-rsm).^2+(zxs-zsm).^2);
      if any(ks)
        xPointInsidePolygon = true;
	% Find first km for which ks is true
	k = km(min(find(ks(km)))); % index to x-point nearest r,z inside rs,zs
      else
        xPointInsidePolygon = false;
	% Find first km for which ks is true
	k = km(1); % index to x-point nearest r,z
      end
      r = nulls_xlim.r(k);
      z = nulls_xlim.z(k);
      ii = nulls_xlim.ii(k,:);
      drdx = nulls_xlim.drdpsi(k,:)*dpsizrdx(ii,:);
      dzdx = nulls_xlim.dzdpsi(k,:)*dpsizrdx(ii,:);
      xB = [rs(1:end-1)' rs(2:end)'];
      yB = [zs(1:end-1)' zs(2:end)'];
      if xPointInsidePolygon
        % Limit x-point displacement in 16 directions
	for j = 1:16
	  pr = cos(j*pi/8);
	  pz = sin(j*pi/8);
	  xA = r+[0 pr];
	  yA = z+[0 pz];
	  [~, ~, flag, tA, tB] = line_intersection(xA, yA, xB, yB);
	  if any(flag > 3 & tA > 0)
            Ai(end+1,:) = (pr*drdx+pz*dzdx)*dxdxcirc(:,iopen);
            bi(end+1,1) = max(tA(flag > 3 & tA > 0));      
	  end
        end
      else % x-point outside polygon limits, find closest distance
        xA = r+diff(zs)'*[0 1];
	yA = z-diff(rs)'*[0 1];
	[rr, zz, flag, tA, tB] = line_intersection(xA, yA, xB, yB);
	rt = [rr(flag>3)' rs];
	zt = [zz(flag>3)' zs];
	[dmin,k] = min(sqrt((rt-r).^2+(zt-z).^2));
	if dmin > 0
	  pr = (rt(k)-r)/dmin;
	  pz = (zt(k)-z)/dmin;
	  xA = r+[0 pr];
	  yA = z+[0 pz];
	  [~, ~, flag, tA, tB] = line_intersection(xA, yA, xB, yB);
	  if sum(flag > 3 & tA > 0) > 1
            Ai(end+1,:) = (pr*drdx+pz*dzdx)*dxdxcirc(:,iopen);
            bi(end+1,1) = max(tA(flag > 3 & tA > 0));      
            Ai(end+1,:) = -(pr*drdx+pz*dzdx)*dxdxcirc(:,iopen);
            bi(end+1,1) = -min(tA(flag > 3 & tA > 0));      
	  end
	end
      end
    end
  end
end
wrongway.r = [];
wrongway.z = [];
if bdef.target.exists
  % Limit flux change at critical points to get correct bdef point
  psiab = psimag - psibry;
  for i = 1:bdefs.count
    r = bdefs.r(i);
    z = bdefs.z(i);
    if isfield(limits,'bdef_dpsibar')
      bdef_dpsibar = limits.bdef_dpsibar;
    else
      bdef_dpsibar = 0.01;
    end
    [~, k] = min((rbbbs(1:nbbbs)-r).^2+(zbbbs(1:nbbbs)-z).^2);
    % db is distance to boundary from critical point
    db = sqrt((rbbbs(k)-r)^2+(zbbbs(k)-z)^2);
    % dt is distance to bdef target from critical point
    dt = sqrt((bdef.target.r-r)^2+(bdef.target.z-z)^2);
    if db < 3*dr & dt > 9*dr
      % This critical point must be avoided
      % Normalized flux must stay >1 if change monotonic from boundary
      if all(psiab*diff(gs_interp2(rg,zg,psizr,...
        r+(rbbbs(k)-r)*linspace(0,1,9),...
	z+(zbbbs(k)-z)*linspace(0,1,9))) > 0)
      dpsidx = bdefs.w(i,:)*dpsizrdx(bdefs.ii(i,:),:);
      psidiff = bdefs.psi(i) - bdef.target.psi;
      psicorr = (dpsidx - bdef.target.dpsidx)*dxlocked;
      Ai(end+1,:) = sign(psiab)*(dpsidx-bdef.target.dpsidx)*dxdxcirc(:,iopen);
      bi(end+1,1) = -sign(psiab)*(psidiff+psicorr+psiab*bdef_dpsibar);
      wrongway.r(end+1) = r;
      wrongway.z(end+1) = z;
      end
    end
  end
  r = bdef.target.r;
  z = bdef.target.z;
  [dum, i] = min(((rbbbs(1:nbbbs)-r)/dr).^2+((zbbbs(1:nbbbs)-z)/dz).^2);
  if [rmaxis-r zmaxis-z]*[rbbbs(i)-r;zbbbs(i)-z] > 0 & ...
     (rbbbs(i)-r)^2+(zbbbs(i)-z)^2 > dr^2
    % Ensure monotonic flux between boundary and bdef target to pull boundary
    n = 1+round(sqrt(dum)*2);
    rs = r + (rbbbs(i)-r)*(1:n)/n;
    zs = z + (zbbbs(i)-z)*(1:n)/n;
    [is, ws] = gs_interp2(rg,zg,psizr,rs,zs,'WEIGHTS');
    for i = 1:n
      dpsidx = ws(i,:)*dpsizrdx(is(i,:),:);
      psidiff = ws(i,:)*psizr(is(i,:)') - bdef.target.psi;
      psicorr = (dpsidx - bdef.target.dpsidx)*dxlocked;
      Ai(end+1,:) = -sign(psiab)*(dpsidx-bdef.target.dpsidx)*dxdxcirc(:,iopen);
      bi(end+1,1) = sign(psiab)*(psidiff+psicorr);
    end
  end
end
if isfield(limits,'cpasma')
  dcpasmadpv = dcpasmadx*dxdxcirc(:,iopen);
  Ai(end+1,:) = -dcpasmadpv;
  bi(end+1,1) = -(limits.cpasma(1)-cpasma);
  Ai(end+1,:) = +dcpasmadpv;
  bi(end+1,1) = +(limits.cpasma(2)-cpasma);  
end
if isfield(limits,'li')
  dlidpv = dlidx*dxdxcirc(:,iopen);
  Ai(end+1,:) = -dlidpv;
  bi(end+1,1) = -(limits.li(1)-li);
  Ai(end+1,:) = +dlidpv;
  bi(end+1,1) = +(limits.li(2)-li);  
end
if isfield(limits,'betap')
  dbetapdpv = dbetapdx*dxdxcirc(:,iopen);
  Ai(end+1,:) = -dbetapdpv;
  bi(end+1,1) = -(limits.betap(1)-betap);
  Ai(end+1,:) = +dbetapdpv;
  bi(end+1,1) = +(limits.betap(2)-betap);  
end
if isfield(limits,'betan')
  dbetandpv = dbetandx*dxdxcirc(:,iopen);
  Ai(end+1,:) = -dbetandpv;
  bi(end+1,1) = -(limits.betan(1)-betan);
  Ai(end+1,:) = +dbetandpv;
  bi(end+1,1) = +(limits.betan(2)-betan);  
end
if isfield(limits,'psibry')
  dpsibrydpv = dpsibrydx*dxdxcirc(:,iopen);
  Ai(end+1,:) = -dpsibrydpv;
  bi(end+1,1) = -(limits.psibry(1)-psibry);
  Ai(end+1,:) = +dpsibrydpv;
  bi(end+1,1) = +(limits.psibry(2)-psibry);  
end
if isfield(limits,'psimag')
  dpsimagdpv = dpsimagdx*dxdxcirc(:,iopen);
  Ai(end+1,:) = -dpsimagdpv;
  bi(end+1,1) = -(limits.psimag(1)-psimag);
  Ai(end+1,:) = +dpsimagdpv;
  bi(end+1,1) = +(limits.psimag(2)-psimag);  
end
if isfield(limits,'psipla')
  dpsipladpv = dpsipladx*dxdxcirc(:,iopen);
  Ai(end+1,:) = -dpsipladpv;
  bi(end+1,1) = -(limits.psipla(1)-psipla);
  Ai(end+1,:) = +dpsipladpv;
  bi(end+1,1) = +(limits.psipla(2)-psipla);  
end
if isfield(limits,'fl')
  if ~exist('fl','var')
    if iteration_counter == 0
      disp('Warning gsdesign: limits.fl not applied since no fl found')
    end
  else
    n = min(size(limits.fl,1),length(fl));
    for i = 1:n
      dfldpv = dfldx(i,:)*dxdxcirc(:,iopen);
      Ai(end+1,:) = -dfldpv;
      bi(end+1,1) = -(limits.fl(i,1)-fl(i));
      Ai(end+1,:) = +dfldpv;
      bi(end+1,1) = +(limits.fl(i,2)-fl(i));
    end
  end
end
if isfield(limits,'bp')
  if ~exist('bp','var')
    if iteration_counter == 0
      disp('Warning gsdesign: limits.bp not applied since no bp found')
    end
  else
    n = min(size(limits.bp,1),length(bp));
    for i = 1:n
      dbpdpv = dbpdx(i,:)*dxdxcirc(:,iopen);
      Ai(end+1,:) = -dbpdpv;
      bi(end+1,1) = -(limits.bp(i,1)-bp(i));
      Ai(end+1,:) = +dbpdpv;
      bi(end+1,1) = +(limits.bp(i,2)-bp(i));
    end
  end
end
if isfield(limits,'rog')
  if ~exist('rog','var')
    if iteration_counter == 0
      disp('Warning gsdesign: limits.rog not applied since no rog found')
    end
  else
    n = min(size(limits.rog,1),length(rog));
    for i = 1:n
      drogdpv = drogdx(i,:)*dxdxcirc(:,iopen);
      Ai(end+1,:) = -drogdpv;
      bi(end+1,1) = -(limits.rog(i,1)-rog(i));
      Ai(end+1,:) = +drogdpv;
      bi(end+1,1) = +(limits.rog(i,2)-rog(i));
    end
  end
end
if isfield(limits,'psic')
  n = min(size(limits.psic,1),nc);
  for i = 1:n
    k = indic(i);
    dpsicdpv = dysdx(k,:)*dxdxcirc(:,iopen);
    Ai(end+1,:) = -dpsicdpv;
    bi(end+1,1) = -(limits.psic(i,1)-ys(k));
    Ai(end+1,:) = +dpsicdpv;
    bi(end+1,1) = +(limits.psic(i,2)-ys(k));
  end
end
if isfield(limits,'psiv')
  n = min(size(limits.psiv,1),nv);
  for i = 1:n
    k = indiv(i);
    dpsivdpv = dysdx(k,:)*dxdxcirc(:,iopen);
    Ai(end+1,:) = -dpsivdpv;
    bi(end+1,1) = -(limits.psiv(i,1)-ys(k));
    Ai(end+1,:) = +dpsivdpv;
    bi(end+1,1) = +(limits.psiv(i,2)-ys(k));
  end
end
if isfield(limits,'frc')
  n = min(size(limits.frc,1),nc);
  for i = 1:n
    dfrcdpv = dfrcdx(i,:)*dxdxcirc(:,iopen);
    Ai(end+1,:) = -dfrcdpv;
    bi(end+1,1) = -(limits.frc(i,1)-frc(i));
    Ai(end+1,:) = +dfrcdpv;
    bi(end+1,1) = +(limits.frc(i,2)-frc(i));
  end
end
if isfield(limits,'fzc')
  n = min(size(limits.fzc,1),nc);
  for i = 1:n
    dfzcdpv = dfzcdx(i,:)*dxdxcirc(:,iopen);
    Ai(end+1,:) = -dfzcdpv;
    bi(end+1,1) = -(limits.fzc(i,1)-fzc(i));
    Ai(end+1,:) = +dfzcdpv;
    bi(end+1,1) = +(limits.fzc(i,2)-fzc(i));
  end
end
ii = any(Ai,2) & ~isinf(bi) & ~isnan(bi);

if license('checkout','optimization_toolbox')
  opts = optimset(...
    'Algorithm','interior-point-convex',...
    'Display','off',...
    'TolX',1e-12,...
    'TolFun',1e-12);
  xqp1 = quadprog((H+H')/2,h,Ai(ii,:),bi(ii),Ae(ie,:),be(ie),[],[],[],opts);
  dxcirc(iopen) = xqp1;
else
  for i = 1:length(be)
    Ai(end+1,:) = +Ae(i,:);
    bi(end+1,1) = +be(i,1);
    Ai(end+1,:) = -Ae(i,:);
    bi(end+1,1) = -be(i,1);
  end
  ii = any(Ai,2) & ~isinf(bi);
  rep = pdipmqpneq2(H,h,Ai(ii,:),bi(ii),100,1e-9,0.97);
  dxcirc(iopen) = rep.x;
end
dx = dxdxcirc*dxcirc;

dev = devdx*dx;

if outofrange & iteration_counter == 0
  disp(['Allowed r,z are: ' num2str(rg(2)) ' < r < ' num2str(rg(nr-1)), ...
    ' & ' num2str(zg(2)) ' < z < ' num2str(zg(nz-1))])
end
