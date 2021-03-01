function M = mcircle2point(r,z,rp,zp,a,nrh,nt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   M = mcircle2point(r,z,rp,zp,a);
%
%  PURPOSE: Return poloidal flux from circular plasma (or other 
%           axisymmetric conductor) carrying a uniform current density
%
%  INPUTS: r,z major radius and height of field point
%          rp,zp major radius and height of plasma center
%          a, minor radius of plasma          
%
%  OUTPUTS: M, flux for 1-A of current, a.k.a. as mutual inductance
%
%  METHOD: Flux is found in look-up table. Relative error always < 1e-6
%          See also mfcircle2point, mpcircle2point
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
%  VERSION %W% %G%
%
%  WRITTEN BY: Anders Welander ON 2015-03-11
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent a1 da na r1 dr nr
persistent c00 c01 c02 c03 c10 c11 c12 c13 c20 c21 c22 c23 c30 c31 c32 c33

if nargin <= 5 % If true then use look-up table with fix for a/rp < 0.04

  % Check if look-up table has been filled in by checking a1
  if isempty(a1)
    filename = which('mcircle2point.mat');
    load(filename)
  end

  r = r./rp; % if only rp was non-scalar, r will now change dimensions
  a = a./rp; % if only rp was non-scalar, a will now change dimensions

  % Interpolation is relatively very fast so just do it even if a is too small
  % A fix comes below for a < 0.04

  fa = (a(:)-a1)/da+1;
  i = floor(fa);
  i(i < 1) = 1;
  i(i > na) = na;
  x = fa-i;

  fr = (r(:)-r1)/dr+1;
  j = floor(fr);
  j(j < 1) = 1;
  j(j > nr) = nr;
  y = fr-j;

  k = i+(j-1)*na;

  % Allow r, rp, a to have arbitrary dimensions as long as no mixing
  if length(r(:)) > length(a(:)) % rp already influenced dimensions of these
    M = r; % Make M have right dimension, in this case a was scalar
  else
    M = a; % Make M have right dimension, here r was same as a or scalar
  end
  x2 = x.*x;
  y2 = y.*y;
  x3 = x.*x2;
  y3 = y.*y2;
  M(:) = rp(:).*(...
    c00(k)     + c10(k).*x     + c20(k).*x2     + c30(k).*x3 + ...
    c01(k).*y  + c11(k).*x.*y  + c21(k).*x2.*y  + c31(k).*x3.*y + ...
    c02(k).*y2 + c12(k).*x.*y2 + c22(k).*x2.*y2 + c32(k).*x3.*y2 + ...
    c03(k).*y3 + c13(k).*x.*y3 + c23(k).*x2.*y3 + c33(k).*x3.*y3);  

  % Phase in a good solution for a < 0.04
  if length(a(:)) == 1
    if a > 0.04
      return
    elseif a > 0.03
      f = polyval([6 -15 10 0 0 0],(a-0.03)/0.01); % Keep this much of solution
    else
      f = 0; % Don't use interpolation result at all
    end
    m = mutinds(1,0,r(:),0);
    rh = r(:)-1;
    if any(abs(rh) < a)
      [mi,mir] = mutinds(1,0,1-a,0);
      [mo,mor] = mutinds(1,0,1+a,0);
      % Let m be polynomial in x=(1-r): c0+c1*x+c2*x^2+c3*x^3
      % c0-c1*a+c2*a^2-c3*a^3 = mi
      % c1-2*c2*a+3*c3*a^2 = mir
      % c0+c1*a+c2*a^2+c3*a^3 = mo
      % c1+2*c2*a+3*c3*a^2 = mor
      % p = inv([-1 1 -1 1; 3 -2 1 0; 1 1 1 1; 3 2 1 0])*[mi; a*mir; mo; a*mor];
      p = [1 1 -1 1;0 -1 0 1; -3 -1 3 -1; 2 1 2 -1]*[mi; a*mir; mo; a*mor];
      for i = 1:length(rh(:))
	if abs(rh(i)) < a
          m(i) = polyval(p,rh(i)/a)/4;
	end
      end
    end
    M(:) = f*M(:)+(1-f)*rp(:).*m;
  else
    % Here length(a(:)) > 1
    for i = 1:length(a(:))
      if a(i) < 0.04
	if a(i) > 0.03
          % Keep fraction f of interpolation result
          f = polyval([6 -15 10 0 0 0],(a(i)-0.03)/0.01);
	else
          f = 0; % Don't use interpolation result at all
	end
	if length(r(:)) == 1
          rh = r-1;
	else
          rh = r(i)-1;
	end
	if abs(rh) < a(i) % Solution for inside the small ring
          [mi,mir] = mutinds(1,0,1-a(i),0);
          [mo,mor] = mutinds(1,0,1+a(i),0);
          p = [1 1 -1 1;0 -1 0 1;-3 -1 3 -1;2 1 2 -1]*...
	    [mi; a(i)*mir; mo; a(i)*mor];
          m = polyval(p,rh/a(i))/4;
	else % Solution for outside the small ring
          m = mutinds(1,0,r(i),0);
	end
	if length(rp(:)) == 1
          M(i) = f*M(i)+(1-f)*rp*m;
	else
          M(i) = f*M(i)+(1-f)*rp(i).*m;
	end
      end
    end
  end

else % Calculate mutuals the slow way

  if ~exist('nrh','var') | isempty(nrh)
    nrh = 100;
  end
  if ~exist('nt','var') | isempty(nt)
    nt = 100;
  end

  [n,k] = max([length(r(:)) length(z(:)) length(rp(:)) length(zp(:)) length(a(:))]);

  % Give M right dimensions
  if k == 1
    M = r;
  elseif k == 2
    M = z;
  elseif k == 3
    M = rp;
  elseif k == 4
    M = zp;
  else
    M = a;
  end
  Mr = M;
  Mz = M;
  Mrr = M;
  Mrz = M;
  Mzz = M;
  Mrrr = M;
  Mrrz = M;
  Mrzz = M;
  Mzzz = M;
  % Convert any scalar input to same length as others
  if length(r(:)) < n
    r = r(1)+zeros(n,1);
  end
  if length(z(:)) < n
    z = z(1)+zeros(n,1);
  end
  if length(rp(:)) < n
    rp = rp(1)+zeros(n,1);
  end
  if length(zp(:)) < n
    zp = zp(1)+zeros(n,1);
  end
  if length(a(:)) < n
    a = a(1)+zeros(n,1);
  end

  for i = 1:n
    % Calculating geometric quantities
    rh = a(i)*(1:nrh)'/nrh; % Minor radii
    t = [1:2:2*nt-1]/nt*pi; % Poloidal angles
    ct = cos(t);
    st = sin(t);
    drh = a(i)/nrh;
    dt = 2*pi/nt;
    rs = rp(i) + rh*ct;
    zs = zp(i) + rh*st;
    Imdr = [(rh(1)/2+drh/6); rh(2:nrh-1); (rh(nrh)/2-drh/6)]*(drh*dt);
    fA = Imdr*ones(1,nt);
    iring = sum(fA(:));
    [m,mr,mz,mrr,mrz,mzz,mrrr,mrrz,mrzz,mzzz] = mutinds(rs,zs,r(i),z(i));
    if any(isnan(m))
      m(isnan(m)) = 0;
      mr(isnan(mr)) = 0;
      mz(isnan(mz)) = 0;
      mrr(isnan(mrr)) = 0;
      mrz(isnan(mrz)) = 0;
      mzz(isnan(mzz)) = 0;
      mrrr(isnan(mrrr)) = 0;
      mrrz(isnan(mrrz)) = 0;
      mrzz(isnan(mrzz)) = 0;
      mzzz(isnan(mzzz)) = 0;
    end
    M(i) = sum(fA(:).*m(:))/iring;
    Mr(i) = sum(fA(:).*mr(:))/iring;
    Mz(i) = sum(fA(:).*mz(:))/iring;
    Mrr(i) = sum(fA(:).*mrr(:))/iring;
    Mrz(i) = sum(fA(:).*mrz(:))/iring;
    Mzz(i) = sum(fA(:).*mzz(:))/iring;
    Mrrr(i) = sum(fA(:).*mrrr(:))/iring;
    Mrrz(i) = sum(fA(:).*mrrz(:))/iring;
    Mrzz(i) = sum(fA(:).*mrzz(:))/iring;
    Mzzz(i) = sum(fA(:).*mzzz(:))/iring;
  end
end
