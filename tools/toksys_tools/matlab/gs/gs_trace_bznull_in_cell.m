function [jx, ix, m, vjx, vix] = ...
  gs_trace_bznull_in_cell(jr, iz, j1, j2, i1, i2, y44, x, y, vj, vi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   [jx,ix,m,vjx,vix] = gs_trace_bznull_in_cell(jr,iz,j1,j2,i1,i2,y44);
%
%  PURPOSE: Find where the contour bz = 0 exits the rectangle j1, j2, i1, i2
%
%  INPUTS: jr, iz, entry point for the contour bz = 0 (both values >=0 & <=1)
%          j1, j2, min and max r of the rectangle (both values >=0 & <=1)
%          i1, i2, min and max z of the rectangle (both values >=0 & <=1)
%          y44, values of flux for jr,iz = -1, 0, 1, 2
%
%  OUTPUTS:  jx, ix = r, z of exit point
%            m, where exit is, 0=x, -1=bottom, +1=top, -4=in, +4=out
%            vjx, vix = vectors of j,i points recorded within the cell
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
%  METHOD:  If there is only one place to exit then pick it, otherwise
%           call this function recursively to walk through in steps

%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander  ON	3/31/14
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2; % For interpolation

if exist('vj','var')
  vjx = vj;
else
  vjx = [];
end

if exist('vi','var')
  vix = vi;
else
  vix = [];
end


% Sacrifice some precision to avoid recursion limit and save time
if i2-i1 < 1e-8 & j2-j1 < 1e-8
  jx = jr;
  ix = iz;
  m = 0;
  return
end


if ~exist('y','var')

  % Calculate roots to br = 0 along edges.

  wr0 = [1 j1 j1^2 j1^3]*mx;
  a = reshape(mx'*[0 0 0 3]'*wr0,1,16)*y44(:);
  b = reshape(mx'*[0 0 2 0]'*wr0,1,16)*y44(:);
  c = reshape(mx'*[0 1 0 0]'*wr0,1,16)*y44(:);
  y.vi = (-b+[-1 1]*sqrt(b^2-4*a*c))/(2*a);

  wr0 = [1 j2 j2^2 j2^3]*mx;
  a = reshape(mx'*[0 0 0 3]'*wr0,1,16)*y44(:);
  b = reshape(mx'*[0 0 2 0]'*wr0,1,16)*y44(:);
  c = reshape(mx'*[0 1 0 0]'*wr0,1,16)*y44(:);
  y.vo = (-b+[-1 1]*sqrt(b^2-4*a*c))/(2*a);

  wz1 = mx'*[0 1 2*i1 3*i1^2]';
  a = reshape(wz1*[0 0 0 1]*mx,1,16)*y44(:);
  b = reshape(wz1*[0 0 1 0]*mx,1,16)*y44(:);
  c = reshape(wz1*[0 1 0 0]*mx,1,16)*y44(:);
  d = reshape(wz1*[1 0 0 0]*mx,1,16)*y44(:);
  y.hl = roots([a b c d])';

  wz1 = mx'*[0 1 2*i2 3*i2^2]';
  a = reshape(wz1*[0 0 0 1]*mx,1,16)*y44(:);
  b = reshape(wz1*[0 0 1 0]*mx,1,16)*y44(:);
  c = reshape(wz1*[0 1 0 0]*mx,1,16)*y44(:);
  d = reshape(wz1*[1 0 0 0]*mx,1,16)*y44(:);
  y.hu = roots([a b c d])';

end

r10 = ~[imag(y.hl) imag(y.vi) imag(y.hu) imag(y.vo)];
k10 = [y.hl >= j1 & y.hl < j2,  y.vi >= i1 & y.vi < i2, ...
       y.hu >= j1 & y.hu < j2,  y.vo >= i1 & y.vo < i2];

brnull = any(r10 & k10);

if ~exist('x','var')

  % Calculate roots to bz = 0 along edges.

  wz0 = mx'*[1 i1 i1^2 i1^3]';
  a = reshape(wz0*[0 0 0 3]*mx,1,16)*y44(:);
  b = reshape(wz0*[0 0 2 0]*mx,1,16)*y44(:);
  c = reshape(wz0*[0 1 0 0]*mx,1,16)*y44(:);
  x.hl = (-b+[-1 1]*sqrt(b^2-4*a*c))/(2*a);

  wz0 = mx'*[1 i2 i2^2 i2^3]';
  a = reshape(wz0*[0 0 0 3]*mx,1,16)*y44(:);
  b = reshape(wz0*[0 0 2 0]*mx,1,16)*y44(:);
  c = reshape(wz0*[0 1 0 0]*mx,1,16)*y44(:);
  x.hu = (-b+[-1 1]*sqrt(b^2-4*a*c))/(2*a);

  wr1 = [0 1 2*j1 3*j1^2]*mx;
  a = reshape(([0 0 0 1]*mx)'*wr1,1,16)*y44(:);
  b = reshape(([0 0 1 0]*mx)'*wr1,1,16)*y44(:);
  c = reshape(([0 1 0 0]*mx)'*wr1,1,16)*y44(:);
  d = reshape(([1 0 0 0]*mx)'*wr1,1,16)*y44(:);
  x.vi = roots([a b c d])';

  wr1 = [0 1 2*j2 3*j2^2]*mx;
  a = reshape(([0 0 0 1]*mx)'*wr1,1,16)*y44(:);
  b = reshape(([0 0 1 0]*mx)'*wr1,1,16)*y44(:);
  c = reshape(([0 1 0 0]*mx)'*wr1,1,16)*y44(:);
  d = reshape(([1 0 0 0]*mx)'*wr1,1,16)*y44(:);
  x.vo = roots([a b c d])';

end

r10 = ~[imag(x.hl) imag(x.vi) imag(x.hu) imag(x.vo)];
k10 = [x.hl >= j1 & x.hl < j2,  x.vi >= i1 & x.vi < i2, ...
       x.hu >= j1 & x.hu < j2,  x.vo >= i1 & x.vo < i2];

f10 = logical(ones(1,10));

if iz == i1 % Entry point is on the lower edge
  [~, p] = min(abs(x.hl-jr));
  f10(p) = 0;
end
if jr == j1 % Entry point is on the inner edge
  [~, p] = min(abs(x.vi-iz));
  f10(p+2) = 0;
end
if iz == i2 % Entry point is on the upper edge
  [~, p] = min(abs(x.hu-jr));
  f10(p+5) = 0;
end
if jr == j2 % Entry point is on the outer edge
  [~, p] = min(abs(x.vo-iz));
  f10(p+7) = 0;
end
if sum(r10 & k10 & f10) == 0 % No place to exit
  jx = jr;
  ix = iz;
  m = 0;
elseif ~brnull & sum(r10 & k10 & f10) == 1 % Only one way to exit
  for p = 1:10
    if r10(p) & k10(p) & f10(p)
      if p < 3 % Going down
	jx = floor(j1) + x.hl(p);
	ix = i1;
	m = -1;
      elseif p < 6 % Going in
	jx = j1;
	ix = floor(i1) + x.vi(p-2);
	m = -4;
      elseif p < 8 % Going up
	jx = floor(j1) + x.hu(p-5);
	ix = i2;
	m = +1;
      else % Going out
	jx = j2;
	ix = floor(i1) + x.vo(p-7);
	m = +4;
      end
    end
  end
  vjx(end+1) = jx;
  vix(end+1) = ix;
else
  % More than one way to exit or brnull, that's complicated
  % Divide rectangle into 4 pieces numbered 1,2,3,4 and step through them
  imid = (i1+i2)/2;
  wz0 = mx'*[1 imid imid^2 imid^3]';
  a = reshape(wz0*[0 0 0 3]*mx,1,16)*y44(:);
  b = reshape(wz0*[0 0 2 0]*mx,1,16)*y44(:);
  c = reshape(wz0*[0 1 0 0]*mx,1,16)*y44(:);
  xhmid = (-b+[-1 1]*sqrt(b^2-4*a*c))/(2*a);
  wz1 = mx'*[0 1 2*imid 3*imid^2]';
  a = reshape(wz1*[0 0 0 1]*mx,1,16)*y44(:);
  b = reshape(wz1*[0 0 1 0]*mx,1,16)*y44(:);
  c = reshape(wz1*[0 1 0 0]*mx,1,16)*y44(:);
  d = reshape(wz1*[1 0 0 0]*mx,1,16)*y44(:);
  yhmid = roots([a b c d])';
  jmid = (j1+j2)/2;
  wr1 = [0 1 2*jmid 3*jmid^2]*mx;
  a = reshape(([0 0 0 1]*mx)'*wr1,1,16)*y44(:);
  b = reshape(([0 0 1 0]*mx)'*wr1,1,16)*y44(:);
  c = reshape(([0 1 0 0]*mx)'*wr1,1,16)*y44(:);
  d = reshape(([1 0 0 0]*mx)'*wr1,1,16)*y44(:);
  xvmid = roots([a b c d])';
  wr0 = [1 jmid jmid^2 jmid^3]*mx;
  a = reshape(mx'*[0 0 0 3]'*wr0,1,16)*y44(:);
  b = reshape(mx'*[0 0 2 0]'*wr0,1,16)*y44(:);
  c = reshape(mx'*[0 1 0 0]'*wr0,1,16)*y44(:);
  yvmid = (-b+[-1 1]*sqrt(b^2-4*a*c))/(2*a);
  jx = jr;
  ix = iz;
  irect = 1 + (iz >= imid) + 2*(jr >= jmid);
  done_tracing = false;
  while ~done_tracing
    x2 = x;
    y2 = y;
    if irect == 1
      x2.vo = xvmid;
      y2.vo = yvmid;
      x2.hu = xhmid;
      y2.hu = yhmid;
      [jx, ix, m, vjx, vix] = ...
        gs_trace_bznull_in_cell(jx, ix, j1, jmid, i1, imid, y44, x2, y2, vjx, vix);
      if m == 1
        irect = 2;
      elseif m == 4
        irect = 3;
      else
        done_tracing = true;
      end
    elseif irect == 2
      x2.vo = xvmid;
      y2.vo = yvmid;
      x2.hl = xhmid;
      y2.hl = yhmid;
      [jx, ix, m, vjx, vix] = ...
        gs_trace_bznull_in_cell(jx, ix, j1, jmid, imid, i2, y44, x2, y2, vjx, vix);
      if m == -1
        irect = 1;
      elseif m == 4
        irect = 4;
      else
        done_tracing = true;
      end
    elseif irect == 3
      x2.vi = xvmid;
      y2.vi = yvmid;
      x2.hu = xhmid;
      y2.hu = yhmid;
      [jx, ix, m, vjx, vix] = ...
        gs_trace_bznull_in_cell(jx, ix, jmid, j2, i1, imid, y44, x2, y2, vjx, vix);
      if m == 1
        irect = 4;
      elseif m == -4
        irect = 1;
      else
        done_tracing = true;
      end
    else % irect = 4
      x2.vi = xvmid;
      y2.vi = yvmid;
      x2.hl = xhmid;
      y2.hl = yhmid;
      [jx, ix, m, vjx, vix] = ...
        gs_trace_bznull_in_cell(jx, ix, jmid, j2, imid, i2, y44, x2, y2, vjx, vix);
      if m == -1
        irect = 3;
      elseif m == -4
        irect = 2;
      else
        done_tracing = true;
      end
    end
  end  
end
