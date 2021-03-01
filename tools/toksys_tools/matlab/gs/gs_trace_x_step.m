function [jx, ix, m, vjx, vix] = ...
  gs_trace_x_step(jr, iz, j1, j2, i1, i2, y44, vj, vi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   [jx,ix,m,vjx,vix] = gs_trace_x_step(jr,iz,j1,j2,i1,i2,y44,vj,vi);
%
%  PURPOSE: Find where the contour bz = 0 exits the rectangle j1, j2, i1, i2
%
%  INPUTS: jr, iz, entry point for the contour bz = 0 (both values >=0 & <=1)
%          j1, j2, min and max r of the rectangle (both values >=0 & <=1)
%          i1, i2, min and max z of the rectangle (both values >=0 & <=1)
%          y44, values of flux for jr,iz = -1, 0, 1, 2
%
%  OUTPUTS:  jx, ix, r, z of exit point
%            m, where exit is, 0=x, -1=bottom, +1=top, -4=in, +4=out
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
%  METHOD:  If there is only one place to exit then pick it, otherwise
%           call this function recursively to walk through in steps

%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander  ON	3/2/14
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
if i2-i1 < 1e-12 & j2-j1 < 1e-12
  jx = jr;
  ix = iz;
  m = 0;
  return
end

% Check for br = 0 along edges.

wr0 = [1 j1 j1^2 j1^3]*mx;
a = reshape(mx'*[0 0 0 3]'*wr0,1,16)*y44(:);
b = reshape(mx'*[0 0 2 0]'*wr0,1,16)*y44(:);
c = reshape(mx'*[0 1 0 0]'*wr0,1,16)*y44(:);
yvi = (-b+[-1 1]*sqrt(b^2-4*a*c))/(2*a);

wr0 = [1 j2 j2^2 j2^3]*mx;
a = reshape(mx'*[0 0 0 3]'*wr0,1,16)*y44(:);
b = reshape(mx'*[0 0 2 0]'*wr0,1,16)*y44(:);
c = reshape(mx'*[0 1 0 0]'*wr0,1,16)*y44(:);
yvo = (-b+[-1 1]*sqrt(b^2-4*a*c))/(2*a);

wz1 = mx'*[0 1 2*i1 3*i1^2]';
a = reshape(wz1*[0 0 0 1]*mx,1,16)*y44(:);
b = reshape(wz1*[0 0 1 0]*mx,1,16)*y44(:);
c = reshape(wz1*[0 1 0 0]*mx,1,16)*y44(:);
d = reshape(wz1*[1 0 0 0]*mx,1,16)*y44(:);
yhl = roots([a b c d])';

wz1 = mx'*[0 1 2*i2 3*i2^2]';
a = reshape(wz1*[0 0 0 1]*mx,1,16)*y44(:);
b = reshape(wz1*[0 0 1 0]*mx,1,16)*y44(:);
c = reshape(wz1*[0 1 0 0]*mx,1,16)*y44(:);
d = reshape(wz1*[1 0 0 0]*mx,1,16)*y44(:);
yhu = roots([a b c d])';

r10 = ~[imag(yhl) imag(yvi) imag(yhu) imag(yvo)];
k10 = [yhl >= j1 & yhl < j2,  yvi >= i1 & yvi < i2, ...
       yhu >= j1 & yhu < j2,  yvo >= i1 & yvo < i2];

brnull = any(r10 & k10);

wz0 = mx'*[1 i1 i1^2 i1^3]';
a = reshape(wz0*[0 0 0 3]*mx,1,16)*y44(:);
b = reshape(wz0*[0 0 2 0]*mx,1,16)*y44(:);
c = reshape(wz0*[0 1 0 0]*mx,1,16)*y44(:);
xhl = (-b+[-1 1]*sqrt(b^2-4*a*c))/(2*a);

wz0 = mx'*[1 i2 i2^2 i2^3]';
a = reshape(wz0*[0 0 0 3]*mx,1,16)*y44(:);
b = reshape(wz0*[0 0 2 0]*mx,1,16)*y44(:);
c = reshape(wz0*[0 1 0 0]*mx,1,16)*y44(:);
xhu = (-b+[-1 1]*sqrt(b^2-4*a*c))/(2*a);

wr1 = [0 1 2*j1 3*j1^2]*mx;
a = reshape(([0 0 0 1]*mx)'*wr1,1,16)*y44(:);
b = reshape(([0 0 1 0]*mx)'*wr1,1,16)*y44(:);
c = reshape(([0 1 0 0]*mx)'*wr1,1,16)*y44(:);
d = reshape(([1 0 0 0]*mx)'*wr1,1,16)*y44(:);
xvi = roots([a b c d])';

wr1 = [0 1 2*j2 3*j2^2]*mx;
a = reshape(([0 0 0 1]*mx)'*wr1,1,16)*y44(:);
b = reshape(([0 0 1 0]*mx)'*wr1,1,16)*y44(:);
c = reshape(([0 1 0 0]*mx)'*wr1,1,16)*y44(:);
d = reshape(([1 0 0 0]*mx)'*wr1,1,16)*y44(:);
xvo = roots([a b c d])';

r10 = ~[imag(xhl) imag(xvi) imag(xhu) imag(xvo)];
k10 = [xhl >= j1 & xhl < j2,  xvi >= i1 & xvi < i2, ...
       xhu >= j1 & xhu < j2,  xvo >= i1 & xvo < i2];
f10 = logical(ones(1,10));
if iz == i1 % Entry point is on the lower edge
  [~, p] = min(abs(xhl-jr));
  f10(p) = 0;
end
if jr == j1 % Entry point is on the inner edge
  [~, p] = min(abs(xvi-iz));
  f10(p+2) = 0;
end
if iz == i2 % Entry point is on the upper edge
  [~, p] = min(abs(xhu-jr));
  f10(p+5) = 0;
end
if jr == j2 % Entry point is on the outer edge
  [~, p] = min(abs(xvo-iz));
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
	jx = floor(j1) + xhl(p);
	ix = i1;
	m = -1;
      elseif p < 6 % Going in
	jx = j1;
	ix = floor(i1) + xvi(p-2);
	m = -4;
      elseif p < 8 % Going up
	jx = floor(j1) + xhu(p-5);
	ix = i2;
	m = +1;
      else % Going out
	jx = j2;
	ix = floor(i1) + xvo(p-7);
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
  jmid = (j1+j2)/2;
  jx = jr;
  ix = iz;
  irect = 1 + (iz >= imid) + 2*(jr >= jmid);
  done_tracing = false;
  while ~done_tracing
    if irect == 1
      [jx, ix, m, vjx, vix] = ...
        gs_trace_x_step(jx, ix, j1, jmid, i1, imid, y44, vjx, vix);
      if m == 1
        irect = 2;
      elseif m == 4
        irect = 3;
      else
        done_tracing = true;
      end
    elseif irect == 2
      [jx, ix, m, vjx, vix] = ...
        gs_trace_x_step(jx, ix, j1, jmid, imid, i2, y44, vjx, vix);
      if m == -1
        irect = 1;
      elseif m == 4
        irect = 4;
      else
        done_tracing = true;
      end
    elseif irect == 3
      [jx, ix, m, vjx, vix] = ...
        gs_trace_x_step(jx, ix, jmid, j2, i1, imid, y44, vjx, vix);
      if m == 1
        irect = 4;
      elseif m == -4
        irect = 1;
      else
        done_tracing = true;
      end
    else % irect = 4
      [jx, ix, m, vjx, vix] = ...
        gs_trace_x_step(jx, ix, jmid, j2, imid, i2, y44, vjx, vix);
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
