function [jx, ix, m, hit] = gs_trace_step(...
  jr, iz, j1, j2, i1, i2, y44, yb, xhl, xvi, xhu, xvo, rzhit);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   [jx,ix,m,hit] = gs_trace_step(...
%  jr,iz,j1,j2,i1,i2,y44,yb,xhl,xvi,xhu,xvo, rzhit);
%
%  PURPOSE: Find where contour y = yb exits rectangle j1, j2, i1, i2
%
%  INPUTS: jr, iz, entry point for the contour y = yb (both values >=0 & <=1)
%          j1, j2, min and max r of the rectangle (both values >=0 & <=1)
%          i1, i2, min and max z of the rectangle (both values >=0 & <=1)
%          y44, values of y for jr,iz = -1, 0, 1, 2
%          yb, value for which contour is sought
%          xhl, (optional) 3 values of j where y = yb for i = j1 (or nan)
%          xvi, (optional) 3 values of i where y = yb for j = j1 (or nan)
%          xhu, (optional) 3 values of j where y = yb for i = j2 (or nan)
%          xvo, (optional) 3 values of i where y = yb for j = j2 (or nan)
%          rzhit, (optional) [r1 r2 z1 z2] of region for hit test
%
%  OUTPUTS:  jx, ix, r, z of exit point
%            m, where exit is, 0=nowhere, -1=bottom, +1=top, -4=in, +4=out
%            hit, boolean true if contour goes through rzhit region
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
%  METHOD:  If there is only one place to exit then pick it, otherwise
%           call this function recursively to get through in several steps

%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander  ON	2/26/14
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2; % For interpolation

if ~exist('xhl','var') | isnan(xhl(1))
  xhl = gs_solve_hermite_cubic([1 i1 i1^2 i1^3]*mx*y44, yb)';
end
if ~exist('xvi','var') | isnan(xvi(1))
  xvi = gs_solve_hermite_cubic([1 j1 j1^2 j1^3]*mx*y44', yb)';
end
if ~exist('xhu','var') | isnan(xhu(1))
  xhu = gs_solve_hermite_cubic([1 i2 i2^2 i2^3]*mx*y44, yb)';
end
if ~exist('xvo','var') | isnan(xvo(1))
  xvo = gs_solve_hermite_cubic([1 j2 j2^2 j2^3]*mx*y44', yb)';
end

r12 = ~[imag(xhl) imag(xvi) imag(xhu) imag(xvo)];
k12 = [xhl >= j1 & xhl < j2,  xvi >= i1 & xvi < i2, ...
       xhu >= j1 & xhu < j2,  xvo >= i1 & xvo < i2];
f12 = logical(ones(1,12));
if iz == i1 % Entry point is on the lower edge
  [~, p] = min(abs(xhl-jr));
  f12(p) = 0;
end
if jr == j1 % Entry point is on the inner edge
  [~, p] = min(abs(xvi-iz));
  f12(p+3) = 0;
end
if iz == i2 % Entry point is on the upper edge
  [~, p] = min(abs(xhu-jr));
  f12(p+6) = 0;
end
if jr == j2 % Entry point is on the outer edge
  [~, p] = min(abs(xvo-iz));
  f12(p+9) = 0;
end

if exist('rzhit','var')
  if length(rzhit) ~= 4
    error('rzhit must have a length of 4')
  end
  rzhit = rzhit(:)';
  if rzhit(1) < j1 | rzhit(1) >= j2 | rzhit(1) >= rzhit(2) | ...
     rzhit(3) < i1 | rzhit(3) >= i2 | rzhit(3) >= rzhit(4)
    error('rzhit region must be inside region [j1 j2 i1 i2]')
  end
  jmid = mean(rzhit(1:2));
  imid = mean(rzhit(3:4));
  test4hit = true;
else
  test4hit = false;
end
hit = false;

if test4hit % Divide into smaller regions where the middle one is rzhit
  js = unique([j1 rzhit(1:2) j2]); % Typically 3 radial divisions
  is = unique([i1 rzhit(3:4) i2]); % Typically 3 vertical divisions
  nr = length(js)-1;
  nz = length(is)-1;
  j0 = max(find(js(1:end-1) <= jr & jr <= js(2:end))); % r-index to first div
  i0 = max(find(is(1:end-1) <= iz & iz <= is(2:end))); % z-index to first div
  jhit = 1 + (rzhit(1) > j1); % r-index to division that is tested for hits
  ihit = 1 + (rzhit(3) > i1); % z-index to division that is tested for hits
  jx = jr; % Will become the r at exit point
  ix = iz; % Will become the z at exit point
  done_tracing = false;
  while ~done_tracing
    [jx, ix, m] = gs_trace_step(...
      jx, ix, js(j0), js(j0+1), is(i0), is(i0+1), y44, yb);
    if j0 == jhit & i0 == ihit
      hit = true;
    end
    if     m == -4
      j0 = j0-1;
    elseif m == -1
      i0 = i0-1;
    elseif m == +1
      i0 = i0+1;
    elseif m == +4
      j0 = j0+1;
    else
      done_tracing = true;
    end
    if j0 < 1 | j0 > nr | i0 < 1 | i0 > nz
      done_tracing = true;
    end
  end  
elseif sum(r12 & k12 & f12) == 1 % There is only one way to exit
  for p = 1:3
    if r12(p) & k12(p) & f12(p) % Going down
      jx = floor(j1) + xhl(p);
      ix = i1;
      m = -1;
    end
    if r12(p+3) & k12(p+3) & f12(p+3) % Going in
      jx = j1;
      ix = floor(i1) + xvi(p);
      m = -4;
    end
    if r12(p+6) & k12(p+6) & f12(p+6) % Going up
      jx = floor(j1) + xhu(p);
      ix = i2;
      m = +1;
    end
    if r12(p+9) & k12(p+9) & f12(p+9) % Going out
      jx = j2;
      ix = floor(i1) + xvo(p);
      m = +4;
    end
  end
elseif sum(r12 & k12 & f12) > 1 % More than one way to exit, that's complicated
  % Divide this rectangle into 4 pieces numbered 1,2,3,4 and step through them
  imid = (i1+i2)/2;
  jmid = (j1+j2)/2;
  xhmid = gs_solve_hermite_cubic([1 imid imid^2 imid^3]*mx*y44, yb)';
  xvmid = gs_solve_hermite_cubic([1 jmid jmid^2 jmid^3]*mx*y44', yb)';
  jx = jr;
  ix = iz;
  irect = 1 + (iz >= imid) + 2*(jr >= jmid);
  done_tracing = false;
  while ~done_tracing
    if irect == 1
      [jx, ix, m] = ...
        gs_trace_step(jx, ix, j1, jmid, i1, imid, y44, yb, xhl, xvi, xhmid, xvmid);
      if m == 1
        irect = 2;
      elseif m == 4
        irect = 3;
      else
        done_tracing = true;
      end
    elseif irect == 2
      [jx, ix, m] = ...
        gs_trace_step(jx, ix, j1, jmid, imid, i2, y44, yb, xhmid, xvi, xhu, xvmid);
      if m == -1
        irect = 1;
      elseif m == 4
        irect = 4;
      else
        done_tracing = true;
      end
    elseif irect == 3
      [jx, ix, m] = ...
        gs_trace_step(jx, ix, jmid, j2, i1, imid, y44, yb, xhl, xvmid, xhmid, xvo);
      if m == 1
        irect = 4;
      elseif m == -4
        irect = 1;
      else
        done_tracing = true;
      end
    else % irect = 4
      [jx, ix, m] = ...
        gs_trace_step(jx, ix, jmid, j2, imid, i2, y44, yb, xhmid, xvmid, xhu, xvo);
      if m == -1
        irect = 3;
      elseif m == -4
        irect = 2;
      else
        done_tracing = true;
      end
    end
  end  
else % No place to exit
  jx = jr;
  ix = iz;
  m = 0;
end
