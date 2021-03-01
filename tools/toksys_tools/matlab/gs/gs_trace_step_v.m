function [jx, ix, m, vjx, vix] = gs_trace_step_v(jr, iz, j1, j2, i1, i2, y44, yb, vj, vi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   [jx,ix,m,vjx,vix] = gs_trace_step_v(jr,iz,j1,j2,i1,i2,y44,yb,vj,vi);
%
%  PURPOSE: Find where the contour y = yb exits the rectangle j1, j2, i1, i2
%
%  INPUTS: jr, iz, entry point for the contour y = yb (both values >=0 & <=1)
%          j1, j2, min and max r of the rectangle (both values >=0 & <=1)
%          i1, i2, min and max z of the rectangle (both values >=0 & <=1)
%          y44, values of y for jr,iz = -1, 0, 1, 2
%          yb, value for which contour is sought
%
%  OUTPUTS:  jx, ix, r, z of exit point
%            m, where exit is, 0=nowhere, -1=bottom, +1=top, -4=in, +4=out
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
%  METHOD:  If there is only one place to exit then pick it, otherwise
%           call this function recursively to step through in 2 steps

%
%  WRITTEN BY:  Anders Welander  ON	2/26/14
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2; % For interpolation

if exist('vj') == 1
  vjx = vj;
else
  vjx = [];
end
if exist('vi') == 1
  vix = vi;
else
  vix = [];
end


psi4 = [1 j1 j1^2 j1^3]*mx*y44';
xvi = gs_solve_hermite_cubic(psi4, yb)';
psi4 = [1 j2 j2^2 j2^3]*mx*y44';
xvo = gs_solve_hermite_cubic(psi4, yb)';
psi4 = [1 i1 i1^2 i1^3]*mx*y44;
xhl = gs_solve_hermite_cubic(psi4, yb)';
psi4 = [1 i2 i2^2 i2^3]*mx*y44;
xhu = gs_solve_hermite_cubic(psi4, yb)';

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

if sum(r12 & k12 & f12) == 1 % There is only one way to exit
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
  vjx(end+1) = jx;
  vix(end+1) = ix;
elseif sum(r12 & k12 & f12) > 1 % More than one way to exit, that's complicated
  % Divide this rectangle into 4 pieces number 1,2,3,4 and step through them
  imid = (i1+i2)/2;
  jmid = (j1+j2)/2;
  jx = jr;
  ix = iz;
  irect = 1 + (iz >= imid) + 2*(jr >= jmid);
  done_tracing = false;
  while ~done_tracing
    if irect == 1
      [jx, ix, m, vjx, vix] = gs_trace_step_v(jx, ix, j1, jmid, i1, imid, y44, yb, vjx, vix);
      if m == 1
        irect = 2;
      elseif m == 4
        irect = 3;
      else
        done_tracing = true;
      end
    elseif irect == 2
      [jx, ix, m, vjx, vix] = gs_trace_step_v(jx, ix, j1, jmid, imid, i2, y44, yb, vjx, vix);
      if m == -1
        irect = 1;
      elseif m == 4
        irect = 4;
      else
        done_tracing = true;
      end
    elseif irect == 3
      [jx, ix, m, vjx, vix] = gs_trace_step_v(jx, ix, jmid, j2, i1, imid, y44, yb, vjx, vix);
      if m == 1
        irect = 4;
      elseif m == -4
        irect = 1;
      else
        done_tracing = true;
      end
    else % irect = 4
      [jx, ix, m, vjx, vix] = gs_trace_step_v(jx, ix, jmid, j2, imid, i2, y44, yb, vjx, vix);
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
