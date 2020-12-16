function [x, y, f, tA, tB] = line_intersection(xA, yA, xB, yB);
%
%  USAGE:   [x, y, f, tA, tB] = line_intersection(xA, yA, xB, yB);
%
%  PURPOSE: Return point of intersection between lines A and B
%
%  INPUTS:  xA, yA = x, y for two points on line A, sizes [1,2]
%           xB, yB = x, y for two points on line B, sizes [1,2]
%           If called with multiple line-pairs (n), sizes are [n,2]
%           For this case any input can still be size [1,2]
%
%  OUTPUTS: x, y = point where lines A and B intersect
%           f = flag regarding segment between specified points
%             0 = A & B are parallel and never intersect
%             1 = A & B are the same line but the segments have no overlap
%             2 = A & B intersect but not between specified points
%             3 = Line A intersected between specified points
%             4 = Line B intersected between specified points
%             5 = Both lines intersect between specified points
%             6 = A & B are the same line and segments partly overlap (not yet implemented)
%             7 = Points xA, yA are the same as xB, yB (not yet implemented)
%           tA, fraction given by: x = xA(:,1) + tA.*(xA(:,2)-xA(:,1))
%           tB, fraction given by: x = xB(:,1) + tB.*(xB(:,2)-xB(:,1))
	
%  VERSION @(#)line_intersection.m	1.2 02/08/15
%
%  WRITTEN BY:  Anders Welander  ON	7/1/14
%
%  MODIFICATION HISTORY: ASW 2/1/15 Using solved 2x2 matrix and vectors
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
  error('4 inputs are required, all with 2 columns') 
end

if size(xA,2) ~= 2 | size(yA,2) ~= 2 | size(xB,2) ~= 2 | size(yB,2) ~= 2
  error('All inputs must have 2 columns') 
end

n = max([size(xA,1), size(yA,1), size(xB,1), size(yB,1)]);
if size(xA,1) > 1 & size(xA,1) ~= n | ...
   size(yA,1) > 1 & size(yA,1) ~= n | ...
   size(xB,1) > 1 & size(xB,1) ~= n | ...
   size(yB,1) > 1 & size(yB,1) ~= n
  error('All inputs must have either 1 row or same larger number of rows') 
end

% inv([a b; c d]) = [d -b;-c a]/(a*d-b*c)
a = xA(:,2)-xA(:,1);
b = xB(:,1)-xB(:,2);
c = yA(:,2)-yA(:,1);
d = yB(:,1)-yB(:,2);
e = a.*d-b.*c;
r1 = xB(:,1)-xA(:,1);
r2 = yB(:,1)-yA(:,1);
tA = ( d.*r1-b.*r2)./e;
tB = (-c.*r1+a.*r2)./e;
x = xA(:,1) + tA.*a;
y = yA(:,1) + tA.*c;
f = 2 + (tA >= 0 & tA <= 1) +  2*(tB >= 0 & tB <= 1);
