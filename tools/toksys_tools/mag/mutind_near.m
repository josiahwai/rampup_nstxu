function mutual = mutind_near(rcurr,zcurr,rdim,zdim,rpsi,zpsi,nsplit)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  mutual = mutind_near(rcurr,zcurr,rdim,zdim,rpsi,zpsi)
%
%  PURPOSE:  Calculate mutual inductance from current flowing in a loop
%	of rectangular cross-section to a filamentary loop which is nearby.
%	This function is meant to handle carefully the evaluation of mutual
%	inductances for loops which are very close to one another.
%
%  INPUT:
%	rcurr= radius of centroid of current source (m)
%	zcurr= z-dimension (height) of centroid of current source (m)
%	rpsi = radius of point at which induced flux will be measured (m)
%	zpsi = z-dimension (height) of point where induced flux is measured (m)
%	rdim = r dimension of rectangle in which current flows (m)
%	zdim = z dimension of rectangle in which current flows (m)
%	nsplit = number of pieces to split each coordinate (z,r) of current
%		 source into for calculation (optional, default=[8,8])
%  Variables rcurr,zcurr, rpsi, and zpsi must be same size (scalar or vector).
%  Variables rdim, zdim must be scalar or same size as rcurr.
%
%  OUTPUT:
%	mut  = calculated mutual inductance (Henries) from the current source
%		at each (rcurr(i),zcurr(i)) to measurement at (rpsi(i),zpsi(i))

%  RESTRICTIONS:
%
%  METHOD:  Break up rectangle into small pieces.  Use mutind to calculate
%  mutuals for pieces not containing measurement point. Use 
%  self-inductance for piece containing measurement point.  Add up all 
%  inductances and divide by number of pieces (i.e. divvy up current evenly
%  into the small pieces).
%
%  WRITTEN BY:  Mike Walker 	ON 	10/27/95
%  UPDATED BY:  Mike Walker	ON	8/1/96 (convert from IDL to MATLAB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%print,rcurr
%print,zcurr
%print,'rdim,zdim = ',rdim,zdim
%print,rpsi
%print,zpsi

ss1 = size(rcurr);
ss2 = size(zcurr);
ss3 = size(rpsi);
ss4 = size(zpsi);
if (ss1 ~= ss2 | ss2 ~= ss3 | ss3 ~= ss4)
   print, ' all r,z vectors must be same size'
   return,-1;
end
if (length(ss1) == 1) 			% single (scalar) current center
   npts = 1;
   mutual = zeros(npts);
elseif (length(ss1) > 1) 	% vector of current centers
   npts = max(ss1);
   mutual = zeros(ss1);
end
%end else if (ss1(0) == 2) 
%   npts = ss1(1)*ss1(2)
%   mutual = zeros(ss1(1),ss1(2))
%end

if nargin < 7, nsplit = [8,8]; end;

if(length(rdim)==1)
   rdim1 = rdim*ones(size(rcurr));
else
   rdim1 = rdim;
end
if(length(zdim)==1)
   zdim1 = zdim*ones(size(zcurr));
else
   zdim1 = zdim;
end
dz = zdim1/nsplit(1);
dr = rdim1/nsplit(2);

for ng = 1:npts

%print,'ng = ',ng
   rcenters = rcurr(ng)-(nsplit(2)/2)*dr(ng)+dr(ng)*[0:nsplit(2)-1]+dr(ng)/2;
   zcenters = zcurr(ng)-(nsplit(1)/2)*dz(ng)+dz(ng)*[0:nsplit(1)-1]+dz(ng)/2;

% handle current elements right at the R=0 axis:
   idx = find(rcenters<0);
   nidx = length(idx);
   if(nidx)
      rcenters = rcenters(idx(end)+1:end);
      nsplit2 = nsplit(2)-nidx;
   else
      nsplit2 = nsplit(2);
   end
   n2 = nsplit(1)*nsplit2;
   rcurrc = zeros(nsplit(1),nsplit2);
   zcurrc = zeros(nsplit(1),nsplit2);
   mtemp = zeros(n2,1);
   onearr = ones(nsplit(1),nsplit2);

   for k=1:nsplit2
      rcurrc(:,k) = rcenters(k)*ones(nsplit(1),1);
      zcurrc(:,k) = zcenters';
   end
%print,rcurrc,zcurrc

   rtemp = reshape(rpsi(ng)*onearr,n2,1);
   ztemp = reshape(zpsi(ng)*onearr,n2,1);
   mtemp = mutind(reshape(rcurrc,n2,1),reshape(zcurrc,n2,1),rtemp,ztemp);

% If the test point is actually within the distributed current element, 
% then need to replace mutual inductance with a self inductance for current
% elements closest to the measurement point.
% (e.g. when source and msmnt grid pt are the same, 4 discretized pts lie
%  at exactly the same distance from the msmt point)

   if ( (rpsi(ng) < rcurr(ng)+rdim1(ng)/2)   &          ...
	(rpsi(ng) > rcurr(ng)-rdim1(ng)/2)   &          ...
	(zpsi(ng) < zcurr(ng)+zdim1(ng)/2)   &          ...
	(zpsi(ng) > zcurr(ng)-zdim1(ng)/2) ) 

      distance = sqrt((rcurrc-rpsi(ng)).^2 + (zcurrc-zpsi(ng)).^2);
      distance = reshape(distance,n2,1);
      mindistance = min(distance);
      minindex = find(abs(distance-mindistance)<1e-6);

      rr = reshape(rcurrc,n2,1);
      mtemp(minindex) = rectl(rr(minindex),dz(ng),dr(ng));

	%help,mtemp
	%print,'mtemp = ',mtemp
	%print,'minindex = ',minindex
	%print,'mtemp(minindex) = ',mtemp(minindex)
	%print,ones(n2) / nsplit^2
   end

   mutual(ng) = mean(mtemp)*1e-6;  % both mutind and rectl give uH

end

