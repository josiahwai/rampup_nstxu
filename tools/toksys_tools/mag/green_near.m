function [br,bz,dbrdr,dbrdz,dbzdr,dbzdz] = ...
		 green_near(rcurr,zcurr,rdim,zdim,rb,zb,nsplit)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX: [br,bz,dbrdr,dbrdz,dbzdr,dbzdz] = ...
%		 green_near(rcurr,zcurr,rdim,zdim,rb,zb,nsplit)
%
%  PURPOSE:  Calculate green function from current flowing in a loop
%       of rectangular cross-section to a filamentary loop which is nearby.
%       This function is meant to handle carefully the evaluation of green
%       functions for loops which are very close to one another.
%
%  INPUT:
%       rcurr= radius of centroid of current source (m)
%       zcurr= z-dimension (height) of centroid of current source (m)
%       rb = radius of point at which induced field will be measured (m)
%       zb = z-dimension (height) of point where induced field is measured (m)
%       rdim = r dimension of rectangle in which current flows (m)
%       zdim = z dimension of rectangle in which current flows (m)
%       nsplit = number of pieces to split each coordinate (z,r) of current
%                source into for calculation (option, default=[8,8])
%  Variables rcurr,zcurr, rb, and zb must be same size (scalar or vector).
%  Variables rdim, zdim must be scalar.
%
%  OUTPUT:
%       br = calculated green function for br (Tesla/Amp)
%       bz = calculated green function for bz (Tesla/Amp)
%
%  RESTRICTIONS:
%
%  METHOD:  Break up rectangle into small pieces.  Use standard method 
%  to calculate
%  values for pieces not containing measurement point. Use some other 
%  approximation (what?) for piece containing measurement point.  Add up all
%  values and divide by number of pieces (i.e. divvy up current evenly
%  into the small pieces).
 
%  WRITTEN BY:  Mike Walker     ON      8/22/97 (from IDL version)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
debug = 0;

% jal 5-23-00 green fixed so will compute br~eps for r=0 instead of NAN
if (min(rcurr) <= 0 | min(rb) < 0) 
 fprintf('green_near input ERROR: all r coordinate values must be positive\n')
end

ss1 = size(rcurr);
ss2 = size(zcurr);
ss3 = size(rb);
ss4 = size(zb);
if (any(ss1~=ss2 | ss2~=ss3 | ss3~=ss4)) 
   fprintf(' all r,z vectors must be same size\n');
   return
end
npts = length(rcurr);
%if (ss1(1) == 1 & ss1(2) == 1)
%   npts = 1;
%elseif (ss1(1) eq 1) 
%   npts = ss1(2);
%elseif (ss1(2) eq 1) 
%   npts = ss1(1);
%end
if (debug), fprintf('npts = %d\n',npts), end
br = zeros(npts,1);
bz = zeros(npts,1);
dbrdr = zeros(npts,1);
dbrdz = zeros(npts,1);
dbzdr = zeros(npts,1);
dbzdz = zeros(npts,1);

if nargin < 7, nsplit = [8,8]; end;

dz = zdim/nsplit(1);
dr = rdim/nsplit(2);

for ng = 1:npts
   rcenters = rcurr(ng)-(nsplit(2)/2)*dr + dr*[0:nsplit(2)-1] + dr/2;
   zcenters = zcurr(ng)-(nsplit(1)/2)*dz + dz*[0:nsplit(1)-1] + dz/2;

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

   for j=1:nsplit2 
      rcurrc(:,j) = rcenters(j)*ones(nsplit(1),1);
      zcurrc(:,j) = zcenters';
   end
   m=((rb(ng)+rcurrc).^2-(rb(ng)-rcurrc).^2)./ ...
			((rb(ng)+rcurrc).^2+(zcurrc-zb(ng)).^2);
   k=sqrt(m);
    
   [kk,ek]=ellipke(m);

%   rtemp = reshape(rb(ng)*onearr,n2,1);
%   ztemp = reshape(zb(ng)*onearr,n2,1);
%   [temp1, temp2] = calc_bgreens(reshape(rcurrc,n2,1), rtemp, ...
%				reshape(zcurrc,n2,1), ztemp, k, ek, kk);
   [brtemp,bztemp,dbrdrtemp,dbrdztemp,dbzdrtemp,dbzdztemp] = ...
	 calc_bgreens(rcurrc,rb(ng)*onearr,zcurrc,zb(ng)*onearr,k,ek,kk);
   brtemp    = reshape(brtemp,n2,1);
   bztemp    = reshape(bztemp,n2,1);
   dbrdrtemp = reshape(dbrdrtemp,n2,1);
   dbrdztemp = reshape(dbrdztemp,n2,1);
   dbzdrtemp = reshape(dbzdrtemp,n2,1);
   dbzdztemp = reshape(dbzdztemp,n2,1);

% If the test point is actually within the distributed current element,
% then need to replace green function with something else for current
% elements closest to the measurement point.
% (e.g. when source and msmnt grid pt are the same, 4 discretized pts lie
%  at exactly the same distance from the msmt point)
% Right now the "something else" is 0.  What is a better choice?

   if ( (rb(ng) < rcurr(ng)+rdim/2)  & ...
        (rb(ng) > rcurr(ng)-rdim/2)  & ...
        (zb(ng) < zcurr(ng)+zdim/2)  & ...
        (zb(ng) > zcurr(ng)-zdim/2) ) 

      distance = sqrt((rcurrc-rb(ng)).^2 + (zcurrc-zb(ng)).^2);
      distance = reshape(distance,n2,1);
      mindistance = min(distance);
      minindex = find(abs(distance-mindistance)<1e-6);

      angle = atan2(zb(ng) - zcurrc,rb(ng)-rcurrc);
      angle = reshape(angle,n2,1);
      minangle = angle(minindex);

% use circular approximation for "close" cell
	a = sqrt(dr*dz/pi);
	mu0 = 4*pi*1e-7;
%	magnitude_B = mu0 * mindistance/(2*pi*a^2)
 	magnitude_B = 0;

%help,minindex,rsame,zsame,brtemp,bztemp,mindistance,minangle,magnitude_B

 	brtemp(minindex) = magnitude_B * cos(minangle);
 	bztemp(minindex) = magnitude_B * sin(minangle);

% WHAT SHOULD THESE BE????

        dbrdrtemp(minindex) = 0;
        dbrdztemp(minindex) = 0;
        dbzdrtemp(minindex) = 0;
        dbzdztemp(minindex) = 0;
   end

%help,a,mu0,magnitude_B,minindex,minangle
%fprintf('brtemp = ',brtemp
%fprintf('bztemp = ',bztemp

   br(ng) = mean(brtemp);
   bz(ng) = mean(bztemp);
   dbrdr(ng) = mean(dbrdrtemp);
   dbrdz(ng) = mean(dbrdztemp);
   dbzdr(ng) = mean(dbzdrtemp);
   dbzdz(ng) = mean(dbzdztemp);
end
