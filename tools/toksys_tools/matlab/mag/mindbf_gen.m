function [mut,br,bz,badindex]= mindbf_gen(ra,rb,za,zb,rdim,zdim,ndist,nsplit)
 %
%  SYNTAX: [mut,br,bz,badindex] = mindbf_gen(ra,rb,za,zb,rdim,zdim,ndist,nsplit)
%
%  PURPOSE:  Calculate mutual inductance, B_r, and B_z from a current source
%   "a" which is a coil carrying distributed current to a measurement point "b".
%   Note that the fields are calculated at point b due to current at point a.
%   Green functions are Henries and Tesla/Amp.  This routine is a generalized
%   version of Dave Humphreys' mindbf.pro which handles carefully all the
%   problems that I could think of having to do with a distributed current
%   source and with points a and b close together.
%
%  INPUT:
%    ra,za = r and z coordinates of current source(s) (meters)
%    rb,zb = r and z coordinates of "measurement point(s)" (meters)
%    rdim = r dimension of rectangular distributed current element source(s) (m)
%    zdim = z dimension of rectangular distributed current element source(s) (m)
%    ndist = number of multiples of rdim in r dimension and zdim in z 
%		dimension away from current element that a measurement point
%		point should be in order to use simple filament-to-filament
%		approximation (optional, default=3)
%    nsplit = number of pieces to split each coordinate (z,r) of current
%                source into for "near" calculation (optional, default=[8,8])
%  (Allows input of vectors of ra,za,rb,zb if size of a and b objects equal or one is scalar.)
%
%  OUTPUT:
%   mut = mutual inductance between 2 loops at locations (ra,za), (rb,zb)
%               (Henries)
%   br = multiplier of current at pt(s) a to produce radial field at pt(s) b 
%               (Tesla/Amp)
%   bz = multiplier of current at pt(s) a to produce vertical field at pt(s) b
%               (Tesla/Amp)
%   badindex = if >= 0, index in input vector(s) for which elliptic integral 
%		calculator could not give good value.  This means values of
%		mut(badindex), br(badindex), and bz(badindex) are NOT VALID. 
%  Note: I think I've fixed it so that a bad value is never returned, but I
%  want to leave this in for awhile so I can check if it's ever >=0.
%
%  RESTRICTIONS:
%  NOTE that this ASSUMES there will be AT MOST ONE instance of points
%  from the a and b vectors which are too close together to calculate responses
%  correctly.
%
%  If vector input, all vectors ra, rb, za, and zb must be same size.  In this
%  case, output corresponds to pairings of elements in sets a and b.
%  For example, mut is calculated between (ra(k),za(k)) and (rb(k),zb(k)), but
%  not between (ra(k),za(k)) and (rb(j),zb(j)) for j different from k.
%
%  NOTE: A number of "divide by 0" messages may be generated on first pass
%  of some calculations.  These are fixed later on by recomputation of 
%  data for "too close" points.
 
%  METHOD:  
%
%  WRITTEN BY:  Dave Humphreys	ON 	9/2/91
%  MODIFIED BY: Mike Walker	ON 	10/30/95 (generalized routine)
%  MODIFIED BY: Mike Walker	ON 	7/31/96 (convert from idl to matlab)
%
%  NOTE that this should be maintained in parallel with the IDL version which
%  is in PCS code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
debug=0;
if nargin < 8, nsplit=[8,8];, end;
if nargin < 7, ndist=3;, end;

if (min(ra) <= 0 | min(rb) <= 0) 
   print, 'mindbf_gen input ERROR: all r coordinate values must be positive'
end

    twopir = 2.*pi*rb;
    m=((rb+ra).^2 - (rb-ra).^2)./((rb+ra).^2 + (za-zb).^2);
    k=sqrt(m);
    [kk,ek] = ellipke(m);

% Handle the exceptions where elliptic integrals couldn't calculate correctly.
% The incorrect value which is assigned here is replaced with a correct value
% in the "too close" calculations below.  This index is still passed out of
% this routine for debugging purposes.
    badindex = find(kk == 1e+38 | isinf(kk));
    if ~isempty(badindex) 
        k(badindex) = .5;		% just to get it to run cleanly
        kk(badindex) = 5;		% just to get it to run cleanly
        ek(badindex) = 2;		% just to get it to run cleanly
 	%print,'mindbf: Calculated value for index ',badindex, ' is INVALID!'
    end

    mut=8.*pi*sqrt(ra.*rb).*(1.0e-7).*((1.-.5*m).*kk-ek)./k;

%print,k,ek,kk,mut

% B-field calculations:
%    ss=size(ra);
%    sss=size(rb);
%    iscal=((ss(0) == 0) & (sss(0) == 0));

% Use elliptic integrals for all points:
    [br,bz] = calc_bgreens(ra, rb, za, zb, k, ek, kk);

% fix calculations for measurement points which are "too close" to current

   check = (rb < ra+ndist*rdim) & (rb > ra-ndist*rdim) & ...
           (zb < za+ndist*zdim) & (zb > za-ndist*zdim);
   mi = find(check);
   if ~isempty(mi) 
      if(length(rb)==1)
         rb1 = rb*ones(size(ra));
         zb1 = zb*ones(size(ra));
      else
         rb1 = rb;
         zb1 = zb;
      end
      if(length(rdim)==1)
         rdim1 = rdim*ones(size(ra));
      else
         rdim1 = rdim;
      end
      if(length(zdim)==1)
         zdim1 = zdim*ones(size(ra));
      else
         zdim1 = zdim;
      end


      if(debug)
	fprintf('close points at mi =  ')
        for j=1:length(mi)
           fprintf('%d ',mi(j))
        end
        fprintf('\n')
	fprintf('ra(mi),rb(mi),za(mi),zb(mi) \n')
        for j=1:length(mi)
           fprintf('%d %d %d %d\n',ra(mi(j)), rb1(mi(j)),za(mi(j)),zb1(mi(j)))
        end
      end
	ss = size(mi);
	npts = ss(1);

	rab = ra(mi);
	zab = za(mi);
	rbb = rb1(mi);
	zbb = zb1(mi);

     	mut(mi) = mutind_near(rab,zab,rdim1(mi),zdim1(mi),rbb,zbb,nsplit);
%  	[brmi,bzmi] = green_near(rab,zab,rdim,zdim,rbb,zbb,nsplit);
        ri = rab-rdim1(mi)/2;
        ro = rab+rdim1(mi)/2;
        zl = zab-zdim1(mi)/2;
        zu = zab+zdim1(mi)/2;
        [brmi,bzmi,fl] = fine1(ri,ro,zl,zu,ones(size(ri)),rbb,zbb);
        brmi = 4*pi*1e-7 * brmi;
        bzmi = 4*pi*1e-7 * bzmi;

% SPECIAL CASE: If zab==zbb, then Br=0 by symmetry.
idx=find(zab==zbb); brmi(idx)=0;

        br(mi) = brmi;
        bz(mi) = bzmi;

% This is only for comparison with EFIT values, where self inductance for
% a rectangular current element is calculated.  For control points, 
% measurement point is always a filament.
if 0 
        ss = size(mi)
        for i=0,ss(1)-1 
	   rab = ra(mi(i))
	   zab = za(mi(i))
	   rbb = rb(mi(i))
	   zbb = zb(mi(i))
           if (rab == rbb & zab == zbb) 
              print,mi(i)
              print,rdim,zdim
              print,'rab=',rab
              print,'zab=',zab
              print,'rbb=',rbb
              print,'zbb=',zbb
              print,'mut(mi(i))=',mut(mi(i))

              mut(mi(i)) = rectl(rab,zdim,rdim)
              print,'mut(mi(i))=',mut(mi(i))
           end
        end

end


    end

              
