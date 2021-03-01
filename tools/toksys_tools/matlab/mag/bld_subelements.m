function [geometry,lnrect,snrect] =  ...
			bld_subelements(cond_data,nrect,allow_filaments)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  [geometry,lnrect,snrect] = ...
%			bld_subelements(cond_data,nrect,allow_filaments)
%
%  PURPOSE:  Split conductors into rectangular subelements.
%
%  INPUT:
%	cond_data = data describing geometry of conductor 
%	   Format of this data has each column = [Z;R;dZ;dR;AC;AC2], where:
%		Z   = z centroids of conductors (m)
%		R   = r centroids of conductors (m)
%		dZ  = height of conductors (m)
%		dR  = width of conductors (m)
%		AC  = angle between coil and horizontal (degrees)
%		AC2 = angle between coil and horizontal (degrees)
%	At most one of AC, AC2 is nonzero and indicates type of parallelogram:
%                                xxxx           
%                            xxxx   x
%                        xxxx       x
%                    xxxx           x       ---                 xxxxxxxxxxxxx   
%                xxxx               x        |                xx         xx     
%       ---  xxxx                   x        |              xx         xx       
%        |   x          X        xxxx      dZ             xx         xx   
%        |   x        (R,Z)  xxxx            |          xx   X     xx           
%      dZ    x           xxxx                |        xx   (R,Z) xx             
%        |   x       xxxx                    |      xx         xx               
%        |   x   xxxx                        |     xx         xx                
%       _|_  xxxx   ^                       --- xxxxxxxxxxxxx    ^              
%                 AC )                                        AC2  )            
%            |----------dR----------|           |---dR------|                   
%
%       nrect = number of rectangles to partition conductors. Interpretation
%		depends on size of this object:
%		scalar: nrect = number of rectangles across largest of conductor
%			minimum dimensions (thickest conductor).  All 
%			partitioning automatically computed from this.  (Gives
%			roughly constant subrectangle size over all conductors.)
%		length vector == 2: nrect specifies FIXED nshort x nlong
%			partitioning, nshort refers to min(dR,dZ)
%		2 x nconductors or nconductors x 2:  nrect specifies FIXED
%			values of (nshort, nlong) for each conductor
%		nconductors length vector: nrect specifies fixed nshort for
%			each conductor, nlong computed by dR/dZ scaling
%	allow_filaments = set to 1 to allow filament (dR=dZ=0) coils, otherwise
%		an error message will be printed if this is detected
%
%  OUTPUT:
%	geometry = 7 x n matrix describing all subelements of all conductors
% 		Each column has the following information:
%   		[z of current element; r of current element; 
%		conductor number; dz of subelement; dr of subelement;
%		row number of element within conductor;
%		column number of element within conductor]
%	lnrect = number of rectangles in longer dimension, for each conductor
  
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	7/26/99
%
%  @(#)bld_subelements.m	1.4 02/12/09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
debug = 0;
dZ  = cond_data(3,:);
dR  = cond_data(4,:);
if min(dZ)<0 | min(dR)<0, wait('ERROR bld_subelements: Logic needs dR>0, dZ>0'), end
nconductors = size(cond_data,2);

if length(nrect(:))==1		% calculate both short and long number rects
   snrect0 = nrect;
   max_thickness = max(min([dR; dZ]));	% "thickness" is shortest dimension
   if (max_thickness==0 && max(max([dR; dZ]))>0)   % all flat geometries
      max_thickness = max(max([dR; dZ])); % use larger dimension
   end   
   if max_thickness~=0
      for k=1:nconductors
         if dR(k)>dZ(k)
            lnrect(k) = max([1, round(dR(k)/max_thickness * snrect0)]);
            snrect(k) = max([1, round(dZ(k)/max_thickness * snrect0)]);
         else
            lnrect(k) = max([1, round(dZ(k)/max_thickness * snrect0)]);
            snrect(k) = max([1, round(dR(k)/max_thickness * snrect0)]);
         end;
      end
   else  % all filaments
      lnrect = ones(1,nconductors);
      snrect = ones(1,nconductors);
   end
elseif length(nrect(:))==2	% all rects specified to be the same
   snrect0 = nrect(1);
   lnrect0 = nrect(2);
   for k=1:nconductors
      lnrect(k) = lnrect0;
      snrect(k) = snrect0;
   end;
else
   sz = size(nrect);
   if min(sz)==1		% vector input ==> short number rects specified
      if(sz(1)>1)	% need for snrect to be row vector
         snrect = nrect';
      else
         snrect = nrect;
      end
      for k=1:nconductors
         if(dR(k)==0 | dZ(k)==0)
            if(~allow_filaments)
               wait('ERROR bld_subelements: non-filament dR or dZ =0')
               return;
            end
            if(~(dZ(k)==0 & dR(k)==0))
               wait(['ERROR bld_subelements: only one of dR,dZ=0, ' ...
			'conductor ' int2str(k)])
               return;
            else
               lnrect(k)=1;
               snrect(k)=1;
            end
         elseif dR(k)>dZ(k)
            lnrect(k) = round(dR(k)/dZ(k) * snrect(k));
         else
            lnrect(k) = round(dZ(k)/dR(k) * snrect(k));
         end;
      end;
   else				% 2 x n or n x 2 data ==> all rects specified
      [mm,mi]=min(sz);
      if mi==1			% 2 x n
	 snrect = nrect(1,:);
	 lnrect = nrect(2,:);
      else			% n x 2
	 snrect = nrect(:,1)';
	 lnrect = nrect(:,2)';
      end
   end
end

% We always index R and Z so that R index changes faster.  The number of 
% elements and their indices are the same with or without AC or AC2.
% Locations of subelement (r,z) centers will change however.

geometry = [];
for k=1:size(cond_data,2)
   Z = cond_data(1,k);
   R = cond_data(2,k);
   AC = cond_data(5,k);
   AC2 = cond_data(6,k);
   dZ  = cond_data(3,k);
   dR  = cond_data(4,k);
   s_array = [0.5:1:snrect(k)]/snrect(k);
   s_index = [1:snrect(k)];
   l_array = [0.5:1:lnrect(k)]/lnrect(k);
   l_index = [1:lnrect(k)];

   if(dZ~=0 | dR~=0)	% if either dZ or dR ~= 0, then split into subelements
      if AC~=0
         dZtot = cond_data(3,k) + dR*tan(AC*pi/180);
         dRtot = dR;
      elseif AC2~=0
         dRtot = cond_data(4,k) + dZ/tan(AC2*pi/180);
         dZtot = dZ;
      else
         dZtot = dZ;
         dRtot = dR;
      end
      if(dR > dZ)
         nZ = snrect(k);
         nR = lnrect(k);
         dR_array = dR*l_array;
         dZ_array = dZ*s_array;
         r_index = l_index;
         z_index = s_index;
      else
         nR = snrect(k);
         nZ = lnrect(k);
         dR_array = dR*s_array;
         dZ_array = dZ*l_array;
         r_index = s_index;
         z_index = l_index;
      end
      r_array = R - dRtot/2 + dR_array;
      z_array = Z - dZtot/2 + dZ_array;
      drsub = dR/nR;
      dzsub = dZ/nZ;

      if(AC~=0) 	% (z_array will vary but r_array will remain fixed)
         rtemp = ones(nZ,nR)*diag(r_array);
         iztemp = diag(z_index)*ones(nZ,nR);
         irtemp = ones(nZ,nR)*diag(r_index);

         z_array = z_array + dR_array(1)*tan(AC*pi/180);

         ztemp = z_array';
         zshift =  dR/nR*tan(AC*pi/180);
         for l=2:nR
            ztemp(:,l) = z_array' + (l-1)*zshift;
         end

         geometry = [geometry  ...
	  [ztemp(:)';rtemp(:)';k*ones(1,nR*nZ);dzsub*ones(1,nR*nZ); ...
			drsub*ones(1,nR*nZ);iztemp(:)';irtemp(:)']];

      elseif(AC2~=0) 	% (r_array will vary but z_array will remain fixed)
         ztemp = diag(z_array)*ones(nZ,nR);
         iztemp = diag(z_index)*ones(nZ,nR);
         irtemp = ones(nZ,nR)*diag(r_index);

         r_array = r_array + dZ_array(1)/tan(AC2*pi/180);

         rtemp = r_array;
         rshift =  dZ/nZ/tan(AC2*pi/180);
         for l=2:nZ
            rtemp(l,:) = r_array + (l-1)*rshift;
         end

         geometry = [geometry  ...
	  [ztemp(:)';rtemp(:)';k*ones(1,nR*nZ);dzsub*ones(1,nR*nZ); ...
			drsub*ones(1,nR*nZ);iztemp(:)';irtemp(:)']];

      else				% just a regular rectangle
         ztemp = diag(z_array)*ones(nZ,nR);
         rtemp = ones(nZ,nR)*diag(r_array);
         iztemp = diag(z_index)*ones(nZ,nR);
         irtemp = ones(nZ,nR)*diag(r_index);
         geometry = [geometry  ...
	  [ztemp(:)';rtemp(:)';k*ones(1,nR*nZ);dzsub*ones(1,nR*nZ); ...
			drsub*ones(1,nR*nZ);iztemp(:)';irtemp(:)']];

      end
   else 		% if dZ=dR=0, can't split into subelements

      if(allow_filaments)  % if filaments are allowed, add to geometry
         geometry = [geometry [Z;R;k;dZ;dR;1;1]];

      else		% otherwise, this is an error
	 fprintf('ERROR in bld_subelements: dR=dZ=0 for conductor %d\n',k);
      end
   end
   if(debug)
      figure(10),clf
      axis equal
      plot_box(cond_data(2,:),cond_data(1,:),cond_data(4,:),cond_data(3,:),'g', ...
                cond_data(5,:),cond_data(6,:))
      hold on
      rs = geometry(2,:);
      zs = geometry(1,:);
      delrs = geometry(5,:);
      delzs = geometry(4,:);
      plot_box(rs,zs,delrs,delzs,'b',zeros(size(rs)),zeros(size(rs)));
      hold off
      wait
   end
end
