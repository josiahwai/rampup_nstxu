function [mutuals,lnrecta,snrecta] = mutind_distrib(acoil_data,bcoil_data, ...
					nrecta,nrectb,make_plot)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX: mutuals=mutind_distrib(acoil_data,bcoil_data, ...
%					nrecta,nrectb,make_plot)
%
%  PURPOSE:  Compute mutual inductance between 2 sets of axisymmetric 
%	conductors, each of whose cross-section is a parallelogram.  Filament
%	conductors (identified by setting dR=dZ=0 in bcoil_data) are allowed
%	for conductor set B.
%
%  INPUT:
%	acoil_data = data describing geometry of coil set A
%	bcoil_data = data describing geometry of coil set B
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
%	nrecta = rectangle partitioning for conductor set A (see header
%		documentation for bld_subelements.m)
%	nrectb = rectangle partitioning for conductor set B (see header
%		documentation for bld_subelements.m)
%	make_plot = set to 1 to plot conductors' cross-sections
%
%  OUTPUT:
%	mutuals = mutual inductance matrix, with rows corresponding to coil
%		set A, columns to coil set B (micro-Henries)
% 
%  RESTRICTIONS:
%
%  METHOD:  Each conductor is partitioned into a bunch of tiny rectangles.  All
%	mutuals are computed for these rectangles, then combined at the end.
%	Mutuals calculation is done in 2 steps:
%	(1) use filament to filament calculation for all subelements
%	(2) replace those subelement calculations which would have inaccurate
%	results, based on geometric information, with more accurate calculation.
 
%  WRITTEN BY:  Mike Walker 	ON 	7/26/99
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<5
   make_plot=0;
end

[geom_a,lnrecta,snrecta] = bld_subelements(acoil_data,nrecta,0);
disp(['lnrecta = ' int2str(lnrecta)])
disp(['snrecta = ' int2str(snrecta)])
nconductors_a = size(acoil_data,2);

[geom_b,lnrectb,snrectb] = bld_subelements(bcoil_data,nrectb,1);
disp(['lnrectb = ' int2str(lnrectb)])
disp(['snrectb = ' int2str(snrectb)])
nconductors_b = size(bcoil_data,2);

if make_plot
   figure(100)
   plot_d3d_geo
   hold on
   plot(geom_b(2,:),geom_b(1,:),'*r')
   plot(geom_a(2,:),geom_a(1,:),'xy')
   hold off
   title('distributed current element centroids: set A(yx), set B(r*)')
end

% mut1(k) represents multiplier of current in conductor k of set A to get
% flux at each subelement of conductor set B.

N = 3;		% number of cells defining near points

kk=1;
nelements_a = size(geom_a,2);
nelements_b = size(geom_b,2);
clear mut1
rb = geom_b(2,:);
zb = geom_b(1,:);

% Code below processes 1 subelement of set A at a time versus all 
% subelements of set B.

for k=1:nconductors_a
   fprintf('Processing set A conductor %d ... \n',k);
   mut1(k,:) = zeros(1,nelements_b);
   for kk=find(geom_a(3,:)==k)
      ra = geom_a(2,kk);
      za = geom_a(1,kk);
      dra = geom_a(5,kk);
      dza = geom_a(4,kk);
      mut = mutind(ra,za,rb,zb);

% find all near points
      ni = ...			% (near index)
	find((rb<ra+N*dra) & (rb>ra-N*dra) & (zb<za+N*dza) & (zb>za-N*dza));

% fix mutuals for all near points (units from mutind_near in Henries)
     if ~isempty(ni)
       mut(ni) = mutind_near(ra*ones(size(ni)),za*ones(size(ni)), ...
                 dra,dza,rb(ni),zb(ni)) * 1e6;
     end
     mut1(k,:) = mut1(k,:) + mut;
   end

% Divide current in bulk conductor (set A) between all subelements in conductor.
   mut1(k,:) = mut1(k,:)/(lnrecta(k)*snrecta(k));
end


mutuals = zeros(nconductors_a,nconductors_b);
for k=1:nconductors_b
   b_index = find(geom_b(3,:)==k);
   lenb = length(b_index);
   if lenb>1
      mutuals(:,k) = (sum(mut1(:,b_index)')/lenb)';
   else
      mutuals(:,k) = mut1(:,b_index);
   end
end

