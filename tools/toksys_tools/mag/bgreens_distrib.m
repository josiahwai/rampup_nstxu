function [brgreens,bzgreens,bpgreens,lnrectc,snrectc] = ...
		bgreens_distrib(coil_data,bprobe_data,nrect,verbose)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX: [brgreens,bzgreens,bpgreens] = ...
%		bgreens_distrib(coil_data,bprobe_data,nrect,verbose)
%
%  PURPOSE:  Compute Greens functions from axisymmetric current carrying 
%	conductors, each of whose cross-section is a parallelogram, to Br,
%	Bz, and angled probe measurement at a set of specified points.  
%	(Use mutind_distrib with bcoil_data = filaments to compute Greens 
%	functions for flux.)  Can be used to compute Br,Bz Greens functions
%	for arbitrary points in space.
%
%  INPUT:
%	coil_data = data describing geometry of coil set
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
%	bprobe_data  = Bprobe data, 5 x number of probes (or any msmnt pts): 
%		bpdata(1:5,:)=[Z;R;Ang;length(m);Ip coeff]
%		(5th row is ignored, so 4 x number probes also works)
%	nrect = rectangle partitioning for conductor set (see header
%		documentation for bld_subelements.m)
%	verbose = print/display verbose info.  0 = no print/display.  Numbers>0
%		give increasingly more information/displays.
%
%  OUTPUT:
%	brgreens = Br Greens function matrix, with rows corresponding to probes,
%		columns to conductors (Tesla/MA)
%	bzgreens = Bz Greens function matrix, with rows corresponding to probes,
%		columns to conductors (Tesla/MA)
%	bpgreens = Greens function matrix to field measured by bprobes, with
%		rows corresponding to probes, columns to conductors (Tesla/MA)
% 
%  RESTRICTIONS:
%
%  METHOD:  Each conductor is partitioned into a bunch of tiny rectangles.  All
%	bgreens are computed for these rectangles, then combined at the end.
%	Calculation is done in 2 steps:
%	(1) use filament to filament calculation for all subelements
%	(2) replace those subelement calculations which would have inaccurate
%	results, based on geometric information, with more accurate calculation.
 
%  WRITTEN BY:  Mike Walker 	ON 	7/26/99
%  jal10-3-03 remove one line commented below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4
   verbose=0;
end

[geom_cond,lnrectc,snrectc] = bld_subelements(coil_data,nrect,0);
if(verbose>1)
   disp(['lnrectc = ' int2str(lnrectc)])
   disp(['snrectc = ' int2str(snrectc)])
end
nconductors = size(coil_data,2);

% convert bpdata format to conductor format
bpdata = bprobe_data;
nprobes = size(bprobe_data,2);

cosangles = cos(bprobe_data(3,:)*pi/180);
sidx = find(abs(cosangles)<0.5);
all_index = 1:nprobes;
cidx = setdiff(all_index,sidx);
bpdata(5,:) = zeros(1,nprobes);
bpdata(6,:) = zeros(1,nprobes);

% Note conversion of probe to trapazoid format increases dimension by one
% this is reason for dZ,dR=2cm for second dimension below (jal)
bpdata(3,cidx) = 0.02*ones(1,length(cidx));	% make all dZ=2cm
bpdata(4,cidx) = abs(bprobe_data(4,cidx).*cos(bprobe_data(3,cidx)*pi/180));
bpdata(5,cidx) = bprobe_data(3,cidx);		% make all AC=angles

bpdata(3,sidx) = abs(bprobe_data(4,sidx).*sin(bprobe_data(3,sidx)*pi/180));
bpdata(4,sidx) = 0.02*ones(1,length(sidx));	% make all dR=2cm
bpdata(6,sidx) = bprobe_data(3,sidx);		% make all AC2=angles


[geom_probe,lnrectp,snrectp] = bld_subelements(bpdata,1,1);
if(verbose>1)
   disp(['lnrectp = ' int2str(lnrectp)])
   disp(['snrectp = ' int2str(snrectp)])
end

if verbose>2
   figure(100)
   plot_d3d_geo
   hold on
   plot(geom_cond(2,:),geom_cond(1,:),'xy')
   for k=1:nprobes
      index = find(geom_probe(3,:)==k);
      plot(geom_probe(2,index),geom_probe(1,index),'*r')
      r = mean(geom_probe(2,index));
      z = mean(geom_probe(1,index));
      axis([r-.2 r+.2 z-.2 z+.2])
      wait
   end
   hold off
   title('distributed current element centroids: set A(yx), set B(r*)')
end

% green(k) represents multiplier of current in conductor k of set A to get
% flux at each subelement of conductor set B.

N = 3;		% number of cells defining near points

kk=1;
nelements_p = size(geom_probe,2);
clear brgreen1 bzgreen1 bpgreen1
rb = geom_probe(2,:);
zb = geom_probe(1,:);

% Code below processes 1 subelement of conductor at a time versus all 
% measurement points.

for k=1:nconductors
   if(verbose>0)
      fprintf('Processing conductor %d ... \n',k);
   end
   brgreen1(k,:) = zeros(1,nelements_p);
   bzgreen1(k,:) = zeros(1,nelements_p);
   for kk=find(geom_cond(3,:)==k)
      ra = geom_cond(2,kk);
      za = geom_cond(1,kk);
      dra = geom_cond(5,kk);
      dza = geom_cond(4,kk);
      [brg,bzg,dbrdrg,dbrdzg,dbzdrg,dbzdzg] = calc_bgreens(ra,rb,za,zb);
      brg = brg*1e6;		% Tesla/MA
      bzg = bzg*1e6;		% Tesla/MA
      dbrdrg = dbrdrg*1e6;	% Tesla/MA/s
      dbrdzg = dbrdzg*1e6;	% Tesla/MA/s
      dbzdrg = dbzdrg*1e6;	% Tesla/MA/s
      dbzdzg = dbzdzg*1e6;	% Tesla/MA/s

% find all near points
      ni = ...			% (near index)
	find((rb<ra+N*dra) & (rb>ra-N*dra) & (zb<za+N*dza) & (zb>za-N*dza));

% fix bgreens for all near points (units from mutind_near in Henries)
      [brtemp,bztemp,dbrdrtemp,dbrdztemp,dbzdrtemp,dbzdztemp] = ...
	 green_near(ra*ones(size(ni)),za*ones(size(ni)),dra,dza,rb(ni),zb(ni));
      brg(ni)    = brtemp*1e6;		% Tesla/MA
      bzg(ni)    = bztemp*1e6;		% Tesla/MA
      dbrdrg(ni) = dbrdrtemp*1e6;	% Tesla/MA
      dbrdzg(ni) = dbrdztemp*1e6;	% Tesla/MA
      dbzdrg(ni) = dbzdrtemp*1e6;	% Tesla/MA
      dbzdzg(ni) = dbzdztemp*1e6;	% Tesla/MA

      brgreen1(k,:) = brgreen1(k,:) + brg;
      bzgreen1(k,:) = bzgreen1(k,:) + bzg;
   end

% Divide current in bulk conductor set between all subelements in conductor.
   brgreen1(k,:) = brgreen1(k,:)/(lnrectc(k)*snrectc(k));
   bzgreen1(k,:) = bzgreen1(k,:)/(lnrectc(k)*snrectc(k));
end

brgreens = zeros(nprobes,nconductors);
bzgreens = brgreens;
bpgreens = brgreens;
for k=1:nprobes
   p_index = find(geom_probe(3,:)==k);
   lenp = length(p_index);
   if lenp>1
      brgreens(k,:) = sum(brgreen1(:,p_index)')/lenp;
      bzgreens(k,:) = sum(bzgreen1(:,p_index)')/lenp;
   else
      brgreens(k,:) = brgreen1(:,p_index)';
      bzgreens(k,:) = bzgreen1(:,p_index)';
%jal 10-3-03 removed      bpgreens(k,:) = bpgreen1(:,p_index)';
   end
end

% Project B-field at probe onto direction of probe measurement.

if nargout>2
   alpha = bprobe_data(3,:)' * pi/180;	% bprobe angles in radians
   for k=1:nconductors
      bpgreens(:,k) = [brgreens(:,k).*cos(alpha) + bzgreens(:,k).*sin(alpha)];
   end
end
