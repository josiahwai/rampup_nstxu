function [brgreens,bzgreens,bpgreens,lnrectc,snrectc] = ...
		bgreens_fine(coil_data,bprobe_data,nrect,verbose)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)bgreens_fine.m	1.4 02/11/08
%
%  SYNTAX: [brgreens,bzgreens,bpgreens] = ...
%		bgreens_fine(coil_data,bprobe_data,nrect,verbose)
%
%  PURPOSE:  Compute Greens functions from axisymmetric current carrying 
%	conductors, each of whose cross-section is a parallelogram, to Br,
%	Bz, and angled probe measurement at a set of specified points.  
%	This version of code is based on low level calculation using fine.m.
%	Can be used to compute Br,Bz Greens functions for arbitrary points in
%	space.
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
%		bprobe_data(1:6,:)=[Z;R;Ang;length(m);Ip coeff;width(m)]
%     (abs(4th row) specifies width if < 0)
%		(5th row is ignored, so 4 x number probes also works)
%     (6th row is optional, assumed 0.02m if not provided)
%	nrect = rectangle partitioning for conductor set (see header
%		documentation for bld_subelements.m)
%       verbose = print/display verbose info.  0 = no print/display.  Numbers>0
%               give increasingly more information/displays.
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
%     If probe width is specified in 6 row bprobe_data AND there is length<0 
%     in row 4, the row 4 specified length prevails. Probe will be treated as
%     flat if width specified in row 4 (6 row convention allows l & w)
%
%  METHOD:  Each conductor is partitioned into a bunch of tiny rectangles.  All
%	bgreens are computed for these rectangles, then combined at the end.
%	Calculation is done in 2 steps:
%	(1) use filament to filament calculation for all subelements
%	(2) replace those subelement calculations which would have inaccurate
%	results, based on geometric information, with more accurate calculation.
% 
%  WRITTEN BY:  Mike Walker 	ON 	7/26/99
%  10/03/03   JAL1   remove one line commented below
%  10/30/07   NWE    length=0 "saddle probes" treated as single point
%  02/08/08   NWE    applies probe cross section width if specified in
%                    bprobe_data(6,:) or length<0 in bprobe_data(4,:).  
%                    
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
bpdata(5,:) = 0;   % ignore Ip coeff
probe_length = bprobe_data(4,:);
probe_angle = bprobe_data(3,:);
if(size(bprobe_data,1) == 6) % x-sections specified
   xsec = bprobe_data(6,:);
else
   xsec = 0.02 * ones(1,nprobes); % default width = 2 cm, see JAL note below
end
% deal with saddle loops
for k=1:nprobes
   if probe_length(k) < 0 % EFIT style saddle loop: -length = x-section
      tmp = xsec(k);
      xsec(k) = abs(probe_length(k)); % move to standard form
      probe_length(k) = tmp;
   end
   if xsec(k) > abs(probe_length(k)) % saddle loop
      tmp = xsec(k);
      xsec(k) = probe_length(k);
      probe_length(k) = tmp;
      probe_angle(k) = probe_angle(k) + 90; % rotate to allow proper decomposition
   end
end

cosangles = cos(probe_angle*pi/180);
sidx = find(abs(cosangles)<0.5);
all_index = 1:nprobes;
cidx = setdiff(all_index,sidx);   
% Note conversion of probe to trapezoid format increases dimension by one
% this is reason for dZ,dR=2cm for second dimension below (jal)
for k = cidx
   bpdata(3,k) = xsec(k);
   bpdata(4,k) = abs(probe_length(k).*cos(probe_angle(k)*pi/180));		
   bpdata(5,k) = probe_angle(k); % make all AC=angles
   bpdata(6,k) = 0;
end
for k = sidx
   bpdata(4,k) = xsec(k);
   bpdata(3,k) = abs(probe_length(k).*sin(probe_angle(k)*pi/180));			
   bpdata(5,k) = 0;
   bpdata(6,k) = probe_angle(k);   % make all AC2=angles
end

if (size(bprobe_data,1) < 6)  % backwards compatibility - treat as lines
   [geom_probe,lnrectp,snrectp] = bld_subelements(bpdata,1,1);
else                         % treat as rectangles
   [geom_probe,lnrectp,snrectp] = bld_subelements(bpdata,max(max(nrect)),1);
end

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

      mu0= 4*pi*1e-7;
      ri = ra-dra/2;
      ro = ra+dra/2;
      zu = za+dza/2;
      zl = za-dza/2;
      lb = length(rb);
%      [brg,bzg,psi] =fine1(ri*ones(lb,1),ro*ones(lb,1),zl*ones(lb,1), ...
%		zu*ones(lb,1),ones(lb,1),rb,zb);
      [brg,bzg,psi] =fine1(ri*ones(1,lb),ro*ones(1,lb),zl*ones(1,lb), ...
		zu*ones(1,lb),ones(1,lb),rb,zb);
      brg = mu0*brg;
      bzg = mu0*bzg;
      psi = mu0*psi;

      if 0
         [psi0,brg0,bzg0] = rect2mag(ra,za,dra,dza,rb,zb);
         figure(10),clf
         subplot(3,1,1)
         plot((brg0-brg)./brg0)
         ylabel('Br error')
         subplot(3,1,2)
         plot((bzg0-bzg)./bzg0)
         ylabel('Bz error')
         subplot(3,1,3)
         plot((psi0-psi)./psi0)
         ylabel('Psi error')
         wait
      end

      brg = brg*1e6;		% Tesla/MA
      bzg = bzg*1e6;		% Tesla/MA

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
