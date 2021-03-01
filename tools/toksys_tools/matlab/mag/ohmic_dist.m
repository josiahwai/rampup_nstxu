function [i_ohmic,fl,br,bz,rbs,zbs,err] = ...
	 ohmic_dist(rbs,zbs,rcc,zcc,drcc,dzcc,mk_bs_uniform,minimize_i)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  USAGE:  
%        [i_ohmic,fl,br,bz,rbs,zbs,err] = ...
%	 ohmic_dist(rbs,zbs,rcc,zcc,drcc,dzcc,mk_bs_uniform,minimize_i)
%
%   Also can be called as script:
%         ohmic_dist
%
%  PURPOSE: Calculates the Best Ohmic current distribution for a set of F coils
%           to a plasma boundary to generate 1Vs in boundary
%
%  INPUTS: [default]
%           
%     rbs, zbs=          R,Z Plasma boundary Points   (m)
%     rcc,zcc,drcc,dzcc= Rectangular Coil center and widths (m)  
%	 (e.g. use the 2nd, 1st, 4th, and 3rd rows of standard fcdata object)
%     mk_bs_uniform=     Makes rbs, zbs a uniform poloidal grid scale [1]
%                        if >1 then makes mk_bs_uniform number of points
%                        if =1 then makes length(rbs) number of uniform points 
%     minimize_i =       =0.1; % Includes currents in opts. to reduce +-I [0]
%                        The value of minimize_i weights how much I smooting 
%                        Value of 0.1 reasonably balances I. err~ 0.25% 
%  OUTPUTS: 
%     i_ohmic =		 Optimum Ohmic Current Dist. for 1Vs on Boundary (A-Turns)
%     fl      =          Flux on Boundary (Vs)
%     br      =          Br-field on Boundary (T)
%     bz      =          Bz-field on Boundary (T)
%     rbs zbs =          R,Z boundary points used in analysis (m)
%                        NOTE: rbs,zbs changed to uniform version of input
%     err     =          Average Error in booundary flux value (Delta_Vs/Vs)
%
%  RESTRICTIONS:
%     uses square coil for all d3d type trapazoidal coils (small effect)
%     Produces Amp-Turns. Use c_ohmic./fcnturn to convert to coil terminal Amp
%     
%
%  METHOD:  
%     Determine optimum Coil current distribution using least square technique to 
%     produce a uniform 1 Vs flux on the boundary

%  WRITTEN BY:  Jim Leuer  ON	7-17-03
%
%  MODIFICATION HISTORY:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Below allow call of routine as script or function
  if nargin <= 0 % READ INPUT FROM CALLING or BASE area
    vnam= 'rbs'; in_script; if ~dumm return; end
    vnam= 'rbs'; in_script; if ~dumm return; end
    vnam= 'zbs'; in_script; if ~dumm return; end
    vnam= 'rcc'; in_script; if ~dumm return; end
    vnam= 'zcc'; in_script; if ~dumm return; end
    vnam= 'drcc'; in_script; if ~dumm return; end
    vnam= 'dzcc'; in_script; if ~dumm return; end
    vnam= 'mk_bs_uniform'; in_script
    vnam= 'minimize_i'; in_script
   end

%  if(nargin<6)
%     disp('ERROR ohmic_dist: Requires minimum of 6 inputs')
%     return;
%  end

% Make input vectors into column vectors, if needed.
  if(size(rcc,1)==1)
     rcc1=rcc';
  else
     rcc1=rcc;
  end
  if(size(zcc,1)==1)
     zcc1=zcc';
  else
     zcc1=zcc;
  end
  if(size(drcc,1)==1)
     drcc1=drcc';
  else
     drcc1=drcc;
  end
  if(size(dzcc,1)==1)
     dzcc1=dzcc';
  else
     dzcc1=dzcc;
  end

  if exist('mk_bs_uniform')~=1, mk_bs_uniform= 1; end

  if exist('minimize_i')~=1,  minimize_i=0; end
  
% efit plasma boundary is very choppy this makes uniform
  if mk_bs_uniform
    dleng= sqrt(diff(rbs).^2 + diff(zbs).^2);
    leng= [0;cumsum(dleng)];
    if mk_bs_uniform >= 2
       lbs= linspace(0,leng(end),mk_bs_uniform)';
    else
       lbs= linspace(0,leng(end),length(rbs))';
    end
    rbs= interp1(leng,rbs,lbs);
    zbs= interp1(leng,zbs,lbs);
%    plot(rbs,zbs,'xr');        plot profile
  end    
  
  
  if rbs(1)==rbs(end) & zbs(1)==zbs(end)  % remove redundant end point
     lbs= lbs(1:end-1);
     rbs= rbs(1:end-1);
     zbs= zbs(1:end-1);
  end
  
  rbso= mean(rbs); % plasma geometric center
  zbso= mean(zbs);
  npt= length(rbs);
  dum= ones(npt,1);
  for ii=1:length(rcc1)
    ri= (rcc1(ii)-0.5*drcc1(ii))*ones(npt,1);
    ro= (rcc1(ii)+0.5*drcc1(ii))*ones(npt,1);
    zl= (zcc1(ii)-0.5*dzcc1(ii))*ones(npt,1);
    zu= (zcc1(ii)+0.5*dzcc1(ii))*ones(npt,1);
    [hr,hz,hf] =fine1(ri,ro,zl,zu,dum,rbs,zbs,0);
    hff(:,ii)=  hf;
    hrr(:,ii)=  hr;
    hzz(:,ii)=  hz;
  end

% construct minimization: (i_ohmic= I*AMU0 to put 1Vs on boundary)
% add in current constraint

   amu0= 4*pi*1.0e-07;
 
  if minimize_i
   amat= [hff; minimize_i*eye(size(rcc1,1))];
   bmat= [ones(size(hff,1),1);zeros(size(rcc1,1),1)];
  else
   amat= hff;
   bmat= ones(size(hff,1),1);
  end
   
  i_ohmic= inv(amat'*amat) * amat' * bmat / amu0; % Opt. Current (A)

  fl= amu0*hff*i_ohmic;       % Flux on boundary (Vs)
  br= amu0*hrr*i_ohmic;       % Br on boundary (T)
  bz= amu0*hzz*i_ohmic;       % Bz on boundary (T)
  bbgauss= sqrt(br.^2+bz.^2)'*1e+4; % - Field on boundary

  err= sum(abs(fl-1))/length(fl);
  if err > .005
    disp(['% CAUTION: Average error in bdry. flux is greater than 0.5% err= ',...
          num2str(err*100), '%']);
    if minimize_i
       disp(['% CAUTION: You should reduce minimize_i ', num2str(minimize_i)]);
    end
  end

  if nargout<=0 % output if called as script
    vnam= 'i_ohmic', ot_script;
    vnam= 'fl', ot_script;
    vnam= 'br', ot_script;
    vnam= 'bz', ot_script;
    vnam= 'rbs', ot_script;
    vnam= 'zbs', ot_script;
    vnam= 'err', ot_script;
  end
 
  return

% ======================================================================
% testing
% ======================================================================

  load_d3denv
  filename= '/u/leuer/efit/d3d/s118526/g118526.04000';
  read_gfile
  id= 1:18; % use only PF coils,
  rcc=  fcdata(2,id)';
  zcc=  fcdata(1,id)';
  drcc= fcdata(4,id)';
  dzcc= fcdata(3,id)';
  rbs= rbbbs;
  zbs= zbbbs;

  idd= [1,2,3,4,5, 8, 9, 7, 6]; % d3d clockwise coil id from center
  idd= [idd idd(9:-1:1)+9]

% No Current Optimimization
  minimize_i= 0;

% work from base area for all type of calls
  ohmic_dist
 [i_ohmic,fl,br,bz,rbs,zbs,err] = ...
	 ohmic_dist(rbs,zbs,rcc,zcc,drcc,dzcc);
 [i_ohmic,fl,br,bz,rbs,zbs,err] = ...
	 ohmic_dist 

% seems to work from base area for all type of calls
	 
  figure(1), clf
  
  plot(1:18,i_ohmic(idd))
  hold on
  grid on
  title('OHNMIC\_DIST.m: Blue=No Current Opt, Red= Current Opt (0.1)')
  xlabel('DIII-D F-coils, Clockwise from F1A,(Shot 118526) ')
  ylabel('DIII-D Coil Currents')  

  figure(2), clf
  plot(lbs,fl)
  hold on
  grid on
  title('OHNMIC\_DIST.m: Blue=No Current Opt, Red= Current Opt (0.1)')
  xlabel('Distance around Plasma Boundary, (Shot 118526)')
  ylabel('Flux on Plasma Boundary')

% With some Current Optimimization
  minimize_i= 0.1;
  ohmic_dist  
  
  figure(1)
  plot(1:18,i_ohmic(idd),'r')
  
  figure(2)
  plot(lbs,fl,'r')


  
