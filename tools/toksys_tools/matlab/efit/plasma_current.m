  function [pcur,jphi,cpasma]= plasma_current(psizr, pprime, ffprim,...
                                              psimag,psibnd,...
                                              rgefit,zgefit,...
					      rbbbs,zbbbs,nbbbs,...
					      iplot,dofast,cpasma)
% Function generates plasma current (pcurrt) from flux (psizr)
% 
% SYNTAX:
%  [pcur,jphi,cpasma]= plasma_current(psizr, pprime, ffprim,...
%                                     psimag,psibnd,...
%                                     rgefit,zgefit,...
%				      rbbbs,zbbbs,nbbbs,...
%				      iplot,dofast)
%
%  PURPOSE: Generate plasma current from flux and current profile information
%           Input is same as read_gfile input or mds_eq => eq_ga_env
%           pcur is almost identical to that obtained using iecurr=2 in EFIT
%           (error <~ 0.2% and max near boundary)
%
%  INPUTS: [default] <optional>
%
%          psizr=   true flux on grid (Wb)
%          pprime=  p' vector  (pprime from read_gfile) (nt/m^2)/(Vs/rad)
%          ffprim=  ff' vector (ffprim in read_gfile)   (mT)^2/(Vs/rad)
%          psimag=  Flux on axis (Wb) [max(max(psizr))]
%          psibnd=  Flux on boundary (Wb)  [average of flux on boundary]
%          rgefit=  Radius vector of grid (m)
%          zgefit=  Z vector of grid (m)
%          rbbbs=   Radius vector of boundary (m)
%          zbbbs=   Z vector of boundary (m)
%          nbbbs=   number of elements in rbbbs zbbbs [length(rbbbs]
%          iplot=   2; % 1= plot boundary, 2= contour jphi; 3= clabel [0]
%          dofast=  1; % increases algorithm speed at risk of including some
%                   private flux in solution. Only extreme shapes like
%                   crescents and beans are expected to need dofast=0.
%                   Default is dofast=0 since speed is still very fast. [0]
% 
%         <cpasma>= Plasma Current (A)
%
%  OUTPUTS:
%          pcur   = plasma current on nr,nz grid [A]
%          jphi=    plasma current density on nr,nz grid [MA/m^2]
%          cpasma=  plasma current = sum(sum(pcur))
%
%          Warning messages displayed if sum(sum(pcur)) and existing cpasma
%          are widely different. cpasma then overwritten.
%
%  NOTE:    ALL UNITS ARE MKS (except jphi is MA/m^2)
%  CAUTION: Overwrites pcur, jphi & cpasma if read_gfile read from gfile
%
%  RESTRICTIONS:
%
%  SEE ALSO:  inside_plasma

%  WRITTEN BY:  Jim Leuer ON 6-3-03	
%
%  Routines Used: inside_plasma
%
%  MODIFICATION HISTORY: 
%  jal30sep10: +cpasma normalization, spline fit; garea per efit
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)plasma_current.m	1.4 07/15/09

% Default:
% Note: cpasma is 12 argument

  if nargin==11
    dofast= 0;
  elseif nargin==10
    iplot=0;
    dofast= 0;
  elseif nargin==9
    nbbb= length(rbbb);
  elseif nargin<=8
    disp('%ERROR plasma_current: Needs at least 9 arguments to work')
    return
  end
  
  if isempty(psimag), 
     disp('% CAUTION: plasma_current: psimag is empty using max on grid')
     psimag= max(max(psizr)),
  end

  if isempty(psibnd) % average on boundary
    disp('% CAUTION: plasma_current: psibnd is empty using ave. on boundary')
    temp = interp2(rgefit,zgefit,psizr,rbbbs,zbbbs);
    idx = find(isnan(temp));
    idx = setdiff(1:length(temp),idx);
    temp = temp(idx);
    psibnd= mean(temp);
  end
  if exist('iplot')~=1, iplot=0; end   % no plotting
  if exist('dofast')~=1, dofast=0; end % 100% guaranteed mapping

% ====================================================
% Find Plasma Grid that is "INSIDE" plasma Boundary:

  inside= inside_plasma(psizr,rgefit,zgefit,rbbbs,zbbbs,nbbbs, psimag,psibnd,...
                        iplot,dofast);  

% ==================================================================
% Compute (see GA-A14490, pg8): Jphi= r*Pprime + F*Fprime/(amu0*r)

  garea=   (rgefit(2)-rgefit(1))*(zgefit(2)-zgefit(1)); % EFIT area def
  areaomu= garea/(4*pi*1e-7);
  [nz,nr]= size(psizr);
  jphi=    zeros(nz,nr);             % initialize to zero current
  pcur=    jphi;                     % Plasma current on grid (initial 0)
  rg=      ones(nz,1)*rgefit';       % R on grid (same size as psizr)
  id=      find(inside);                   % work only on inside plasma points
  if(isempty(id))	% handle badly converged EFIT
     cpasma = 0;
     return;
  end
  rr=      rg(id);                         % reduced set of R's Inside plasma
  psi=     linspace(psimag,psibnd,length(pprime))'; % from center to edge
  pp=      interp1(psi,pprime,psizr(id),'spline');  % Pprime
  ffp=     interp1(psi,ffprim,psizr(id),'spline');  % F*Fprime
  pcur(id)=  garea*rr.*pp + areaomu*ffp./rr;  % Toroidal current density (inside)
  cpas= sum(pcur(id));
  if exist('cpasma')==1 & ~isempty(cpasma)
    if abs((cpas-cpasma)/cpasma) > 0.01
     disp(['%CAUTION: plasma_current sum(pcur)-cpasma>1.0% want,got: ',...
          num2str(cpasma),' ',num2str(cpas)])
     disp(['%         Re-Normalizing current distribution to cpasma: ', ...
          num2str(cpasma)])
    end
     pcur(id)= pcur(id)*cpasma/cpas;
     cpas= cpasma;
  end

  jphi(id)=  pcur(id)/(garea*1e+6);           % Convert to MA/m^2 (daves Convention)
  cpasma= cpas; % set plasma current to sum of grid current

  if iplot>=2
    hold on
    [c,h]= contour(rgefit,zgefit,jphi);
    if iplot>=3
      clabel(c,h)
    end
    hold on
    plot(rbbbs(1:nbbbs),zbbbs(1:nbbbs),'k');
    hold on
    axis image
    xlabel('R [m]')
    ylabel('Z [m]')
    title([' Plasma Current Density [MA/m^2], Ip= ',...
           num2str(1e-6*cpasma), ' MA']);  
  end
  
  return

% =============================================================================
% Testing (FOR MAIN TESTING SEE Leuer's efit area: test_plasma_current)

%  filename= '/u/leuer/efit/diiid/g110214.01740'; % DIIID file with: iecurr=2
  filename= '/u/leuer/efit/diiid/shot113697/g113697.03000'; %
  read_gfile
  cpasma0= cpasma;
  if exist('pcurrt')
    pcurrt0= pcurrt; % save gfile currents for later error analysis
  end

  iplot=3;
  dofast=0;
  
% do plasma_current

  figure(3)
  clf
  iplot=3;
  [pcur,jphi,cpasma]= plasma_current(psizr, pprime, ffprim,...
                                     psimag,psibnd,...
                                     rgefit,zgefit,...
				     rbbbs,zbbbs,nbbbs,...
				     iplot,dofast);
% do some checks to see error
 if exist('pcurrt0')==1
  figure(7)
  clf
  colormap('default');
  cmap= colormap;
  cmap= [cmap(1:end-1,:);[1 1 1]];
  colormap(cmap);
  err1= zeros(size(pcurrt0));
  id= find(pcurrt0~=0);
  err1(id)= (pcur(id)-pcurrt0(id))./pcurrt0(id);  
  contourf(rgefit,zgefit,100*err1)
  hold on
  plot(rbbbs(1:nbbbs),zbbbs(1:nbbbs),'k');
  hold on
  axis image  
  title(['%Error= 100*(pcur-pcurrt0)/pcurrt0, Max value= ',...
         num2str(100*max(max(abs(err1)))),' %'])
  colorbar('vertical')
  
% Results: Max error for 98549.04000: 0.11% toward edge Most error <.04%
 end

