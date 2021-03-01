function [Iq, t, Iqn] = getIq(shot,q,t1,t2,efit_source,tokamak)
%
%  USAGE: [Iq, t, Iqn] = getIq(shot,q,t1,t2,efit_source,tokamak)
%
%  PURPOSE: get current within q-surface
%
%  INPUTS: shot: shot number
%          q: q-values for surfaces (default 2)
%          t1, t2: Specification of the time samples (t)
%            t1 <= t <= t2, if length(t1)==1 & length(t2)==1
%            t = t1,        if length(t1)>=1 & length(t2)==0
%          efit_source: (default efit02)
%          tokamak: (default d3d)
%
%  OUTPUTS: Iq: current within q-contours
%            t: times [sec]
%          Iqn: Iq/cpasma, i.e. fraction of current within q-surfaces

%
%  RESTRICTIONS: 

%
%  METHOD: 
%	
%  VERSION @(#)getIq.m	1.3 08/30/11
%
%  WRITTEN BY:  Anders Welander  ON	5/17/11
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  twopi = 2*pi;
  mu0 = .4e-6*pi;
  
  % Defaults
  if ~exist('q','var')
    q = 2;
  end
  if ~exist('t1','var')
    t1 = [];
  end
  if ~exist('t2','var')
    t2 = [];
  end
  if ~exist('efit_source','var')
    efit_source = 'efit02';
  end
  if ~exist('tokamak','var')
    tokamak = 'd3d';
  end
  
  % Read the efits
  equilibria = read_eq(shot,[],efit_source,tokamak);
  
  % Select time samples
  if     length(t1) == 0 & length(t2) == 0
    iselect = 1:length(equilibria.time);
  elseif length(t1) == 1 & length(t2) == 0
    [dum, j] = min(abs(equilibria.time-t1)); iselect = j;
  elseif length(t1) == 1 & length(t2) == 1
    iselect = find(t1<=equilibria.time & equilibria.time<=t2);
  elseif length(t1) > 1
    iselect = find(ismember(equilibria.time,t1));
  end
  t = equilibria.time(iselect);
  
  % Make sure the outputs are defined even if length(t) or nq = 0
  nq = length(q);
  Iq = zeros(nq,length(t));
  Iqn = zeros(nq,length(t));  
  
  % Process the equilibria
  for jeq = 1:length(iselect)
    % Get this eq
    if length(equilibria.time) == 1
      eq = equilibria;
      if isfield(eq,'gdata')
        eq = eq.gdata;
      end
    else
      eq = equilibria.gdata(iselect(jeq));
    end
    struct_to_ws(eq);
    % Find the flux for the specified q values, psiq
    psigrid = linspace(psimag,psibry,nw)'; % Wb
    for j = 1:nq
      if q(j)<max(qpsi) & q(j)>min(qpsi)
	k = nw-1; % Searching for q(j) value in qpsi array
	while (qpsi(k)-q(j))*(qpsi(k+1)-q(j)) > 0 & k > 1
	  k=k-1;
	end
	k2 = min(nw-2,k+1); while (qpsi(k2+1)-qpsi(k2))*(qpsi(k2+2)-qpsi(k2+1))>0 & k2<nw-2, k2 = k2+1; end
	k1 = min(nw-2,k); while (qpsi(k1+2)-qpsi(k1+1))*(qpsi(k1+1)-qpsi(k1))>0 & k1>1, k1 = k1-1; end
	psiq(j) = spline(qpsi(k1:k2),psigrid(k1:k2),q(j));
      else
	psiq(j) = NaN;
      end
    end
    % Divide plasma into pies
    abbbs = angle(rbbbs-rmaxis+i*(zbbbs-zmaxis));
    [a,ib] = sort(abbbs(1:nbbbs-1)); % To avoid phase jumps
    rc = zeros(nw,nbbbs-1);
    zc = zeros(nw,nbbbs-1);
    for j = 1:nbbbs-1
      rc(:,j) = linspace(rmaxis,rbbbs(ib(j)),nw)';
      zc(:,j) = linspace(zmaxis,zbbbs(ib(j)),nw)';
    end
    psic = interp2(rg,zg,psizr,rc,zc,'cubic');
    rhoc = sqrt((rc-rmaxis).^2+(zc-zmaxis).^2);
    pprimec = spline(psigrid,pprime,psic);
    ffprimc = spline(psigrid,ffprim,psic);
    jphic = rc.*pprimec+ffprimc/mu0./rc;
    da = [a(1)-a(end)+twopi; diff(a)];
    da = (da+da([2:end 1]))/2; % sum(da)=twopi, each pie-cut has half the angles to its neighbors
    Ip(jeq) = 0;
    Iq(1:nq,jeq) = 0;
    for j = 1:nbbbs-1
      ii = 1; % Will hold indices of monotonically decreasing flux
      for k=2:nw
        if psic(k,j)<psic(ii(end),j), ii(end+1)=k; end
      end
      Asector = da(j)*rhoc(:,j).^2/2;
      Isector(2:nw,1) = cumsum((jphic(1:end-1,j)+jphic(2:end,j)).*diff(Asector))/2;
      Ip(jeq) = Ip(jeq)+Isector(end);
      Iq(:,jeq) = Iq(:,jeq) + spline(psic(ii,j),Isector(ii),psiq(:));
    end
    Iqn(:,jeq) = Iq(:,jeq)/Ip(jeq);
    Iq(:,jeq) = Iqn(:,jeq)*cpasma;
    cpasmas(jeq) = cpasma;
  end
  
  
