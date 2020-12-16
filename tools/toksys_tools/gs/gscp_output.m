%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gscp_output
%
%  PURPOSE: Populate output vector with contents controlled by index_in_y
%
%  INPUTS: index_in_y and workspace after gscp_analysis_response
%
%  OUTPUTS:  lae.y
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander ON 2015-03-04
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(index_in_y,'psizr')
  lae.y(index_in_y.psizr) = psizr;
end
if isfield(index_in_y,'ic')
  lae.y(index_in_y.ic) = ic;
end
if isfield(index_in_y,'iv')
  lae.y(index_in_y.iv) = iv;
end
if isfield(index_in_y,'sp')
  lae.y(index_in_y.sp) = sp;
end
if isfield(index_in_y,'sf')
  lae.y(index_in_y.sf) = sf;
end
if isfield(index_in_y,'er')
  lae.y(index_in_y.er) = er;
end
if isfield(index_in_y,'fl')
  fl = mlc*ic+mlv*iv+mpl'*pcurrt(:);
  lae.y(index_in_y.fl) = fl;
end
if isfield(index_in_y,'lv')
  lv = mhc*ic+mhv*iv+mph'*pcurrt(:);
  lae.y(index_in_y.lv) = lv;
end
if isfield(index_in_y,'bp')
  bp = gbc*ic+gbv*iv+gpb'*pcurrt(:);
  lae.y(index_in_y.bp) = bp;
end
if isfield(index_in_y,'rog')
  rog = rldata(:,indic)*ic+rldata(:,indiv)*iv+rldata(:,nc+nv+1)*cpasma;
  lae.y(index_in_y.rog) = rog;
end
if isfield(index_in_y,'psidp')
  psidp = mdc*ic+mdv*iv+mpd'*pcurrt(:);
  lae.y(index_in_y.psidp) = psidp;
end
if isfield(index_in_y,'brdp')
  brdp = grdc*ic+grdv*iv+grpd'*pcurrt(:);
  lae.y(index_in_y.brdp) = brdp;
end
if isfield(index_in_y,'bzdp')
  bzdp = gzdc*ic+gzdv*iv+gzpd'*pcurrt(:);
  lae.y(index_in_y.bzdp) = bzdp;
end
if isfield(index_in_y,'btdp')
  btdp = gs_interp2(rg,zg,fpolzr,rdp,zdp)./rdp;
  lae.y(index_in_y.btdp) = btdp;
end
if isfield(index_in_y,'gapdp')
  gscp_gapdp
  lae.y(index_in_y.gapdp) = gapdp;
end
if isfield(index_in_y,'gaps')
  gscp_gaps
  lae.y(index_in_y.gaps) =  gaps;
end
if isfield(index_in_y,'ys')
  lae.y(index_in_y.ys) = ys;
end
kxpoints = 1:nxpoints; % Vector for sorting x-points that go into y
if xpointorder == 2
  [~, kxpoints(1:n)] = sort(nulls.z(1:n));
end
if isfield(index_in_y,'rx')
  rxoutput = zeros(nxpoints,1); % nan's can't be used in Simulink so using zeros
  n = min(nxpoints,nulls.count);
  rxoutput(1:n) = nulls.r(kxpoints(1:n));
  lae.y(index_in_y.rx) = rxoutput;
end
if isfield(index_in_y,'zx')
  zxoutput = zeros(nxpoints,1); % nan's can't be used in Simulink so using zeros
  n = min(nxpoints,nulls.count);
  zxoutput(1:n) = nulls.z(kxpoints(1:n));
  lae.y(index_in_y.zx) = zxoutput;
end
if isfield(index_in_y,'psix')
  psixoutput = zeros(nxpoints,1); % nan's can't be used in Simulink so using zeros
  n = min(nxpoints,nulls.count);
  psixoutput(1:n) = nulls.psi(kxpoints(1:n));
  lae.y(index_in_y.psix) = psixoutput;
end
if isfield(index_in_y,'rcur')
  lae.y(index_in_y.rcur) = rcur;
end
if isfield(index_in_y,'zcur')
  lae.y(index_in_y.zcur) = zcur;
end
if isfield(index_in_y,'cpasma')
  lae.y(index_in_y.cpasma) = cpasma;
end
if isfield(index_in_y,'aminor')
  lae.y(index_in_y.aminor) = a0;
end
if isfield(index_in_y,'nbbbs')
  lae.y(index_in_y.nbbbs) = nbbbs_max;
end
if isfield(index_in_y,'rbbbs')
  lae.y(index_in_y.rbbbs) = rbbbs;
end
if isfield(index_in_y,'zbbbs')
  lae.y(index_in_y.zbbbs) = zbbbs;
end
if isfield(index_in_y,'rmaxis')
  lae.y(index_in_y.rmaxis) = rmaxis;
end
if isfield(index_in_y,'zmaxis')
  lae.y(index_in_y.zmaxis) = zmaxis;
end
if isfield(index_in_y,'psimag')
  lae.y(index_in_y.psimag) = psimag;
end
if isfield(index_in_y,'rbdef')
  lae.y(index_in_y.rbdef) = rbdef;
end
if isfield(index_in_y,'zbdef')
  lae.y(index_in_y.zbdef) = zbdef;
end
if isfield(index_in_y,'psibry')
  lae.y(index_in_y.psibry) = psibry;
end
if isfield(index_in_y,'li')
  lae.y(index_in_y.li) = li;
end
if isfield(index_in_y,'betap')
  lae.y(index_in_y.betap) = betap;
end
if isfield(index_in_y,'Wth')
  lae.y(index_in_y.Wth) = Wth;
end
if isfield(index_in_y,'psipla')
  lae.y(index_in_y.psipla) = psipla;
end
if isfield(index_in_y,'bp2flx')
  lae.y(index_in_y.bp2flx) = bp2flx;
end
if isfield(index_in_y,'psiplaapp')
  lae.y(index_in_y.psiplaapp) = psiplaapp;
end
if isfield(index_in_y,'lconn') % Connection length calculation goes here
  lae.y(index_in_y.lconn) = 1000;
end
if isfield(index_in_y,'Ltot')
  lae.y(index_in_y.Ltot) = Ltot;
end
if isfield(index_in_y,'Atot')
  lae.y(index_in_y.Atot) = Atot;
end
if isfield(index_in_y,'Vtot')
  lae.y(index_in_y.Vtot) = Vtot;
end
if isfield(index_in_y,'rhot')
  lae.y(index_in_y.rhot) = rhot;
end
if isfield(index_in_y,'jtav')
  lae.y(index_in_y.jtav) = jtav;
end
if isfield(index_in_y,'Vres')
  lae.y(index_in_y.Vres) = Vres;
end
for i = 1:length(user_signal)
  if isfield(index_in_y,user_signal(i).name)
    lae.y(getfield(index_in_y,user_signal(i).name)) = ...
      user_signal(i).y0 + user_signal(i).dydic*ic + ...
      user_signal(i).dydiv*iv + user_signal(i).dydpcurrt*pcurrt(:);
  end
end
if isfield(index_in_y,'fluxerror')
  lae.y(index_in_y.fluxerror) = max(abs(psiherr)/(psibry-psimag));
end
if isfield(index_in_y,'time')
  lae.y(index_in_y.time) = time;
end
