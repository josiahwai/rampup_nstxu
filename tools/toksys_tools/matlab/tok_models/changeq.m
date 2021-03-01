function [eq, eqx] = changeq(eq0,eq_target,tok_data_struct,options,idoplot)
%
%  USAGE:   [eq, eqx] = changeq(eq0,eq_target,tok_data_struct,options,idoplot)
%
%  PURPOSE: Change an equilibrium from eq0 to eq_target
%           Targeted variables are: (rbbbs,zbbbs), cc, cpasma, betap, li
%           Bus constraint such as the VFI on DIII-D can be imposed on cc
%
%  INPUTS:  eq0: Original equilibrium structure on toksys format
%           eq_target: Equilibrium structure containing target quantities 
%             Specify NaN for variables to adjust in order to achieve target
%             example: if cc can be adjusted then set eq_target.cc(:) = NaN;
%           tok_data_struct: standard toksys object for tokamak description
%           options.converrmax = upper limit on flux error in resulting eq
%           options.maxiter: maximum convergence iterations, default = 25;
%           options.cccirc: standard toksys assignment of circuit numbers to coils,
%           options.bus_code: indices are 1 for coils on bus, 0 for others
%              For DIII-D do: PP_objs =  get_PP_objs(shot)
%                             options.bus_code = [0 0 PP_objs.bus_code]
%           idoplot: plot each iteration for idoplot seconds, default 0 = no plot
%
%  OUTPUTS: eq: structure containing new converged equilibrium on tok_data_struct grid
%           eqx: structure containing extra information (see eqx.descriptions)
%	
%  RESTRICTIONS: MAXIMUM grid size is 65x65
%                Only these quantities are updated for the new equilibrium:
%    		   rbbbs, zbbbs, jphi, psizr, pprime, ffprim, pres, fpol
%                  cc, psimag, psibry, cpasma, rmaxis, zmaxis, rg, zg, dr, dz
%                  nw, nh, psirz, ssibry, ssimag, qpsi
%
%  METHOD:  Linear solution is found by gspert and converged with convergeq
	
%  VERSION @(#)changeq.m	1.4 01/25/12
%
%  WRITTEN BY:  Anders Welander  ON	5/6/11
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  mu0 = 0.4e-6*pi;

  options.iconstraints = 1;

  if ~exist('idoplot','var')
    idoplot = 0;
  end
  
  struct_to_ws(tok_data_struct);
  dr = rg(2)-rg(1);
  dz = zg(2)-zg(1);
  
  % If we want to control the boundary then specify it as iso flux points
  if isempty(find(isnan(eq_target.rbbbs+eq_target.zbbbs)))
    options.iso = [eq_target.rbbbs(1:eq_target.nbbbs) eq_target.zbbbs(1:eq_target.nbbbs)];
    control_boundary = 1;
    eq0.nbbbs = eq_target.nbbbs;
    eq0.rbbbs = eq_target.rbbbs;
    eq0.zbbbs = eq_target.zbbbs;
  else
    control_boundary = 0;
  end
  
  [response, eqx] = gspert(eq0,tok_data_struct,options);
  
  % Matrix to convert from cc to ic (ic is aka cc0 in other codes)
  eq = eq0;
  ncc = length(eq.cc);
  for j = 1:ncc
    eq.cc = eq.cc*0;
    eq.cc(j) = 1;
    equil_I = cc_efit_to_tok(tok_data_struct,eq);
    iccc(:,j) = equil_I.cc0t;
  end
  eq = eq0;
  
  dcc = eq_target.cc - eq0.cc;
  
  % Transform dcc,iccc for circuits if cccirc exists
  if isfield(options,'cccirc')
    cccirc = options.cccirc;
    ncc = max(abs(cccirc));
    iccc2 = zeros(nc,ncc);
    dcc2 = zeros(ncc,1);
    for j = 1:ncc % j = number of circuit
      k = find(abs(cccirc)==j); % k = numbers of the coils in ciruit j      
      n = length(k);
      for l = k
        s = sign(cccirc(l));
        dcc2(j) = dcc2(j) + s*dcc(l)/n;
	iccc2(:,j) = iccc2(:,j) + s*iccc(:,l);
      end
    end
    dcc = dcc2;
    iccc = iccc2;
  end  
  
  icx = find(isnan(dcc)); % Indices to coils that can be used to achieve target 
  ncx = length(icx);
  
  dip = eq_target.cpasma - eq0.cpasma;
  
  if ~isfield(eq0,'betap'), eq0.betap = eqx.betap; end
  if ~isfield(eq_target,'betap'), eq_target.betap = eqx.betap; end
  dbetap = eq_target.betap - eq0.betap;
  
  if ~isfield(eq0,'li'), eq0.li = eqx.li; end
  if ~isfield(eq_target,'li'), eq_target.li = eqx.li; end
  dli = eq_target.li - eq0.li;
  
  % Calculate changes needed to acheive eq_target based on linear responses
  if control_boundary & ncx>0 % if true, control boundary by changing cc(icx)
    if eqx.ilimited
      disoflux = eqx.isopsi-eq0.psibry;
      for j = 1:eq0.nbbbs
        if ~isnan(dip), disoflux(j) = disoflux(j)+(response.disopsidip(j)-response.dpsibrydip)*dip; end
        if ~isnan(dli), disoflux(j) = disoflux(j)+(response.disopsidli(j)-response.dpsibrydli)*dli; end
        if ~isnan(dbetap), disoflux(j) = disoflux(j)+(response.disopsidbetap(j)-response.dpsibrydbetap)*dbetap; end
        disofluxdcc(j,:) = (response.disopsidis(j,1:nc)-response.dpsibrydis(1:nc))*iccc;
      end
      dcc(icx) = -pinv([disofluxdcc(:,icx); eye(ncx)])*[disoflux; zeros(ncx,1)];
    else
      [d2, ix] = min((eq0.rbbbs-eqx.r0).^2+(eq0.zbbbs-eqx.z0).^2);
      disoflux = eqx.isopsi-eqx.isopsi(ix);
      for j = 1:eq0.nbbbs
        if ~isnan(dip), disoflux(j) = disoflux(j)+(response.disopsidip(j)-response.disopsidip(ix))*dip; end
        if ~isnan(dli), disoflux(j) = disoflux(j)+(response.disopsidli(j)-response.disopsidli(ix))*dli; end
        if ~isnan(dbetap), disoflux(j) = disoflux(j)+(response.disopsidbetap(j)-response.disopsidbetap(ix))*dbetap; end
        disofluxdcc(j,:) = (response.disopsidis(j,1:nc)-response.disopsidis(ix,1:nc))*iccc;
      end
      dbzxdcc = response.disobzdis(ix,1:nc)*iccc;
      dbrxdcc = response.disobrdis(ix,1:nc)*iccc;
      bzx = eqx.isobz(ix);
      brx = eqx.isobr(ix);
      if ~isnan(dip), bzx = bzx+response.disobzdip(ix)*dip; brx = brx+response.disobrdip(ix)*dip; end
      if ~isnan(dli), bzx = bzx+response.disobzdli(ix)*dli; brx = brx+response.disobrdli(ix)*dli; end
      if ~isnan(dbetap), bzx = bzx+response.disobzdbetap(ix)*dbetap; brx = brx+response.disobrdbetap(ix)*dbetap; end
      wx = 100; % how hard to control Br, Bz to zero at target x-point
      dcc(icx) = -pinv([disofluxdcc(:,icx); wx*dbzxdcc(icx); wx*dbrxdcc(icx); eye(ncx)])*[disoflux; wx*bzx; wx*brx; zeros(ncx,1)];      
    end
  end
  
  % Change the equilibrium by calculated perturbation
  dic = iccc*dcc;
  dcphi = response.dcphidis(:,1:nc)*dic+response.dcphidip(:)*dip+response.dcphidbetap(:)*dbetap+response.dcphidli(:)*dli;
  dpprime = response.dpprimedis(:,1:nc)*dic+response.dpprimedip(:)*dip+response.dpprimedbetap(:)*dbetap+response.dpprimedli(:)*dli;
  dffprim = response.dffprimdis(:,1:nc)*dic+response.dffprimdip(:)*dip+response.dffprimdbetap(:)*dbetap+response.dffprimdli(:)*dli;
  dpsimag = response.dpsimagdis(:,1:nc)*dic+response.dpsimagdip(:)*dip+response.dpsimagdbetap(:)*dbetap+response.dpsimagdli(:)*dli;
  dpsibry = response.dpsibrydis(:,1:nc)*dic+response.dpsibrydip(:)*dip+response.dpsibrydbetap(:)*dbetap+response.dpsibrydli(:)*dli;
  eq.psizr = eq0.psizr + gscalc(reshape(dcphi/1e6/dr/dz,nz,nr),rg,zg,mpp)+reshape(mpc*dic,nz,nr);
  eq.pprime = eq0.pprime+dpprime;
  eq.ffprim = eq0.ffprim+dffprim;
  eq.psimag = eq0.psimag+dpsimag;
  eq.psibry = eq0.psibry+dpsibry;
  eq.rbbbs = eq_target.rbbbs;
  eq.zbbbs = eq_target.zbbbs;
  psibarzr = (eq.psizr-eq.psimag)/(eq.psibry-eq.psimag);
  
  % Grid points that can be considered in the plasma have a maximum initial rho compared to bbbs
  rhobbbs = sqrt((eq.rbbbs(2:eq.nbbbs)-eq.rmaxis).^2+(eq.zbbbs(2:eq.nbbbs)-eq.zmaxis).^2);
  thbbbs = angle((eq.rbbbs(2:eq.nbbbs)-eq.rmaxis)+sqrt(-1)*(eq.zbbbs(2:eq.nbbbs)-eq.zmaxis));
  rhogg = sqrt((rgg-eq.rmaxis).^2+(zgg-eq.zmaxis).^2);
  thgg = angle((rgg-eq.rmaxis)+sqrt(-1)*(zgg-eq.zmaxis));
  rhomaxgg = interp1([thbbbs-2*pi; thbbbs; thbbbs+2*pi],[rhobbbs; rhobbbs; rhobbbs],thgg);
  
  ivac = find(psibarzr>1 | zgg>max(eq.zbbbs) | zgg<min(eq.zbbbs) | rhogg>rhomaxgg+dr);
  psigrid = linspace(eq.psimag,eq.psibry,nr)'; % Wb
  Pprime = spline(psigrid,eq.pprime,eq.psizr);
  FFprim = spline(psigrid,eq.ffprim,eq.psizr);
  eq.jphi = (rgg.*Pprime+FFprim/mu0./rgg)/1e6;
  eq.jphi(ivac) = 0;
  eq.cpasma = eq0.cpasma+dip;
  
  % If cccirc was specified then transform dcc back to original coil vector
  if isfield(options,'cccirc')
    dcc2 = eq0.cc;
    for j = 1:length(cccirc)
      dcc2(j) = sign(cccirc(j))*dcc(abs(cccirc(j)));
    end
    dcc = dcc2;
    icx = find(isnan(eq_target.cc));
  end
  eq.cc = eq0.cc+dcc;
  
  % Now converge the changed equilibrium, (gives finite grid and nonlinear effects)
  options.iconstraints = 3;
  options.icx = icx;
  [eq, eqx] = convergeq(eq,tok_data_struct,options,idoplot);
  
  % Add back garbage bbbs elements that were taken out by convergeq
  if length(eq0.rbbbs) > length(eq.rbbbs)
    eq.rbbbs(length(eq0.rbbbs)) = 0;
    eq.zbbbs(length(eq0.rbbbs)) = 0;
  end
