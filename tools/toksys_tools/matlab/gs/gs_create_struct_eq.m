%
%  USAGE:   gs_create_struct_eq
%
%  PURPOSE:  Create a structure eq with an equilibrium description of same Toksys format
%            as is returned by read_mds_eqdsk.m and with additional fields that are used
%            by write_afile.m & write_gfile.m
%
%  INPUTS:   Several equilibrium quantities produced by codes gs_*.m
%
%  OUTPUTS:  eq, structure with equilibrium on Toksys format
%
	
%  METHOD:  
%
%  NOTES:   Grid cells that are partly covered by plasma contain a current
%           j_edge*covered area [A]. Applies to pcurrt only.
%
	
%  VERSION @(#)gs_create_struct_eq.m	1.2 10/16/13
%
%  WRITTEN BY:  Anders Welander  ON	10/08/13
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    clear ds da

if ~exist('psizr_app','var')
gs_eq_analysis
end

    ds.psizr = 'Poloidal flux on grid [Wb]';
    eq.psizr = psizr;

    ds.psizr_app = 'Poloidal flux made by conductor currents [ic;iv]';
    eq.psizr_app = psizr_app;

    ds.psizr_pla = 'Poloidal flux made by plasma currents pcurrt';
    eq.psizr_pla = psizr_pla;

    ds.psipla = 'Current-weighted average flux in plasma';
    eq.psipla = psipla;

    ds.psic = 'Poloidal flux at coils';
    eq.psic = mcc*ic+mcv*iv+mpc'*pcurrt(:);

    ds.psiv = 'Poloidal flux at vessel elements';
    eq.psiv = mcv'*ic+mvv*iv+mpv'*pcurrt(:);

    if exist('frc','var')
    ds.frc = 'Radial forces on coils';
    eq.frc = frc;
    end

    if exist('fzc','var')
    ds.fzc = 'Vertical forces on coils';
    eq.fzc = fzc;
    end

    ds.cc = 'PF coil currents in same units as the template that was used to design this equilibrium';
    eq.cc = piccc*ic;

    ds.ic = 'PF coil currents in Toksys units of terminal Amperes';
    eq.ic = ic;

    ds.iv = 'Vessel currents in Toksys units';
    eq.iv = iv;

    ds.cpasma = 'Plasma current [A]';
    eq.cpasma = cpasma;

    ds.jphi = 'Toroidal current density on grid [MA/m^2]';  
    eq.jphi = (rgg.*pprimezr+ffprimzr/mu0./rgg)/1e6;

    ds.pcurrt = 'Toroidal current within grid rectangles including partially covered cells [A]';  
    eq.pcurrt = pcurrt;

    ds.rcur ='Radius of current centroid [m]';
    eq.rcur = rcur;

    ds.zcur = 'Height of current centroid [m]';
    eq.zcur = zcur;

    ds.nbbbs = 'Number of boundary points';
    eq.nbbbs = nbbbs;

    ds.rbbbs = 'Radius of boundary points [m]';
    eq.rbbbs = rbbbs(nbbbs:-1:1);

    ds.zbbbs = 'Height of boundary points [m]';
    eq.zbbbs = zbbbs(nbbbs:-1:1);
    
    dbbbs = zeros(nbbbs_max,1);
    for j = 1:nbbbs
      dbbbs(j) = sqrt(min((rl-rbbbs(j)).^2+(zl-zbbbs(j)).^2));
    end

    ds.dbbbs = 'Nearest distance to limiter from points [m]';
    eq.dbbbs = dbbbs(nbbbs:-1:1);

    ds.psimag = 'Flux at the magnetic axis [Wb]';
    eq.psimag = psimag;

    ds.rmaxis ='Radius of magnetic axis [m]';
    eq.rmaxis = rmaxis;

    ds.zmaxis = 'Height of magnetic axis [m]';
    eq.zmaxis = zmaxis;

    ds.psibry = 'Flux on the boundary [Wb]';
    eq.psibry = psibry;

    ds.rbdef = 'Radius of the point that defines the boundary [m]';
    eq.rbdef = rbdef;

    ds.zbdef = 'Height of the point that defines the boundary [m]';
    eq.zbdef = zbdef;

    ds.nw = 'Number of radial points on the grid';
    eq.nw = nr;

    ds.nh = 'Number of vertical points on the grid';
    eq.nh = nz;

    ds.rg = 'Radial positions of grid points [m]';
    eq.rg = rg;

    ds.zg = 'Vertical positions of grid points [m]';
    eq.zg = zg;

    ds.dr = 'Radial spacing between grid points [m]';
    eq.dr = dr;

    ds.dz = 'Vertical spacing between grid points [m]';
    eq.dz = dz;

    ds.xdim = 'Width of grid [m]';
    eq.xdim = rg(end)-rg(1);

    ds.zdim = 'Height of grid [m]';
    eq.zdim = zg(end)-zg(1);

    ds.zmid = 'Vertical position of center of grid [m]';
    eq.zmid = (zg(1)+zg(end))/2;

    ds.limitr = 'Number of limiter points in xlim, ylim';
    eq.limitr = nlim;

    ds.xlim = 'Radial locaion of limiter points [m]';
    eq.xlim = Rlim;

    ds.ylim = 'Vertical location of limiter points [m]';
    eq.ylim = Zlim;

    ds.bzero = 'Applied toroidal field at radius rzero [T]';
    eq.bzero = bzero;

    ds.rzero = 'Radius where toroidal field bzero is given [m]';
    eq.rzero = rzero;

    ds.fpol = ['Poloidal current function at ' num2str(nr) ' flux values from axis to boundary [Tm]'];
    eq.fpol = fpol;

    ds.ffprim = 'Derivative of fpol^2/2 w.r.t flux per radian [2*pi*T]';
    eq.ffprim = ffprim;

    ds.pres = ['Pressure at ' num2str(nr) ' flux values from axis to boundary [N/m^2]'];
    eq.pres = pres;

    ds.pprime = 'Derivative of pressure w.r.t. flux per radian [2*pi*N/Wb]';
    eq.pprime = pprime;

    ds.li = 'Normalized plasma inductance';
    eq.li = li;

    ds.betap = 'Plasma pressure / poloidal magnetic pressure';
    eq.betap = betap;

    ds.betat = 'Plasma pressure / toroidal magnetic pressure';
    eq.betat = betat;

    ds.betan = 'Normalized betat';
    eq.betan = betan;
    
%    ds.Vzr = 'Volumes of plasma within grid elements';
%    eq.Vzr = Vzr;

    ds.Vtot = 'Total plasma volume';
    eq.Vtot = Vtot;

    ds.preszr = 'Pressure on grid elements';
    eq.preszr = preszr;

    ds.Wzr = 'Thermal energy within grid elements';
    eq.Wzr = Wzr;

    ds.Wth = 'Total thermal energy';
    eq.Wth = Wth;

    ds.Cl = 'Length of boundary contour';
    eq.Cl = Cl;
    
    ds.Bp2zr = 'square of poloidal field within grid elements';
    eq.Bp2zr = Bp2zr;

    ds.bp2flx = '(mu0*cpasma/Cl)^2';
    eq.bp2flx = bp2flx;

    ds.Bp2V = 'sum(sum(Vzr.*Bp2zr))';
    eq.Bp2V = Bp2V;
    
    ds.psirz = '-psizr''/2/pi, efitviewer plots contours with this variable';
    eq.psirz = -psizr'/2/pi;

    ds.ssibry = '-psibry/2/pi, efitviewer plots contours with this variable';
    eq.ssibry = -psibry/2/pi;

    ds.ssimag = '-psimag/2/pi, efitviewer plots contours with this variable';
    eq.ssimag = -psimag/2/pi;
    
    if exist('nulls','var')
      ds.nulls = 'Information about x-points within the limiter';
      eq.nulls = nulls;
    end

    if exist('rx1','var')
      ds.rx1 = 'Radius of closest x point found below plasma [m]';
      eq.rx1 = rx1;
    end

    if exist('zx1','var')
      ds.zx1 = 'Height of closest x point found below plasma [m]';
      eq.zx1 = zx1;
    end

    if exist('psix1','var')
      ds.psix1 = 'Flux at x1 [Wb]';
      eq.psix1 = psix1;
    end

    if exist('rx2','var')
      ds.rx2 = 'Radius of closest x point found above plasma [m]';
      eq.rx2 = rx2;
    end

    if exist('zx2','var')
      ds.zx2 = 'Height of closest x point found above plasma [m]';
      eq.zx2 = zx2;
    end

    if exist('psix2','var')
      ds.psix2 = 'Flux at x2 [Wb]';
      eq.psix2 = psix2;
    end
    
    rsurf = (min(rbbbs)+max(rbbbs))/2;
    ds.rsurf = 'Radius of geometric center of boundary [m]';
    eq.rsurf = rsurf;

    ds.zsurf = 'Height of geometric center of boundary [m]';
    eq.zsurf = (min(zbbbs)+max(zbbbs))/2;

    ds.aminor = 'Minor radius of boundary [m]';
    eq.aminor = aminor;

    ds.bminor = 'Minor height of boundary [m]';
    eq.bminor = (max(zbbbs)-min(zbbbs))/2;

    ds.elong = 'Elongation of boundary';
    eq.elong = eq.bminor/eq.aminor;

    [rout, k] = max(rbbbs);
    zout = zbbbs(k);
    [ztop, k] = max(zbbbs);
    rtop = rbbbs(k);
    [rinn, k] = min(rbbbs);
    zinn = zbbbs(k);
    [zbot, k] = min(zbbbs);
    rbot = rbbbs(k);

    doutl = (rsurf-rbot)/aminor;
    doutu = (rsurf-rtop)/aminor;

    ds.doutl = 'Lower triangularity';
    eq.doutl = doutl;

    ds.doutu = 'Upper triangularity';
    eq.doutu = doutu;

    da.uday = 'Date the equilibrium was created';
    uday = datestr(date,'dd-mmm-yyyy');
    uday = uday([1:7 10:11]);
    ea.uday = uday;

    da.mfvers = 'Something about version';
    ea.mfvers = datestr(now,23);

    da.ktime = 'Number of time slices';
    ea.ktime = 1;

    da.jflag = 'Some efit flag';
    ea.jflag = 1;

    da.lflag = 'Some efit flag';
    ea.lflag = 0;
    
if 0
    thbdef = angle(rbdef-rmaxis + 1i*(zbdef-zmaxis));
    da.limloc = 'Shape descriptor, possible values are: IN , OUT, TOP, BOT, DN , SNT, SNB, MAR';
    if lim
      if     thbdef*180/pi > -45 & thbdef*180/pi <=  45
	ea.limloc = 'OUT';
      elseif thbdef*180/pi >  45 & thbdef*180/pi <= 135
	ea.limloc = 'TOP';
      elseif thbdef*180/pi >  -135 & thbdef*180/pi <= -45
	ea.limloc = 'BOT';
      else
	ea.limloc = 'IN ';
      end
    elseif ix1 & ix2
      if     drsep > +0.0025
	ea.limloc = 'SNT';
      elseif drsep < -0.0025
	ea.limloc = 'SNB';
      else
	ea.limloc = 'DN ';
      end
    elseif ix1
      ea.limloc = 'SNB';
    elseif ix2
      ea.limloc = 'SNT';
    end
end

    da.oleft = 'Smallest distance to limiter from boundary sector 135<theta<225 [cm]';
    ea.oleft = inf;
    for j = 1:nbbbs
      if thbbbs(j)*180/pi > 135 | thbbbs(j)*180/pi < -135
        ea.oleft = min(100*dbbbs(j),ea.oleft);
      end
    end

    da.oright = 'Smallest distance to limiter from boundary sector -45<theta<45 [cm]';
    ea.oright = inf;
    for j = 1:nbbbs
      if thbbbs(j)*180/pi > -45 & thbbbs(j)*180/pi < 45
        ea.oright = min(100*dbbbs(j),ea.oright);
      end
    end

    da.otop = 'Smallest distance to limiter from boundary sector 45<theta<135 [cm]';
    ea.otop = inf;
    for j = 1:nbbbs
      if thbbbs(j)*180/pi > 45 & thbbbs(j)*180/pi < 135
        ea.otop = min(100*dbbbs(j),ea.otop);
      end
    end

    da.obott = 'Smallest distance to limiter from boundary sector -135<theta<-45 [cm]';
    ea.obott = inf;
    for j = 1:nbbbs
      if thbbbs(j)*180/pi > -135 & thbbbs(j)*180/pi < -45
        ea.obott = min(100*dbbbs(j),ea.obott);
      end
    end

    if exist('rx1','var')
    da.rseps1 = 'Radial position of lower x-point [cm]';
    ea.rseps1 = eq.rx1*100;
    end

    if exist('zx1','var')
    da.zseps1 = 'Vertical position of lower x-point [cm]';
    ea.zseps1 = eq.zx1*100;
    end

    if exist('rx2','var')
    da.rseps2 = 'Radial position of upper x-point [cm]';
    ea.rseps2 = eq.rx2*100;
    end

    if exist('zx2','var')
    da.zseps2 = 'Vertical position of upper x-point [cm]';
    ea.zseps2 = eq.zx2*100;
    end

    da.betapd = 'Diamagnetic betap';
    ea.betapd = betap;

    da.wplasmd = 'Diamagnetic thermal energy [J]';
    ea.wplasmd = Wth;

    da.wplasm = 'Thermal energy [J]';
    ea.wplasm = Wth;

if 0
    if x1exists & x1inside
      dminlx = 100*sqrt(min((rl-rx1).^2+(zl-zx1).^2));
    else
      dminlx = 1e3;
    end

    da.dminlx = 'Nearest distance from lower x-point to limiter [cm]';
    ea.dminlx = dminlx;

    if x2exists & x2inside
      dminux = 100*sqrt(min((rl-rx2).^2+(zl-zx2).^2));
    else
      dminux = 1e3;
    end

    da.dminux = 'Nearest distance from upper x-point to limiter [cm]';
    ea.dminux = dminux;
end

    if exist('fluxexp','var')
      ds.fluxexp = 'Flux expansion';
      eq.fluxexp = fluxexp;
    end

    if exist('rfluxexp','var')
      ds.rfluxexp = 'Radius where flux expansion calculated';
      eq.rfluxexp = rfluxexp;
    end

    if exist('zfluxexp','var')
      ds.zfluxexp = 'Height where flux expansion calculated';
      eq.zfluxexp = zfluxexp;
    end

    if exist('rfluxexp2','var')
      ds.rfluxexp2 = 'Radius for outer contour at flux expansion point';
      eq.rfluxexp2 = rfluxexp2;
    end

    if exist('zfluxexp2','var')
      ds.zfluxexp2 = 'Height for outer contour at flux expansion point';
      eq.zfluxexp2 = zfluxexp2;
    end

    if exist('dfluxexp','var')
      ds.dfluxexp = 'midplane distance in calculation of flux expansion';
      eq.dfluxexp = dfluxexp;
    end
    
    if exist('rcont','var')
      ds.rcont = 'Radius of flux surface points';
      eq.rcont = rcont;
    end
    if exist('zcont','var')
      ds.zcont = 'Height of flux surface points';
      eq.zcont = zcont;
    end
    if exist('V','var')
      ds.V = 'Volume contained within flux surfaces';
      eq.V = V;
    end
    if exist('A','var')
      ds.A = 'Area contained within flux surfaces';
      eq.A = A;
    end
    if exist('L','var')
      ds.L = 'Integral of 1/R within flux surfaces';
      eq.L = L;
    end
    if exist('W','var')
      ds.W = 'Thermal energy within flux surfaces';
      eq.W = W;
    end
    if exist('B','var')
      ds.B = 'Profile of betap';
      eq.B = B;
    end
    if exist('I','var')
      ds.I = 'Toroidal current within flux surfaces';
      eq.I = I;
    end
    if exist('T','var')
      ds.T = 'Toroidal flux within flux surfaces';
      eq.T = T;
    end
    if exist('jtav','var')
      ds.jtav = 'Contour-averaged current density (dI/dA)';
      eq.jtav = jtav;
    end
    if exist('rhot','var')
      ds.rhot = 'Square-root of normalized toroidal flux, sqrt(T/T(nr))';
      eq.rhot = rhot;
    end
    if exist('qpsi','var')
      ds.qpsi = 'q values (dT/dpsi)';
      eq.qpsi = qpsi;
    end

    ds.psibar = 'Normalized poloidal flux';
    eq.psibar = linspace(0,1,nr)';

    ds.rhop = 'Square-root of normalized poloidal flux';
    eq.rhop = sqrt(linspace(0,1,nr)');

    ds.adata = 'Variables needed by write_afile to write an a-file';
    eq.adata = ea;

    fields_from_template = strvcat('fcturn','turnfc','fcid','ecturn','ecid','idx_efit_to_tok');
    fields_from_template_descr = strvcat(...
      'Multiplier of EFIT input currents, if = 1 brsp is in terminal [A]',...
      'Multiplier used in generation of EFIT greens tables; if = 1 brsp is in A-turn',...
      'Coil identifier used in construction of greens tables',...
      'Individual turns in ecoil(s)',...
      'Coil identifier used in construction of greens tables',...
      'Indices to reorder coils from efit to Toksys order');
    for j = 1:size(fields_from_template,1)
      s = deblank(fields_from_template(j,:));
      ss = deblank(fields_from_template_descr(j,:));
      if exist(s,'var') & ~isempty(eval(s))
	ds = setfield(ds,s,ss);
	eq = setfield(eq,s,eval(s));
      end
    end
    if ~isfield(eq,'ecid')
      eq.ecid = [];
    end

    if exist('psizr_app','var') & ~isempty(psizr_app)
      ds.psizr_app = 'flux generated by coil and vessel currents: [ic; iv]';
      eq.psizr_app = psizr_app;
    end

    if exist('psizr_pla','var') & ~isempty(psizr_pla)
      ds.psizr_pla = 'flux generated by plasma currents: pcurrt';
      eq.psizr_pla = psizr_pla;
    end

    ds.psizr_err = 'flux error, psizr_err = psizr-psizr_pla-psizr_app;';
    eq.psizr_err = psizr_err;

    if exist('fl','var') & ~isempty(fl)
      ds.fl = 'Calculated flux loop signals';
      eq.fl = fl;
    end

    if exist('bp','var') & ~isempty(bp)
      ds.bp = 'Calculated magnetic probe signals';
      eq.bp = bp;
    end

    if exist('rog','var') & ~isempty(rog)
      ds.rog = 'Calculated rogowski signals';
      eq.rog = rog;
    end
    
    if exist('archive_design_data','var') & archive_design_data
      ds.design_data = 'equilibrium was created by gsdesign, specifics are here';
      dsx = [];
      if exist('targets','var')
        eq.design_data.targets = targets;
        dsx.targets = 'design targets';
      end
      if exist('weights','var')
        eq.design_data.weights = weights;
        dsx.weights = 'weights for design targets';
      end
      if exist('locks','var')
        eq.design_data.locks = locks;
        dsx.locks = 'locked quantities (nan = not locked)';
      end
      if exist('limits','var')
        eq.design_data.limits = limits;
        dsx.limits = 'lower and upper limits';
      end
      if exist('ev','var')
        eq.design_data.ev = ev;
        dsx.ev = 'error vector = weights.*errors in targets';
      end
      if exist('strev','var')
        eq.design_data.strev = strev;
        dsx.strev = 'descriptions of ev in tex format';
	for i = 1:size(strev,1)
	  str = deblank(strev(i,:));
	  str = strrep(str,'_','');
	  str = strrep(str,'\','');
	  str = strrep(str,'partial','d');
	  eq.design_data.evnames{i} = str;
	end
        dsx.evnames = 'descriptions of ev in plain text';
      end
      eq.design_data.descriptions = dsx;
    end
    
    if exist('y','var') & ~isempty(y)
      ds.y = 'The outputs y from the gs function';
      eq.y = y;
      if exist('index_in_y','var') & ~isempty(index_in_y)
        ds.index_in_y = 'indices in y for different diagnostics';
        eq.index_in_y = index_in_y;
      end
    end
    
    eq.descriptions = ds;
