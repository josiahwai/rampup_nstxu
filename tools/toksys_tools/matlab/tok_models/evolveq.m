function [y, e, eq, eqx] = evolveq(u,C,tok,eq0)
%
%  USAGE:   [y, e, eq, eqx] = evolveq(u,C,tok,eq0)
%
%  PURPOSE: Evolve an equilibrium. 
%           Initialize by calling with all arguments:  y = evolveq(u,C,tok,eq0)
%           Evolve by calling with new inputs in u:    y = evolveq(u)
%
%  INPUTS: u,  structure with inputs that affect the equilibrium:
%               cc0t: coil currents such that psizr_app = mpc*cc0t
%               vc0t: vessel currents such that psizr_app = mpv*vc0t
%                 ip: plasma current
%                 li: normalized inductance
%              betap: poloidal beta
%              If an input is missing then no change of the quantity is made.
%          C,  output matrix such that y = C*[cc0t; vc0t; pcurrt], 
%              where pcurrt is plasma currents within grid elements [A]
%              Example: C = [[mcc mcv mpc']; [mcv' mvv mpv']] outputs flux at conductors
%        tok,  Toksys description of the tokamak
%        eq0,  initial equilibrium
%
%  OUTPUTS:   y, outputs defined by C
%             e, structure with all persistent variables
%            eq, equilibrium (returning eq increases execution time)
%           eqx. last analyzed equilibrium and some simulation parameters
%	
%  METHOD:  The equilibrium and its response matrix are persistent variables.
%           When the routine is called, the equilibrium is updated by the linear response to the
%             changes that have occurred in u, and a new output, y is calculated.
%           Large changes can be made in one step, evolveq estimates magnitude of nonlinear
%            response and divides large changes into several smaller internally.
%           When uncertainties in the linear response have accrued, a correction is made.
%           When the response matrix becomes inaccurate a new is calculated.

%  NOTES:   eq describes an equilibrium in the same way as an EFIT. However, additional fields
%           are added by evolveq. Also, grid cells that are partly covered by plasma contain a
%           j_edge*covered area [A]. This affects pcurrt only (not jphi). psizr_pla = Mgg*pcurrt
	
%  VERSION @(#)evolveq.m	1.2 02/08/13
%
%  WRITTEN BY:  Anders Welander  ON	11/19/12
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Remember C
  persistent Cmat
  
  % Remember the tokamak
  persistent nr nz ngg rg zg dr dz rgg zgg mpc mpv mps mpp xpj Rlim Zlim nlim rl zl nl wls iil ilimgg klimgg piccc nc nv
  
  % Remember helping variables
  persistent mu0 twopi R13 mx neighbors nkn psikn T0 T sp0 sf0 sg0 c0 c1 c2 c3 iknotg psibar 
  persistent v1r v2r v3r v4r v1z v2z v3z v4z ii ni igridedge dP dixdx disdx Ua Ub
  persistent p1_0 p2_0 p3_0 p4_0 f1_0 f2_0 f3_0 f4_0 g1_0 g2_0 g3_0 g4_0
  
  % Remember reference equilibrium
  persistent psimag0 psibry0 jpla0 iplasma0 iknot0 djdpsi0 sqsum_djdpsi0

  % Remember predicted changes since last reference equilibrium
  persistent jpla0_pred jpla0_est
  
  % Remember last analyzed equilibrium
  persistent cpasma rmaxis zmaxis psimag r0 z0 psibry x1exists x1inside rx1 zx1 x2exists x2inside rx2 zx2 rzero bzero Pprime FFprim P
  persistent iplasma rbbbs zbbbs zbot ztop psipla li betap jphi pcurrt ii0 w0 iia wa ilimited ix1 ix2 nbbbs 
  persistent fpol ffprim pres pprime psix1 psix2 thbbbs Wth
  persistent psizr_err Vzr Bp2zr Cl Bp2V Vtot bp2flx
  
  % Remember present equilibrium
  persistent ic iv psizr sp sf yprev psibarzr
  
  % Remember plasma response
  persistent Ti dpsizrdx dpsimagdx dpsibrydx djphidpsizr dpcurrtdx dcpasmadx dlidx dbetapdx
  persistent dydx dxdu drbbbsdx dzbbbsdx djphidx djpla0dx dj1 djx
  
  % For testing
  persistent old_value new_value predicted_change dips dlis dbetaps derrors Cls Bp2Vs Vtots bp2flxs r0s Wths
  
  % For EFITs
  persistent fcturn turnfc fcid ecturn ecid idx_efit_to_tok 
  
  
  
  
%                                                                                                                                              
%                                                                                                                                            
%     IIIIIIIIII                  iiii          tttt            iiii                    lllllll   iiii                                       
%     I::::::::I                 i::::i      ttt:::t           i::::i                   l:::::l  i::::i                                      
%     I::::::::I                  iiii       t:::::t            iiii                    l:::::l   iiii                                       
%     II::::::II                             t:::::t                                    l:::::l                                              
%       I::::Innnn  nnnnnnnn    iiiiiiittttttt:::::ttttttt    iiiiiii   aaaaaaaaaaaaa    l::::l iiiiiii zzzzzzzzzzzzzzzzz    eeeeeeeeeeee    
%       I::::In:::nn::::::::nn  i:::::it:::::::::::::::::t    i:::::i   a::::::::::::a   l::::l i:::::i z:::::::::::::::z  ee::::::::::::ee  
%       I::::In::::::::::::::nn  i::::it:::::::::::::::::t     i::::i   aaaaaaaaa:::::a  l::::l  i::::i z::::::::::::::z  e::::::eeeee:::::ee
%       I::::Inn:::::::::::::::n i::::itttttt:::::::tttttt     i::::i            a::::a  l::::l  i::::i zzzzzzzz::::::z  e::::::e     e:::::e
%       I::::I  n:::::nnnn:::::n i::::i      t:::::t           i::::i     aaaaaaa:::::a  l::::l  i::::i       z::::::z   e:::::::eeeee::::::e
%       I::::I  n::::n    n::::n i::::i      t:::::t           i::::i   aa::::::::::::a  l::::l  i::::i      z::::::z    e:::::::::::::::::e 
%       I::::I  n::::n    n::::n i::::i      t:::::t           i::::i  a::::aaaa::::::a  l::::l  i::::i     z::::::z     e::::::eeeeeeeeeee  
%       I::::I  n::::n    n::::n i::::i      t:::::t    tttttt i::::i a::::a    a:::::a  l::::l  i::::i    z::::::z      e:::::::e           
%     II::::::IIn::::n    n::::ni::::::i     t::::::tttt:::::ti::::::ia::::a    a:::::a l::::::li::::::i  z::::::zzzzzzzze::::::::e          
%     I::::::::In::::n    n::::ni::::::i     tt::::::::::::::ti::::::ia:::::aaaa::::::a l::::::li::::::i z::::::::::::::z e::::::::eeeeeeee  
%     I::::::::In::::n    n::::ni::::::i       tt:::::::::::tti::::::i a::::::::::aa:::al::::::li::::::iz:::::::::::::::z  ee:::::::::::::e  
%     IIIIIIIIIInnnnnn    nnnnnniiiiiiii         ttttttttttt  iiiiiiii  aaaaaaaaaa  aaaalllllllliiiiiiiizzzzzzzzzzzzzzzzz    eeeeeeeeeeeeee  
%                                                                                                                                            
%                                                                                                                                             
    
  % Initial loop - create variables to assist with boundary tracing and calculation of response matrix
  
  if nargin == 4
    
    Cmat = C;
    
    struct_to_ws(tok);
    struct_to_ws(eq0);

    ngg = nr*nz;

    mu0 = 0.4e-6*pi;
    twopi = 2*pi;
    R13 = (1+i*sqrt(3))^2/4; % Increases angle in complex plane by 2*pi/3	

    mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2; % Used for interpolation
    neighbors = reshape([-1-nz -1 -1+nz -1+2*nz;-nz 0 nz 2*nz;1-nz 1 1+nz 1+2*nz;2-nz 2 2+nz 2+2*nz],1,16);

    nkn = 1;
    psikn = linspace(0,1,nkn+1); % psibar for knot junctions, default is equally spaced values
    % c1, c2, c3 are determined based on boundary conditions.
    % Hence these are new linear combination of old coeffs for each new interval
    x = [1 1 1; 1 2 3; 0 2 6];
    e = [1 2 3; 0 1 2; 0 0 1];
    c0 = zeros(nkn,nkn+2); c0(nkn,nkn+(0:2)) = -1;
    c1 = zeros(nkn,nkn+2); c1(nkn,nkn) = 1;
    c2 = zeros(nkn,nkn+2); c2(nkn,nkn+1) = 1;
    c3 = zeros(nkn,nkn+2); c3(nkn,nkn+2) = 1;
    for j = nkn-1:-1:1
      c = inv(psikn(j+1).^e.*x)*[1;0;0];
      c0(j,j) = 1;
      c1(j,j:nkn+2) = c(1)*(c0(j+1,j:nkn+2)-c0(j,j:nkn+2));
      c1(j,j:nkn+2) = c1(j,j:nkn+2)+c1(j+1,j:nkn+2);
      c2(j,j:nkn+2) = c(2)*(c0(j+1,j:nkn+2)-c0(j,j:nkn+2));
      c2(j,j:nkn+2) = c2(j,j:nkn+2)+c2(j+1,j:nkn+2);
      c3(j,j:nkn+2) = c(3)*(c0(j+1,j:nkn+2)-c0(j,j:nkn+2));
      c3(j,j:nkn+2) = c3(j,j:nkn+2)+c3(j+1,j:nkn+2);
    end   
    
    psibar = linspace(0,1,nr)';
    iknotg = psibar*0+1;
    f = psibar*ones(1,nkn+2);
    for j = nkn:-1:1
      iknotg(psibar <= psikn(j+1),1) = j;
    end
    sp0 = pinv([c1(iknotg,:)+2*c2(iknotg,:).*f+3*c3(iknotg,:).*f.^2])*pprime*(psibry-psimag)/twopi;
    sf0 = pinv([c1(iknotg,:)+2*c2(iknotg,:).*f+3*c3(iknotg,:).*f.^2])*ffprim*(psibry-psimag)/twopi;
    sg0 = pinv([c1(iknotg,:)+2*c2(iknotg,:).*f+3*c3(iknotg,:).*f.^2])*((1-psibar).*ffprim)*(psibry-psimag)/twopi;
    p0_0 = c0*sp0;    p1_0 = c1*sp0;    p2_0 = c2*sp0;    p3_0 = c3*sp0;
    f0_0 = c0*sf0;    f1_0 = c1*sf0;    f2_0 = c2*sf0;    f3_0 = c3*sf0;
    g0_0 = c0*sg0;    g1_0 = c1*sg0;    g2_0 = c2*sg0;    g3_0 = c3*sg0;
    sp = sp0;
    sf = sf0;
    sg = sg0;
    
    psi50 = zeros(nz,nr);
    Vzr = zeros(nz,nr); % Amount of plasma volume within grid element
    dP = zeros(nz,nr); % Changes in pressure
    P = zeros(nz,nr); % Pressure on grid
    Pprime = zeros(nz,nr); % pprime on grid
    FFprim = zeros(nz,nr); % ffprim on grid
    Ua = zeros(ngg,1); % Will hold RHS of equations for modified axis flux due to axis displacement
    Ub = zeros(ngg,1); % Will hold RHS of equations for modified boundary flux due to displacement
    disdx = [eye(nc+nv) zeros(nc+nv,6)]; % How conductor currents respond to x-variables   
    dixdx = [eye(nc+nv+3) zeros(nc+nv+3,3)]; % How conductor currents,a,b,e respond to x-variables

    % The following are used to solve for tz when tr=0.5
    v1z = reshape(mx'*[0 0 0 1]'*[1 0.5 0.25 0.125]*mx,1,16);
    v2z = reshape(mx'*[0 0 1 0]'*[1 0.5 0.25 0.125]*mx,1,16);
    v3z = reshape(mx'*[0 1 0 0]'*[1 0.5 0.25 0.125]*mx,1,16);
    v4z = reshape(mx'*[1 0 0 0]'*[1 0.5 0.25 0.125]*mx,1,16);

    % The following are used to solve for tr when tz=0.5
    v1r = reshape(mx'*[1 0.5 0.25 0.125]'*[0 0 0 1]*mx,1,16);
    v2r = reshape(mx'*[1 0.5 0.25 0.125]'*[0 0 1 0]*mx,1,16);
    v3r = reshape(mx'*[1 0.5 0.25 0.125]'*[0 1 0 0]*mx,1,16);
    v4r = reshape(mx'*[1 0.5 0.25 0.125]'*[1 0 0 0]*mx,1,16);

    x1exists = false;
    x1inside = false;
    x2exists = false;
    x2inside = false;

    if size(limdata,2)>size(limdata,1)
      limdata=limdata';
    end
    if min(limdata(:,1))<min(limdata(:,2))
      limdata=limdata(:,[2 1]);
    end
    if limdata(1,1) == limdata(end,1) & limdata(1,2) == limdata(end,2)
      Rlim = limdata(:,1);
      Zlim = limdata(:,2);
    else
      Rlim = [limdata(:,1); limdata(1,1)];
      Zlim = [limdata(:,2); limdata(1,2)];
    end
    nlim = length(Rlim);
    rl = Rlim(1);
    zl = Zlim(1);
    dcmin = min(dr,dz);
    for j = 1:nlim-1
      d = sqrt((Rlim(j+1)-Rlim(j))^2+(Zlim(j+1)-Zlim(j))^2);
      n = 1+round(2*d/dcmin);
      rl = [rl(1:end-1) linspace(Rlim(j),Rlim(j+1),n)];
      zl = [zl(1:end-1) linspace(Zlim(j),Zlim(j+1),n)];
    end
    nl = length(rl);
    iil = zeros(nl,16);
    wls = zeros(nl,16);
    for j = 1:length(rl)
      kr0 = min(nr-3,max(1,floor((rl(j)-rg(1))/dr)));
      kz1 = min(nz-2,max(2,ceil((zl(j)-zg(1))/dz)));
      k = kr0*nz+kz1;
      iil(j,:) = k+neighbors;
      tr = (rl(j)-rgg(k))/dr;
      tz = (zl(j)-zgg(k))/dz;
      wls(j,:) = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
    end

    % ilimgg is a matrix for the grid containing indices of limdata that intereect rg+[0 dr], zg=[0 dz]
    ilimgg = zeros(nz,nr);
    ilimgg( 1,:) = -1;
    ilimgg(nz,:) = -1;
    ilimgg(:, 1) = -1;
    ilimgg(:,nr) = -1;
    for j = 1:nlim-1
      kra = (Rlim(j+0)-rg(1))/dr;
      kza = (Zlim(j+0)-zg(1))/dz;
      krb = (Rlim(j+1)-rg(1))/dr;
      kzb = (Zlim(j+1)-zg(1))/dz;
      sr = sign(krb-kra);
      sz = sign(kzb-kza);
      % Calculate intersections of Rlim([j j+1]), Zlim([j j+1]) with rg
      if kra == krb
	fr = [];
      else
	fr = interp1([kra krb],[0 1],floor(kra):sr:floor(krb));
	fr = fr(fr>0 & fr<1);
      end
      % Calculate intersections of Rlim([j j+1]), Zlim([j j+1]) with zg
      if kza == kzb
	fz = [];
      else
	fz = interp1([kza kzb],[0 1],floor(kza):sz:floor(kzb));
	fz = fz(fz>0 & fz<1);
      end
      f = [0 fr fz];
      [f, indf] = sort(f);
      dk = 0*f+sr*nz;
      dk(1) = 0;
      dk(indf>length(fr)+1) = sz;
      k = nz*floor(kra)+floor(kza)+1+cumsum(dk);
      ilimgg(k) = j;
    end
    not_found_all_vac_gg = true;
    while not_found_all_vac_gg
      not_found_all_vac_gg = false;
      for j = 2:nz-1
	for k = 2:nr-1
	  if ilimgg(j,k) == 0 & (ilimgg(j-1,k)==-1 | ilimgg(j+1,k)==-1 | ilimgg(j,k-1)==-1 | ilimgg(j,k+1)==-1)
	    ilimgg(j,k) = -1;
	    not_found_all_vac_gg = true;
	  end
	end
      end      
    end
    klimgg = find(ilimgg > -1);

    % Equations will be for delstar(psi) for inner grid points, so r.h.s. must be changed
    ii = (2:nz-1)'*ones(1,nr-2)+ones(nz-2,1)*(1:nr-2)*nz; ii = ii(:); % Indices to the interior of the grid
    iedges = [2:nz-1 ngg+(2-nz:-1) nz+1:nz:ngg-nz 2*nz:nz:ngg-nz]';
    icorners = [1 nz ngg-nz+1 ngg]';
    isides = [iedges; icorners];
    mps = [mpc mpv];
    for j = 1:size(mps,2)
      mps(ii,j) = mps(ii,j)*(-2/dr^2-2/dz^2)+mps(ii-1,j)/dz^2+mps(ii+1,j)/dz^2+...
	mps(ii-nz,j).*(1/dr^2+1/2./rgg(ii)/dr)+mps(ii+nz,j).*(1/dr^2-1/2./rgg(ii)/dr);
    end

    % Finite difference equations and mutuals equations for the edge
    T0 = sparse(ii,ii, -2/dr^2-2/dz^2,           ngg, ngg)+... % Equations for change in delstar: 
	 sparse(ii,ii-nz,1/dr^2+1/2./rgg(ii)/dr, ngg, ngg)+... % delstar (dpsi_pla+dpsi_app) becomes:
	 sparse(ii,ii+nz,1/dr^2-1/2./rgg(ii)/dr, ngg, ngg)+... %  (delstar+mu0*R*djdpsi) dpsi = delstar(dpsi_app)
	 sparse(ii,ii-1, 1/dz^2,                                              ngg, ngg)+...
	 sparse(ii,ii+1, 1/dz^2,                                              ngg, ngg)+...
	 sparse(isides, isides, 1,                                            ngg, ngg);    % The 1 in (1-mpp*di/dpsi)

    % xpj gives flux at edge given jphi [MA/m2] on grid.
    ni = 2*(nz+nr)-4; % There are ni grid points at the edge
    % They are ordered as bottom row, top row, then bottom-top, bottom-top from left to right
    igridedge = [1:nz:ngg nz:nz:ngg]; % define bottom row and top row now and the rest in loop below
    xpj = zeros(ni,ngg);
    xpj(1:nr,:) = 1e6*dr*dz*mpp';
    xpj(nr+1:2*nr,:) = 1e6*dr*dz*mpp((1:ngg)+mod(nz-(1:ngg),nz)-mod((1:ngg)-1,nz),:)';
    k = length(igridedge);
    for j = 2:nz-1
      izshift = (1+abs(j-(1:nz))-(1:nz))'*ones(1,nr);
      k = k+1;
      igridedge(k) = j;
      xpj(k,:) = 1e6*dr*dz*mpp((1:ngg)+izshift(1:ngg),1)';
      k = k+1;
      igridedge(k) = j+(nr-1)*nz;
      xpj(k,:) = 1e6*dr*dz*mpp((1:ngg)+izshift(1:ngg),end)';
    end

    % iccc makes a toksys ic from an efit cc
    ncc = length(cc);
    dum = eq0;
    for j = 1:ncc
      dum.cc = eq0.cc*0;
      dum.cc(j) = 1;
      equil_I = cc_efit_to_tok(tok,dum);
      iccc(:,j) = equil_I.cc0t;
    end
    piccc = pinv(iccc);
    equil_I = cc_efit_to_tok(tok,eq0);
    ic = equil_I.cc0t;
    iv = equil_I.vc0t;
    
    dj1 = inf(ngg,1);
    djx = inf(ngg,1);
    iplasma0 = (1:ngg)';
    psibarzr = inf(ngg,1);
    djdpsi0 = zeros(ngg,1);
    sqsum_djdpsi0 = 1;
    jphi_estimate = inf(nz,nr);
    jphi0 = zeros(nz,nr);
    sum_djphi = zeros(nz,nr);
    
    jpla0_est = ones(ngg,1);
    jpla0_pred = inf(ngg,1);
    
    dcpasmadx  = zeros(1,nc+nv+3);
    dlidx      = zeros(1,nc+nv+3);
    dbetapdx   = zeros(1,nc+nv+3);
    drmaxisdx  = zeros(1,nc+nv+3);
    dzmaxisdx  = zeros(1,nc+nv+3);
    dr0dx      = zeros(1,nc+nv+3);
    dz0dx      = zeros(1,nc+nv+3);
    djphidx    = zeros(ngg,nc+nv+3);
    djphitotdx = zeros(ngg,nc+nv+3);
    dpcurrtdx  = zeros(ngg,nc+nv+3);
    
  end % Done initializing
  
  
  
  
  not_done = true;
  iteration_counter = 0;
  
  while not_done
  
    % Decide whether new plasma response matrix is needed
    djdpsi = (dj1+djx.*psibarzr(iplasma0))/(psibry-psimag)^2;
    sqsum_diff_djdpsi = sum((djdpsi-djdpsi0).^2);
    rel_lin_err = sqrt(sqsum_diff_djdpsi/sqsum_djdpsi0);
    if rel_lin_err > 0.05
      new_Ti = true;
    else
      new_Ti = false;
    end
  
    % Decide whether correction of flux error is needed
    rel_totlin_err = sqrt(sum((jpla0_est-jpla0_pred).^2)/sum(jpla0_est.^2));
    % jphi0 is the jphi of the reference equilibrium
    % sum_djphi is the total linearly predicted change since the reference equilibrium
    % jphi_estimate is a more precise estimate of the current jphi than jphi0+sum_djphi
  




    %     
    %     
    %                  AAA                                                 lllllll                                                              
    %                 A:::A                                                l:::::l                                                              
    %                A:::::A                                               l:::::l                                                              
    %               A:::::::A                                              l:::::l                                                              
    %              A:::::::::A         nnnn  nnnnnnnn      aaaaaaaaaaaaa    l::::lyyyyyyy           yyyyyyyzzzzzzzzzzzzzzzzz    eeeeeeeeeeee    
    %             A:::::A:::::A        n:::nn::::::::nn    a::::::::::::a   l::::l y:::::y         y:::::y z:::::::::::::::z  ee::::::::::::ee  
    %            A:::::A A:::::A       n::::::::::::::nn   aaaaaaaaa:::::a  l::::l  y:::::y       y:::::y  z::::::::::::::z  e::::::eeeee:::::ee
    %           A:::::A   A:::::A      nn:::::::::::::::n           a::::a  l::::l   y:::::y     y:::::y   zzzzzzzz::::::z  e::::::e     e:::::e
    %          A:::::A     A:::::A       n:::::nnnn:::::n    aaaaaaa:::::a  l::::l    y:::::y   y:::::y          z::::::z   e:::::::eeeee::::::e
    %         A:::::AAAAAAAAA:::::A      n::::n    n::::n  aa::::::::::::a  l::::l     y:::::y y:::::y          z::::::z    e:::::::::::::::::e 
    %        A:::::::::::::::::::::A     n::::n    n::::n a::::aaaa::::::a  l::::l      y:::::y:::::y          z::::::z     e::::::eeeeeeeeeee  
    %       A:::::AAAAAAAAAAAAA:::::A    n::::n    n::::na::::a    a:::::a  l::::l       y:::::::::y          z::::::z      e:::::::e           
    %      A:::::A             A:::::A   n::::n    n::::na::::a    a:::::a l::::::l       y:::::::y          z::::::zzzzzzzze::::::::e          
    %     A:::::A               A:::::A  n::::n    n::::na:::::aaaa::::::a l::::::l        y:::::y          z::::::::::::::z e::::::::eeeeeeee  
    %    A:::::A                 A:::::A n::::n    n::::n a::::::::::aa:::al::::::l       y:::::y          z:::::::::::::::z  ee:::::::::::::e  
    %   AAAAAAA                   AAAAAAAnnnnnn    nnnnnn  aaaaaaaaaa  aaaallllllll      y:::::y           zzzzzzzzzzzzzzzzz    eeeeeeeeeeeeee  
    %                                                                                   y:::::y                                                 
    %                                                                                  y:::::y                                                  
    %                                                                                 y:::::y                                                   
    %                                                                                y:::::y                                                    
    %                                                                               yyyyyyy                                                     
    %                                                                                                                                           
    %                                                                                                                                           
    %                                                                                                                                                                             
    %                                                                                    bbbbbbbb                                                                                 
    %                                                               iiii  lllllll   iiii b::::::b                                 iiii                                            
    %                                                              i::::i l:::::l  i::::ib::::::b                                i::::i                                           
    %                                                               iiii  l:::::l   iiii b::::::b                                 iiii                                            
    %                                                                     l:::::l         b:::::b                                                                                 
    %       eeeeeeeeeeee       qqqqqqqqq   qqqqquuuuuu    uuuuuu  iiiiiii  l::::l iiiiiii b:::::bbbbbbbbb    rrrrr   rrrrrrrrr  iiiiiii uuuuuu    uuuuuu     mmmmmmm    mmmmmmm   
    %     ee::::::::::::ee    q:::::::::qqq::::qu::::u    u::::u  i:::::i  l::::l i:::::i b::::::::::::::bb  r::::rrr:::::::::r i:::::i u::::u    u::::u   mm:::::::m  m:::::::mm 
    %    e::::::eeeee:::::ee q:::::::::::::::::qu::::u    u::::u   i::::i  l::::l  i::::i b::::::::::::::::b r:::::::::::::::::r i::::i u::::u    u::::u  m::::::::::mm::::::::::m
    %   e::::::e     e:::::eq::::::qqqqq::::::qqu::::u    u::::u   i::::i  l::::l  i::::i b:::::bbbbb:::::::brr::::::rrrrr::::::ri::::i u::::u    u::::u  m::::::::::::::::::::::m
    %   e:::::::eeeee::::::eq:::::q     q:::::q u::::u    u::::u   i::::i  l::::l  i::::i b:::::b    b::::::b r:::::r     r:::::ri::::i u::::u    u::::u  m:::::mmm::::::mmm:::::m
    %   e:::::::::::::::::e q:::::q     q:::::q u::::u    u::::u   i::::i  l::::l  i::::i b:::::b     b:::::b r:::::r     rrrrrrri::::i u::::u    u::::u  m::::m   m::::m   m::::m
    %   e::::::eeeeeeeeeee  q:::::q     q:::::q u::::u    u::::u   i::::i  l::::l  i::::i b:::::b     b:::::b r:::::r            i::::i u::::u    u::::u  m::::m   m::::m   m::::m
    %   e:::::::e           q::::::q    q:::::q u:::::uuuu:::::u   i::::i  l::::l  i::::i b:::::b     b:::::b r:::::r            i::::i u:::::uuuu:::::u  m::::m   m::::m   m::::m
    %   e::::::::e          q:::::::qqqqq:::::q u:::::::::::::::uui::::::il::::::li::::::ib:::::bbbbbb::::::b r:::::r           i::::::iu:::::::::::::::uum::::m   m::::m   m::::m
    %    e::::::::eeeeeeee   q::::::::::::::::q  u:::::::::::::::ui::::::il::::::li::::::ib::::::::::::::::b  r:::::r           i::::::i u:::::::::::::::um::::m   m::::m   m::::m
    %     ee:::::::::::::e    qq::::::::::::::q   uu::::::::uu:::ui::::::il::::::li::::::ib:::::::::::::::b   r:::::r           i::::::i  uu::::::::uu:::um::::m   m::::m   m::::m
    %       eeeeeeeeeeeeee      qqqqqqqq::::::q     uuuuuuuu  uuuuiiiiiiiilllllllliiiiiiiibbbbbbbbbbbbbbbb    rrrrrrr           iiiiiiii    uuuuuuuu  uuuummmmmm   mmmmmm   mmmmmm
    %                                   q:::::q                                                                                                                                   
    %                                   q:::::q                                                                                                                                   
    %                                  q:::::::q                                                                                                                                  
    %                                  q:::::::q                                                                                                                                  
    %                                  q:::::::q                                                                                                                                  
    %                                  qqqqqqqqq                                                                                                                                  
    %                                                                                                                                                                             
    %     

    % calculate axis and boundary positions and then flux error to remove it
    
    new_eq_analysis = rel_totlin_err > 1e-3 | new_Ti;
    
    if new_eq_analysis

      % Magnetic axis
      twopirbrzmax = inf; % Maximum of 2*pi*R*Br, 2*pi*R*Bz at calculated null position
      j = 9;
      while j>0 & twopirbrzmax>1e-10
	kr0 = min(nr-3,max(1,floor((rmaxis-rg(1))/dr)));
	kz1 = min(nz-2,max(2,ceil((zmaxis-zg(1))/dz)));
	k = kr0*nz+kz1;
	iia = k+neighbors; % iia indexes 16 points around magnetic axis
	pp = psizr(iia);
	tr = (rmaxis-rgg(k))/dr;
	tz = (zmaxis-zgg(k))/dz;
	wa = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
	war = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
	waz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
	warr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]/dr^2*mx,1,16);
	wazz = reshape(([0 0 2 6*tz]/dz^2*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
	warz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
	ashift_Ba = -inv([pp*warr' pp*warz';pp*warz' pp*wazz']);
	c = ashift_Ba*[pp*war'; pp*waz'];
	rmaxis = rmaxis+c(1);
	zmaxis = zmaxis+c(2);
	if norm(c) > dr
	  [dum,k] = min(psizr(iplasma)/(psibry-psimag));
	  rmaxis = rgg(iplasma(k));
	  zmaxis = zgg(iplasma(k));
	end
	j = j-1;
	twopirbrzmax = max(abs([pp*war' pp*waz']));
      end
      psimag = wa*pp';

      psizr_r = [zeros(nz,1) psizr(:,3:end)-psizr(:,1:end-2) zeros(nz,1)]/2/dr;
      psizr_z = [zeros(1,nr);psizr(3:end,:)-psizr(1:end-2,:);zeros(1,nr)]/2/dz;
      
      % Check for LOWER NULL
      if ~(x1exists & x1inside)
	j = find(zbbbs < zmaxis);
	[dum, k] = min(interp2(rg,zg,psizr_r.^2+psizr_z.^2,rbbbs(j),zbbbs(j),'spline'));
	rx1 = rbbbs(j(k));
	zx1 = zbbbs(j(k));
      end

      twopirbrzmax = inf; % Maximum of 2*pi*R*Br, 2*pi*R*Bz at calculated null position
      j = 9;
      while j > 0 & twopirbrzmax > 1e-10 % Try zooming in on x-point with Newton Rhapson.
	j = j-1;
	% Find indices and weights for grid points around the x-point
	kr0 = min(nr-3,max(1,floor((rx1-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
	kz1 = min(nz-2,max(2,ceil((zx1-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
	k = kr0*nz+kz1;
	iix1 = k+neighbors; % indexes 16 points around x point
	pp = psizr(iix1)';
	tr = (rx1-rgg(k))/dr;
	tz = (zx1-zgg(k))/dz;
	wx1 = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
	wx1r = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
	wx1z = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
	wx1rr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]/dr^2*mx,1,16);
	wx1zz = reshape(([0 0 2 6*tz]/dz^2*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
	wx1rz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
	x1shift_Bx = -inv([wx1rr*pp wx1rz*pp;wx1rz*pp wx1zz*pp]);
	cx1 = x1shift_Bx*[wx1r*pp; wx1z*pp];
	rx1 = rx1+cx1(1);
	zx1 = zx1+cx1(2);
	twopirbrzmax = max(abs([wx1r*pp wx1z*pp]));
      end
      psix1 = wx1*pp;

      if twopirbrzmax > 1e-10 || zx1 > zmaxis-dz
	x1exists = false; % No x-point found
      else
	x1exists = true; % x-point found
      end
      if x1exists % X-point found. Now check if it is inside the limiter
	if ilimgg(k) == 0 % In this case it is clearly inside
	  x1inside = true; % Flag that x1 is a real (and important) x-point	  
	elseif ilimgg(k) == -1 % That means clearly outside
	  x1inside = false;
	else
	  x1inside = true;
	  for j = max(1,ilimgg(k)-5):min(nlim-1,ilimgg(k)+5);
            mlinex = [Rlim(j+1)-Rlim(j) rmaxis-rx1; Zlim(j+1)-Zlim(j) zmaxis-zx1];
	    if rcond(mlinex) > 0
	      kg = inv(mlinex)*[rmaxis-Rlim(j); zmaxis-Zlim(j)];
	      if kg(1)>=0 & kg(1)<=1 & kg(2)>0 & kg(2)<1 % x1 outside limiter
		x1inside = false;
	      end
	    end
	  end
	end
      end

      if x1exists & x1inside
	zbot = zx1; % Min z for where plasma can be
      else
	psix1 = inf/(psibry-psimag);
	zbot = min(zbbbs)-dz;
      end
      psibarx1 = (psix1-psimag)/(psibry-psimag);


      % Check for UPPER NULL
      if ~(x2exists & x2inside)
	j = find(zbbbs > zmaxis);
	[dum, k] = min(interp2(rg,zg,psizr_r.^2+psizr_z.^2,rbbbs(j),zbbbs(j),'spline'));
	rx2 = rbbbs(j(k));
	zx2 = zbbbs(j(k));
      end

      twopirbrzmax = inf; % Maximum of 2*pi*R*Br, 2*pi*R*Bz at calculated null position
      j = 9;
      while j > 0 & twopirbrzmax > 1e-10 % Try zooming in on x-point with Newton Rhapson.
	j = j-1;
	% Find indices and weights for grid points around the x-point
	kr0 = min(nr-3,max(1,floor((rx2-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
	kz1 = min(nz-2,max(2,ceil((zx2-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
	k = kr0*nz+kz1;
	iix2 = k+neighbors; % indexes 16 points around x point
	pp = psizr(iix2)';
	tr = (rx2-rgg(k))/dr;
	tz = (zx2-zgg(k))/dz;
	wx2 = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
	wx2r = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
	wx2z = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
	wx2rr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]/dr^2*mx,1,16);
	wx2zz = reshape(([0 0 2 6*tz]/dz^2*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
	wx2rz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
	x2shift_Bx = -inv([wx2rr*pp wx2rz*pp;wx2rz*pp wx2zz*pp]);
	cx2 = x2shift_Bx*[wx2r*pp; wx2z*pp];
	rx2 = rx2+cx2(1);
	zx2 = zx2+cx2(2);
	twopirbrzmax = max(abs([wx2r*pp wx2z*pp]));
      end
      psix2 = wx2*pp;

      if twopirbrzmax > 1e-10 || zx2 < zmaxis+dz
	x2exists = false; % No x-point found
      else
	x2exists = true; % x-point found
      end
      if x2exists % X-point found. Now check if it is inside the limiter
	if ilimgg(k) == 0 % In this case it is clearly inside
	  x2inside = true; % Flag that x2 is a real (and important) x-point	  
	elseif ilimgg(k) == -1 % That means clearly outside
	  x2inside = false;
	else
	  x2inside = true;
	  for j = max(1,ilimgg(k)-5):min(nlim-1,ilimgg(k)+5);
            mlinex = [Rlim(j+1)-Rlim(j) rmaxis-rx2; Zlim(j+1)-Zlim(j) zmaxis-zx2];
	    if rcond(mlinex) > 0
	      kg = inv(mlinex)*[rmaxis-Rlim(j); zmaxis-Zlim(j)];
	      if kg(1)>=0 & kg(1)<=1 & kg(2)>0 & kg(2)<1 % x2 outside limiter
		x2inside = false;
	      end
	    end
	  end
	end
      end

      if x2exists & x2inside
	ztop = zx2; % Max z for where plasma can be
      else
	psix2 = inf/(psibry-psimag);
	ztop = max(zbbbs)+dz;
      end
      psibarx2 = (psix2-psimag)/(psibry-psimag);


      % Find TOUCH POINT candidate
      for j = 1:nl % Calculate for all of limiter to find touch and strike points
	psilim(j) = wls(j,:)*psizr(iil(j,:))';
      end

      k = find(rl>rg(1) & rl<rg(end) & zl>zg(1) & zl<zg(end));
      if x1exists & x1inside
	k = intersect(k,find(zl > zx1));
      end
      if x2exists & x2inside
	k = intersect(k,find(zl < zx2));
      end

      [dum, l] = min(psilim(k)/(psibry-psimag)); k1 = k(l);
      rlim = rl(k1); zlim = zl(k1);
      k2 = k1+1;
      if k2 > length(rl)
        k2 = 2;
      end
      ulim = [rl(k2)-rl(k1) , zl(k1)-zl(k2)];
      nulim = norm(ulim);
      ulim = ulim/nulim;
      kr0 = min(nr-3,max(1,floor((rlim-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
      kz1 = min(nz-2,max(2,ceil((zlim-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
      k = kr0*nz+kz1;
      iilim = k+neighbors;
      tr = (rlim-rgg(k))/dr;
      tz = (zlim-zgg(k))/dz;
      wl = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      wld = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]*mx*ulim(1)/dr+...
                     ([0 1 2*tz 3*tz^2]*mx)'*ulim(2)/dz*[1 tr tr^2 tr^3]*mx,1,16);
      wlb = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]*mx*ulim(1)^2/dr^2+...
        	   ([0 0 2 6*tz]*mx)'*ulim(2)^2/dz^2*[1 tr tr^2 tr^3]*mx+...
		 2*([0 1 2*tz 3*tz^2]*mx)'*ulim(2)/dz*[0 1 2*tr 3*tr^2]*mx*ulim(1)/dr,1,16);
      psilimbis = psizr(iilim)*wlb';
      dlim = -psizr(iilim)*wld'/psilimbis*ulim;
      for j = 1:3
	if norm(dlim)<nulim
	  rlim = rlim+dlim(1); zlim = zlim+dlim(2);
	  kr0 = min(nr-3,max(1,floor((rlim-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
	  kz1 = min(nz-2,max(2,ceil((zlim-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
	  k = kr0*nz+kz1; iilim = k+neighbors;
	  tr = (rlim-rgg(k))/dr; tz = (zlim-zgg(k))/dz;
	  wl = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
	  wld = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]*mx*ulim(1)/dr+...
                	 ([0 1 2*tz 3*tz^2]*mx)'*ulim(2)/dz*[1 tr tr^2 tr^3]*mx,1,16);
	  wlb = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]*mx*ulim(1)^2/dr^2+...
                       ([0 0 2 6*tz]*mx)'*ulim(2)^2/dz^2*[1 tr tr^2 tr^3]*mx+...
		     2*([0 1 2*tz 3*tz^2]*mx)'*ulim(2)/dz*[0 1 2*tr 3*tr^2]*mx*ulim(1)/dr,1,16);
	  psilimbis = psizr(iilim)*wlb';
	  dlim = -psizr(iilim)*wld'/psilimbis*ulim;
	end
      end
      psitouch = psizr(iilim)*wl';
      psibartouch = (psitouch-psimag)/(psibry-psimag);

      % Which point defines the boundary?
      ilimited = 0; ix1 = 0; ix2 = 0; % Flags for what defines boundary
      if psibartouch < min(psibarx1,psibarx2)
	ilimited = 1;
	r0 = rlim;
	z0 = zlim;
	w0 = wl;
	ii0 = iilim;
	psibry = psitouch;
      else
	if psibarx1 < psibarx2
	  ix1 = 1;
	  r0 = rx1;
	  z0 = zx1;
	  w0 = wx1;
	  ii0 = iix1;
	  psibry = psix1;
	else
	  ix2 = 1;
	  r0 = rx2;
	  z0 = zx2;
	  w0 = wx2;
	  ii0 = iix2;
	  psibry = psix2;
	end
      end
      psibarzr = (psizr-psimag)/(psibry-psimag);

      % Trace boundary
      % When tz = 0: psi = [1 tr tr^2 tr^3]*mx*psizr(k+[-nz 0 nz 2*nz]')
      k = find((psizr(klimgg)-psibry).*(psizr(klimgg+nz)-psibry) < 0);
      k = klimgg(k);
      d = psizr(k)-psibry;
      c = (psizr(k+nz)-psizr(k-nz))/2;
      b = (2*psizr(k-nz)-5*psizr(k)+4*psizr(k+nz)-psizr(k+2*nz))/2;
      a = (-psizr(k-nz)+3*psizr(k)-3*psizr(k+nz)+psizr(k+2*nz))/2;
      aa = 2*b.^3-9*a.*b.*c+27*a.^2.*d;
      bb = sqrt(aa.^2-4*(b.^2-3*a.*c).^3);
      xx = d./(psizr(k)-psizr(k+nz));

      q1 = ((aa+bb)/2).^(1/3);
      q2 = ((aa-bb)/2).^(1/3);

      j = find(imag(q1+q2)>1e-6);	q2(j) = R13*q2(j);
      j = find(imag(q1+q2)>1e-6);	q2(j) = R13*q2(j);

      x1 = -(b+real(q1+q2))/3./a;
      x2 = -(b-(1+i*sqrt(3))/2*q1-(1-i*sqrt(3))/2*q2)/3./a;
      x3 = -(b-(1-i*sqrt(3))/2*q1-(1+i*sqrt(3))/2*q2)/3./a;

      for j = 1:length(xx)
	if     abs(imag(x1(j))) < 1e-6 & x1(j) > 0 & x1(j) < 1
	  xx(j) = abs(x1(j));
	elseif abs(imag(x2(j))) < 1e-6 & x2(j) > 0 & x2(j) < 1
	  xx(j) = abs(x2(j));
	elseif abs(imag(x3(j))) < 1e-6 & x3(j) > 0 & x3(j) < 1
	  xx(j) = abs(x3(j));
	end
      end

      rbbbs = rgg(k)+xx*dr;
      zbbbs = zgg(k);
      ib = k;
      flag_for_r_or_z_edge = k*0+1; % 1 for tz=0 and 2 for tr=0

      % When tr = 0: psi = [1 tz tz^2 tz^3]*mx*psizr(k+[-1 0 1 2]')
      k = find((psizr(klimgg)-psibry).*(psizr(klimgg+1)-psibry) < 0);
      k = klimgg(k);
      d = psizr(k)-psibry;
      c = (psizr(k+1)-psizr(k-1))/2;
      b = (2*psizr(k-1)-5*psizr(k)+4*psizr(k+1)-psizr(k+2))/2;
      a = (-psizr(k-1)+3*psizr(k)-3*psizr(k+1)+psizr(k+2))/2;
      aa = 2*b.^3-9*a.*b.*c+27*a.^2.*d;
      bb = sqrt(aa.^2-4*(b.^2-3*a.*c).^3);
      xx = d./(psizr(k)-psizr(k+1));

      q1 = ((aa+bb)/2).^(1/3);
      q2 = ((aa-bb)/2).^(1/3);

      j = find(imag(q1+q2)>1e-6);	q2(j) = R13*q2(j);
      j = find(imag(q1+q2)>1e-6);	q2(j) = R13*q2(j);

      x1 = -(b+real(q1+q2))/3./a;
      x2 = -(b-(1+i*sqrt(3))/2*q1-(1-i*sqrt(3))/2*q2)/3./a;
      x3 = -(b-(1-i*sqrt(3))/2*q1-(1+i*sqrt(3))/2*q2)/3./a;

      for j = 1:length(xx)
	if     abs(imag(x1(j))) < 1e-6 & x1(j) > 0 & x1(j) < 1
	  xx(j) = abs(x1(j));
	elseif abs(imag(x2(j))) < 1e-6 & x2(j) > 0 & x2(j) < 1
	  xx(j) = abs(x2(j));
	elseif abs(imag(x3(j))) < 1e-6 & x3(j) > 0 & x3(j) < 1
	  xx(j) = abs(x3(j));
	end
      end

      rbbbs = [rbbbs; rgg(k)];
      zbbbs = [zbbbs; zgg(k)+xx*dz];
      ib = [ib; k];
      flag_for_r_or_z_edge = [flag_for_r_or_z_edge; k*0+2];
      
      k = find(zbbbs > zbot & zbbbs < ztop);
      rbbbs = rbbbs(k);
      zbbbs = zbbbs(k);
      ib = ib(k);
      flag_for_r_or_z_edge = flag_for_r_or_z_edge(k);

      rbbbs = [rbbbs; r0];
      zbbbs = [zbbbs; z0];
      ib = [ib; ii0(6)];
      flag_for_r_or_z_edge(end+1) = 0;

      [thbbbs, k] = sort(angle(rbbbs-rmaxis+i*(zbbbs-zmaxis)));
      rbbbs = rbbbs(k);
      zbbbs = zbbbs(k);
      ib = ib(k);
      ibdef = find(k==max(k));
      flag_for_r_or_z_edge = flag_for_r_or_z_edge(k);

      rbbbs(end+1) = rbbbs(1);
      zbbbs(end+1) = zbbbs(1);
      thbbbs(end+1) = thbbbs(1)+2*pi;
      ib(end+1) = ib(1);
      flag_for_r_or_z_edge(end+1) = flag_for_r_or_z_edge(1);
      rhobbbs = sqrt((rbbbs-rmaxis).^2+(zbbbs-zmaxis).^2);
      nbbbs = length(rbbbs);

      for j = 1:nbbbs
	k = ib(j);
	iip = k+neighbors';
	tr = (rbbbs(j)-rgg(k))/dr;
	tz = (zbbbs(j)-zgg(k))/dz;
	if flag_for_r_or_z_edge(j) == 1 % Save a bit of time by only calculating the needed derivative
	  dpsibbbsdr(j,1) = reshape(mx'*[1 tz tz^2 tz^3]'*[0 1 2*tr 3*tr^2]/dr*mx,1,16)*psizr(iip);
	else
	  dpsibbbsdz(j,1) = reshape(mx'*[0 1 2*tz 3*tz^2]'/dz*[1 tr tr^2 tr^3]*mx,1,16)*psizr(iip);
	end
	wb(j,1:16) = reshape(mx'*[1 tz tz^2 tz^3]'*[1 tr tr^2 tr^3]*mx,1,16);
      end

      [rbbbsmin, irbbbsmin] = min(rbbbs(1:nbbbs-1));
      [rbbbsmax, irbbbsmax] = max(rbbbs(1:nbbbs-1));
      [zbbbsmin, izbbbsmin] = min(zbbbs(1:nbbbs-1));
      [zbbbsmax, izbbbsmax] = max(zbbbs(1:nbbbs-1));
      thgg = angle(rgg-rmaxis+i*(zgg-zmaxis));
      rhbg = spline(thbbbs(1:nbbbs-1),rhobbbs(1:nbbbs-1),thgg);
      rhgg = sqrt((rgg-rmaxis).^2+(zgg-zmaxis).^2);

      % Find grid points in the plasma, and which are new, and which were removed
      iplasma = find(rhgg<=rhbg+dr & psibarzr<=1 & zgg<=max(zbbbs) & zgg>=min(zbbbs));
      ivac = setdiff(1:ngg,iplasma);
      np = length(iplasma);
      
      x = psibarzr(iplasma);

      iknot = ceil(interp1(psikn,0:nkn,x)); % indices to spline knots
      iknot(x<=0) = 1;  % In case there is a problem with finding the axis
      iknot(x>1) = nkn; % In case there is a problem with finding the boundary defining point

      psigrid = linspace(psimag,psibry,nr)';

      p0 = c0*sp; p1 = c1*sp; p2 = c2*sp; p3 = c3*sp;
      pres = p0(iknotg)+p1(iknotg).*psibar+p2(iknotg).*psibar.^2+p3(iknotg).*psibar.^3; pres(end) = 0;
      pprime = (p1(iknotg)+p2(iknotg)*2.*psibar+p3(iknotg)*3.*psibar.^2)*twopi/(psibry-psimag);

      f0 = c0*sf; f1 = c1*sf; f2 = c2*sf; f3 = c3*sf;
      halffpolsquared = f0(iknotg)+f1(iknotg).*psibar+f2(iknotg).*psibar.^2+f3(iknotg).*psibar.^3+rzero^2*bzero^2/2;
      k = find(halffpolsquared < 0); halffpolsquared(k) = 0;
      fpol = sqrt(2*halffpolsquared);
      ffprim = (f1(iknotg)+f2(iknotg)*2.*psibar+f3(iknotg)*3.*psibar.^2)*twopi/(psibry-psimag);

      Pprime(iplasma) = (p1(iknot)+p2(iknot)*2.*x+p3(iknot)*3.*x.^2)*twopi/(psibry-psimag);
      Pprime(ivac) = 0;
      
      FFprim(iplasma) = (f1(iknot)+f2(iknot)*2.*x+f3(iknot)*3.*x.^2)*twopi/(psibry-psimag);
      FFprim(ivac) = 0;
      
      jphi = (rgg.*Pprime+FFprim/mu0./rgg)/1e6;

      % Calculate what fraction of cells are covered by plasma and how that responds to dpsizr
      irbbbs = (rbbbs-rg(1))/dr+1; % rbbbs in grid number units
      izbbbs = (zbbbs-zg(1))/dz+1; % zbbbs in grid number units
      nxe = 0; % Index of last calculated crossing of boundary into new cell
      a = []; b = []; c = []; d = [];
      for j = 1:nbbbs-1
	sz = sign(izbbbs(j+1)-izbbbs(j));
	kz = round(izbbbs(j))+sz*0.5:sz:round(izbbbs(j+1)); % z edges in grid number units
	lkz = length(kz);
	ikz = 0; % Next index in kz to process
	sr = sign(irbbbs(j+1)-irbbbs(j));
	kr = round(irbbbs(j))+sr*0.5:sr:round(irbbbs(j+1)); % r edges in grid number units
	lkr = length(kr);
	ikr = 0; % Next index in kr to process
	while ikr<lkr | ikz<lkz % There are edges to cross
	  nxe = nxe+1;
	  if ikr < lkr
	    f2r = (kr(ikr+1)-irbbbs(j))/(irbbbs(j+1)-irbbbs(j)); % Fraction of next re,ze that comes from rbbbs(j+1),zbbbs(j+1)
	  else
	    f2r = inf;
	  end
	  if ikz < lkz
	    f2z = (kz(ikz+1)-izbbbs(j))/(izbbbs(j+1)-izbbbs(j)); % Fraction of next re,ze that comes from rbbbs(j+1),zbbbs(j+1)
	  else
	    f2z = inf;
	  end
	  if f2r < f2z % Next crossing is over r edge
	    ikr = ikr+1;
	    f2e(nxe) = f2r;
	    fe(nxe) = 2*sr; % Flag that we crossed over r edge
	  else % Next crossing is over z edge
	    ikz = ikz+1;
	    f2e(nxe) = f2z;
	    fe(nxe) = sz; % Flag that we crossed over z edge
	  end
	  R = (1-f2e(nxe))*rbbbs(j)+f2e(nxe)*rbbbs(j+1);
	  Z = (1-f2e(nxe))*zbbbs(j)+f2e(nxe)*zbbbs(j+1);
	  % calculate correction of the coordinate R,Z
	  kr0 = min(nr-3,max(1,floor((R-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
	  kz1 = min(nz-2,max(2,ceil((Z-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
	  k = kr0*nz+kz1;
	  ie(nxe) = k; % Remember indices for grid points that are inside & below points (re,ze)
	  iip = k+neighbors;
	  pp = psizr(iip)';
	  if f2r < f2z % Next crossing is over r edge, so find exactly right z position
	    a(nxe,1) = v1z*pp;
	    b(nxe,1) = v2z*pp;
	    c(nxe,1) = v3z*pp;
	    d(nxe,1) = v4z*pp-psibry;
	    re(nxe,1) = R;
	  else % Next crossing is over z edge, so find exactly right r position
	    a(nxe,1) = v1r*pp;
	    b(nxe,1) = v2r*pp;
	    c(nxe,1) = v3r*pp;
	    d(nxe,1) = v4r*pp-psibry;
	    ze(nxe,1) = Z;
	  end
	end
      end
      ie = ie(1:nxe);

      aa = 2*b.^3-9*a.*b.*c+27*a.^2.*d;
      bb = sqrt(aa.^2-4*(b.^2-3*a.*c).^3);

      q1 = ((aa+bb)/2).^(1/3);
      q2 = ((aa-bb)/2).^(1/3);

      j = find(imag(q1+q2)>1e-6);	q2(j) = R13*q2(j);
      j = find(imag(q1+q2)>1e-6);	q2(j) = R13*q2(j);

      x1 = -(b+real(q1+q2))/3./a;
      x2 = -(b-(1+i*sqrt(3))/2*q1-(1-i*sqrt(3))/2*q2)/3./a;
      x3 = -(b-(1-i*sqrt(3))/2*q1-(1+i*sqrt(3))/2*q2)/3./a;

    xx = ones(nxe,1)/2;;
    for j = 1:length(xx)
      if     abs(imag(x1(j))) < 1e-6 & x1(j) > -0.5 & x1(j) < 1.5
	xx(j) = abs(x1(j));
      elseif abs(imag(x2(j))) < 1e-6 & x2(j) > -0.5 & x2(j) < 1.5
	xx(j) = abs(x2(j));
      elseif abs(imag(x3(j))) < 1e-6 & x3(j) > -0.5 & x3(j) < 1.5
	xx(j) = abs(x3(j));
      end
    end

      for j = 1:nxe
	pp = psizr(ie(j)+neighbors');
	if abs(fe(j)) == 1
	  tz = 0.5;
	  tr = xx(j);
	  re(j) = rgg(ie(j))+tr*dr;
	else
          tr = 0.5;
          tz = xx(j);
	  ze(j) = zgg(ie(j))+tz*dz;
	end
	dpsiedr(j,1) = reshape(mx'*[1 tz tz^2 tz^3]'*[0 1 2*tr 3*tr^2]/dr*mx,1,16)*pp;
	dpsiedz(j,1) = reshape(mx'*[0 1 2*tz 3*tz^2]'/dz*[1 tr tr^2 tr^3]*mx,1,16)*pp;
	we(j,1:16) = reshape(mx'*[1 tz tz^2 tz^3]'*[1 tr tr^2 tr^3]*mx,1,16);
      end

      % Check for neighboring positions in rec,zec where the order in poloidal angle switched because of the corrections
      dthe = diff(angle(re(1:nxe)-rmaxis+i*(ze(1:nxe)-zmaxis)));
      k = find(abs(dthe) < pi/2 & dthe < 0);
      fec = fe;
      for j = k'
	dum = [re(j) ze(j) fe(j) ie(j) dpsiedr(j) dpsiedz(j) we(j,:)];
	re(j) = re(j+1);
	ze(j) = ze(j+1);
	fe(j) = fe(j+1);
	ie(j) = ie(j+1);
	dpsiedr(j) = dpsiedr(j+1);
	dpsiedz(j) = dpsiedz(j+1);
	we(j,:) = we(j+1,:);
	re(j+1) = dum(1);
	ze(j+1) = dum(2);
	fe(j+1) = dum(3);
	ie(j+1) = dum(4);
	dpsiedr(j+1) = dum(5);
	dpsiedz(j+1) = dum(6);
	we(j+1,:) = dum(7:22);
      end

      % At this point we have: 
      %   re, ze, coordinates where boundary goes into new grid rectangle
      %   fe, flags +1=entered up (+z) into rectange, -1=down (-z) into, +2=outward (+r) into, -2=inward (-r) into grid rectange

      % Calculate index into collapsed grid of cells that are exited at re, ze
      if fe(nxe) == -1
	icell(nxe) = 1+ceil((ze(nxe)-zg(1))/dz)+round((re(nxe)-rg(1))/dr)*nz;
      elseif fe(nxe) == 2
	icell(nxe) = 1+round((ze(nxe)-zg(1))/dz)+floor((re(nxe)-rg(1))/dr)*nz;
      elseif fe(nxe) == +1
	icell(nxe) = ceil((ze(nxe)-zg(1))/dz)+round((re(nxe)-rg(1))/dr)*nz;
      elseif fe(nxe) == -2
	icell(nxe) = 1+round((ze(nxe)-zg(1))/dz)+ceil((re(nxe)-rg(1))/dr)*nz;
      end
      jm = nxe;
      for j = 1:nxe
	if     fe(jm) == +1
	  icell(j) = icell(jm)+1;
	elseif fe(jm) == -1
	  icell(j) = icell(jm)-1;
	elseif fe(jm) == +2
	  icell(j) = icell(jm)+nz;
	elseif fe(jm) == -2
	  icell(j) = icell(jm)-nz;
	end
	jm = j;
      end
      icell = icell(1:nxe);
      
      % Calculate coordinates rc, zc of grid rectangle corners that are inside the plasma and share r or z with re, ze
      for j = 1:nxe
	if     fe(j) == +1
	  rc(j) = rg(1+round((re(j)-rg(1))/dr))-dr/2;
	  zc(j) = ze(j);
	elseif fe(j) == -1
	  rc(j) = rg(1+round((re(j)-rg(1))/dr))+dr/2;
	  zc(j) = ze(j);
	elseif fe(j) == +2
	  zc(j) = zg(1+round((ze(j)-zg(1))/dz))+dz/2;
	  rc(j) = re(j);
	elseif fe(j) == -2
	  zc(j) = zg(1+round((ze(j)-zg(1))/dz))-dz/2;
	  rc(j) = re(j);
	end
      end

      % Calculate Ae = areas of grid cells that are covered by plasma, and their derivatives w.r.t. re, ze
      % Here dAedre1 is dAe(j) w.r.t. re(j-1) and dAedre2 is dAe(j) w.r.t. dre(j)
      Ae = zeros(1,nxe);
      dAedre1 = zeros(1,nxe);
      dAedre2 = zeros(1,nxe);
      dAedze1 = zeros(1,nxe);
      dAedze2 = zeros(1,nxe);
      for j = 1:nxe
	if j > 1, jm = j-1; else jm = nxe; end
	if abs(fe(j)) == 1 & abs(fe(jm)) == 1 % Coming in and out over a z edge
	  Ae(j) = abs(re(j)+re(jm)-rc(j)-rc(jm))*dz/2;
	  dAedre1(j) = dz*sign(re(jm)-rc(jm))/2;
	  dAedre2(j) = dz*sign(re(j)-rc(j))/2;
	 elseif abs(fe(j)) == 2 & abs(fe(jm)) == 2 % Coming in and out over an r edge
	  Ae(j) = abs(ze(j)+ze(jm)-zc(j)-zc(jm))*dr/2;
	  dAedze1(j) = dr*sign(ze(jm)-zc(jm))/2;
	  dAedze2(j) = dr*sign(ze(j)-zc(j))/2;
	else % Coming in over r or z edge and leaving over the other kind
	  if norm([rc(j)-rc(jm) zc(j)-zc(jm)]) < 1e-12 % If they share the same corner, we have a triangle of plasma
	    Ae(j) = norm([re(j)-rc(j) ze(j)-zc(j)])*norm([re(jm)-rc(jm) ze(jm)-zc(jm)])/2;
	    if abs(fe(j)) == 1 % re(j), ze(j) is on a z edge, so re(jm), ze(jm) is an r edge
	      dAedze1(j) = abs(re(j)-rc(j))*sign(ze(jm)-zc(jm))/2;
	      dAedre2(j) = abs(ze(jm)-zc(jm))*sign(re(j)-rc(j))/2;
	    else % re(j), ze(j) is on an r edge
	      dAedre1(j) = abs(ze(j)-zc(j))*sign(re(jm)-rc(j))/2;
	      dAedze2(j) = abs(re(jm)-rc(j))*sign(ze(j)-zc(j))/2;
	    end
	  else % We have a triangle of vacuum, calculate the common corner rcvac, zcvac
	    if abs(fe(j)) == 1 % Shift corner j in R direction
              rcvac = rc(j)+dr*sign(re(j)-rc(j));
	      zcvac = zc(j);
	      dAedze1(j) = -abs(re(j)-rcvac)*sign(ze(jm)-zcvac)/2;
	      dAedre2(j) = -abs(ze(jm)-zcvac)*sign(re(j)-rcvac)/2;
	    else % Shift corner j in R direction
              rcvac = rc(j);
	      zcvac = zc(j)+dz*sign(ze(j)-zc(j));
	      dAedre1(j) = -abs(ze(j)-zcvac)*sign(re(jm)-rcvac)/2;
	      dAedze2(j) = -abs(re(jm)-rcvac)*sign(ze(j)-zcvac)/2;
	    end
	    Ae(j) = dr*dz-norm([re(j)-rcvac ze(j)-zcvac])*norm([re(jm)-rcvac ze(jm)-zcvac])/2;
	  end
	end
      end

      % Current density modification due to partial coverage of grid cell [MA/m2]
      jphixe = zeros(nz,nr);

      % How current is modified within rectangles due to boundary motion due to change of psizr
      djphixedpsizr = zeros(nxe,ngg); % [MA/m2]/Wb
      djphixeds = zeros(nxe,3);
      
      Ac = Ae; % Will become area correction due to partial cell coverage, so negative for cells where jphi ~= 0

      for j = 1:nxe
	if j > 1, jm = j-1; else jm = nxe; end
	R = rgg(icell(j));
	je = (R*pprime(end)+ffprim(end)/mu0/R);
	if jphi(icell(j)) ~= 0 % Subtract current according to fraction of vacuum within cell
	  Ac(j) = Ae(j)-dr*dz;
	  jj = (Ae(j)-dr*dz)/dr/dz*je/1e6; % How to modify jphi
          djphixeds(j,1) = (Ae(j)-dr*dz)/dr/dz/1e6*(p1_0(end)+2*p2_0(end)+3*p3_0(end))*twopi/(psibry-psimag)*R;
          djphixeds(j,2) = (Ae(j)-dr*dz)/dr/dz/1e6*(f1_0(end)+2*f2_0(end)+3*f3_0(end))*twopi/(psibry-psimag)/mu0/R;
          djphixeds(j,3) = (Ae(j)-dr*dz)/dr/dz/1e6*(g1_0(end)+2*g2_0(end)+3*g3_0(end))*twopi/(psibry-psimag)/mu0/R;
	else
	  jj = Ae(j)/dr/dz*je/1e6; % How to modify jphi
          djphixeds(j,1) = Ae(j)/dr/dz/1e6*(p1_0(end)+2*p2_0(end)+3*p3_0(end))*twopi/(psibry-psimag)*R;
          djphixeds(j,2) = Ae(j)/dr/dz/1e6*(f1_0(end)+2*f2_0(end)+3*f3_0(end))*twopi/(psibry-psimag)/mu0/R;
          djphixeds(j,3) = Ae(j)/dr/dz/1e6*(g1_0(end)+2*g2_0(end)+3*g3_0(end))*twopi/(psibry-psimag)/mu0/R;
	end
	jphixe(icell(j)) = jj;
	iie = ie(j)+neighbors;
	djphixedpsizr(j,iie) = -je*(dAedre2(j)/dpsiedr(j)+dAedze2(j)/dpsiedz(j))*we(j,:)/dr/dz/1e6;
	iie = ie(jm)+neighbors;
	djphixedpsizr(j,iie) = djphixedpsizr(j,iie)-je*(dAedre1(j)/dpsiedr(jm)+dAedze1(j)/dpsiedz(jm))*we(jm,:)/dr/dz/1e6;
	djphixedpsizr(j,ii0) = djphixedpsizr(j,ii0)+jj*w0/(psimag-psibry)+...
	  je*(dAedre1(j)/dpsiedr(jm)+dAedze1(j)/dpsiedz(jm)+dAedre2(j)/dpsiedr(j)+dAedze2(j)/dpsiedz(j))*w0/dr/dz/1e6;
	djphixedpsizr(j,iia) = djphixedpsizr(j,iia)-jj*wa/(psimag-psibry);
      end

      jphitot = jphi+jphixe;
      pcurrt = 1e6*dr*dz*jphitot;
      psipla = psizr(:)'*jphitot(:)/sum(sum(jphitot));

      % FLUX ERROR
      psizr_app = reshape(mpc*u.cc0t+mpv*u.vc0t,nz,nr);
      psizr_pla = gscalc(jphitot,rg,zg,mpp);
      psizr_err = psizr-psizr_pla-psizr_app;

      Ue = psizr_err(:);
      Ue(ii) = Ue(ii)*(-2/dr^2-2/dz^2)+Ue(ii-1)/dz^2+Ue(ii+1)/dz^2+...
	Ue(ii-nz).*(1/dr^2+1/2./rgg(ii)/dr)+Ue(ii+nz).*(1/dr^2-1/2./rgg(ii)/dr);

      % Various quantities
      cpasma = 1e6*dr*dz*sum(sum(jphitot));
      Vzr(ivac) = 0;
      Vzr(iplasma) = 2*pi*dr*dz*rgg(iplasma);
      Vzr(icell) = 2*pi*rgg(icell).*Ae(1:nxe);
      Vtot = sum(sum(Vzr));
      P(iplasma) = p0(iknot)+p1(iknot).*x+p2(iknot).*x.^2+p3(iknot).*x.^3;
      P(ivac) = 0;
      Wzr = dr*dz*3*pi*rgg.*P;
      Wzr(icell) = Ae*3*pi.*rgg(icell).*P(icell);
      Wth = sum(sum(Wzr));
      psizr_r = [zeros(nz,1) psizr(:,3:end)-psizr(:,1:end-2) zeros(nz,1)]/2/dr;
      psizr_z = [zeros(1,nr);psizr(3:end,:)-psizr(1:end-2,:);zeros(1,nr)]/2/dz;
      Bp2zr = (psizr_r.^2+psizr_z.^2)/twopi^2./rgg.^2;
      DX = diff(rbbbs)';
      DY = diff(zbbbs)';
      cpieces = sqrt(DX.^2+DY.^2);
      Cl = sum(cpieces); % Contour length
      bp2flx = (mu0*cpasma/Cl)^2;
      Bp2V = sum(sum(Vzr.*Bp2zr));      
      betap = 4/3*mu0*Wth/Vtot/bp2flx;
      li = Bp2V/Vtot/bp2flx; % Volume-averaged Bp^2 / bp2flx



      %                                                                                                                                                                
      %                                                                                                                                                              
      %             CCCCCCCCCCCCC                  lllllll                                       lllllll                           tttt                              
      %          CCC::::::::::::C                  l:::::l                                       l:::::l                        ttt:::t                              
      %        CC:::::::::::::::C                  l:::::l                                       l:::::l                        t:::::t                              
      %       C:::::CCCCCCCC::::C                  l:::::l                                       l:::::l                        t:::::t                              
      %      C:::::C       CCCCCC  aaaaaaaaaaaaa    l::::l     ccccccccccccccccuuuuuu    uuuuuu   l::::l   aaaaaaaaaaaaa  ttttttt:::::ttttttt        eeeeeeeeeeee    
      %     C:::::C                a::::::::::::a   l::::l   cc:::::::::::::::cu::::u    u::::u   l::::l   a::::::::::::a t:::::::::::::::::t      ee::::::::::::ee  
      %     C:::::C                aaaaaaaaa:::::a  l::::l  c:::::::::::::::::cu::::u    u::::u   l::::l   aaaaaaaaa:::::at:::::::::::::::::t     e::::::eeeee:::::ee
      %     C:::::C                         a::::a  l::::l c:::::::cccccc:::::cu::::u    u::::u   l::::l            a::::atttttt:::::::tttttt    e::::::e     e:::::e
      %     C:::::C                  aaaaaaa:::::a  l::::l c::::::c     cccccccu::::u    u::::u   l::::l     aaaaaaa:::::a      t:::::t          e:::::::eeeee::::::e
      %     C:::::C                aa::::::::::::a  l::::l c:::::c             u::::u    u::::u   l::::l   aa::::::::::::a      t:::::t          e:::::::::::::::::e 
      %     C:::::C               a::::aaaa::::::a  l::::l c:::::c             u::::u    u::::u   l::::l  a::::aaaa::::::a      t:::::t          e::::::eeeeeeeeeee  
      %      C:::::C       CCCCCCa::::a    a:::::a  l::::l c::::::c     cccccccu:::::uuuu:::::u   l::::l a::::a    a:::::a      t:::::t    tttttte:::::::e           
      %       C:::::CCCCCCCC::::Ca::::a    a:::::a l::::::lc:::::::cccccc:::::cu:::::::::::::::uul::::::la::::a    a:::::a      t::::::tttt:::::te::::::::e          
      %        CC:::::::::::::::Ca:::::aaaa::::::a l::::::l c:::::::::::::::::c u:::::::::::::::ul::::::la:::::aaaa::::::a      tt::::::::::::::t e::::::::eeeeeeee  
      %          CCC::::::::::::C a::::::::::aa:::al::::::l  cc:::::::::::::::c  uu::::::::uu:::ul::::::l a::::::::::aa:::a       tt:::::::::::tt  ee:::::::::::::e  
      %             CCCCCCCCCCCCC  aaaaaaaaaa  aaaallllllll    cccccccccccccccc    uuuuuuuu  uuuullllllll  aaaaaaaaaa  aaaa         ttttttttttt      eeeeeeeeeeeeee  
      %                                                                                                                                                              
      %                                                                                                                                                              
      %                                                                                                                                                              
      %                                                                                                                                                              
      %                                                                                                                                                              
      %                                                                                                                                                              
      %                                                                                                                                                              
      %                                                                                                                                                          
      %                                                                                                                                                          
      %                                                                                                                                                          
      %                                                                                                                                                          
      %                                                                                                                                                          
      %                                                                                                                                                          
      %     rrrrr   rrrrrrrrr       eeeeeeeeeeee        ssssssssss   ppppp   ppppppppp      ooooooooooo   nnnn  nnnnnnnn        ssssssssss       eeeeeeeeeeee    
      %     r::::rrr:::::::::r    ee::::::::::::ee    ss::::::::::s  p::::ppp:::::::::p   oo:::::::::::oo n:::nn::::::::nn    ss::::::::::s    ee::::::::::::ee  
      %     r:::::::::::::::::r  e::::::eeeee:::::eess:::::::::::::s p:::::::::::::::::p o:::::::::::::::on::::::::::::::nn ss:::::::::::::s  e::::::eeeee:::::ee
      %     rr::::::rrrrr::::::re::::::e     e:::::es::::::ssss:::::spp::::::ppppp::::::po:::::ooooo:::::onn:::::::::::::::ns::::::ssss:::::se::::::e     e:::::e
      %      r:::::r     r:::::re:::::::eeeee::::::e s:::::s  ssssss  p:::::p     p:::::po::::o     o::::o  n:::::nnnn:::::n s:::::s  ssssss e:::::::eeeee::::::e
      %      r:::::r     rrrrrrre:::::::::::::::::e    s::::::s       p:::::p     p:::::po::::o     o::::o  n::::n    n::::n   s::::::s      e:::::::::::::::::e 
      %      r:::::r            e::::::eeeeeeeeeee        s::::::s    p:::::p     p:::::po::::o     o::::o  n::::n    n::::n      s::::::s   e::::::eeeeeeeeeee  
      %      r:::::r            e:::::::e           ssssss   s:::::s  p:::::p    p::::::po::::o     o::::o  n::::n    n::::nssssss   s:::::s e:::::::e           
      %      r:::::r            e::::::::e          s:::::ssss::::::s p:::::ppppp:::::::po:::::ooooo:::::o  n::::n    n::::ns:::::ssss::::::se::::::::e          
      %      r:::::r             e::::::::eeeeeeee  s::::::::::::::s  p::::::::::::::::p o:::::::::::::::o  n::::n    n::::ns::::::::::::::s  e::::::::eeeeeeee  
      %      r:::::r              ee:::::::::::::e   s:::::::::::ss   p::::::::::::::pp   oo:::::::::::oo   n::::n    n::::n s:::::::::::ss    ee:::::::::::::e  
      %      rrrrrrr                eeeeeeeeeeeeee    sssssssssss     p::::::pppppppp       ooooooooooo     nnnnnn    nnnnnn  sssssssssss        eeeeeeeeeeeeee  
      %                                                               p:::::p                                                                                    
      %                                                               p:::::p                                                                                    
      %                                                              p:::::::p                                                                                   
      %                                                              p:::::::p                                                                                   
      %                                                              p:::::::p                                                                                   
      %                                                              ppppppppp                                                                                   
      %                                                                                                                                                          
      %         

      % Test if new plasma respnse matrix, Ti need be calculated
      if new_Ti

	ks = 1:nc+nv+6; % Means we will calculate all responses
	jphi0 = jphi;
	jphi_estimate = jphi;

	% Second derivatives
	Pbis = zeros(nz,nr); Pbis(iplasma) = (2*p2(iknot)+6*p3(iknot).*x)*twopi^2/(psibry-psimag)^2;
	Gbis = zeros(nz,nr); Gbis(iplasma) = (2*f2(iknot)+6*f3(iknot).*x)*twopi^2/(psibry-psimag)^2;
	djphidpsizr = (rgg.*Pbis+Gbis/mu0./rgg)/twopi/1e6;

	% Some variables to facilitate making equations
	x = psibarzr(iplasma);
	a = 1-x;
	R = rgg(iplasma);
	jpla = 1e6*jphi(iplasma)/(psimag-psibry);
	
	% dj1, djx, iplasma0 are uaed to calculate djdpsi
	dj1 = (2*R.*p2(iknot) + 2*f2(iknot)/mu0./R)*twopi;
	djx = (6*R.*p3(iknot) + 6*f3(iknot)/mu0./R)*twopi;
	iplasma0 = iplasma;
	djdpsi = (dj1+djx.*psibarzr(iplasma0))/(psibry-psimag)^2;
	djdpsi0 = djdpsi;
	sqsum_djdpsi0 = sum(djdpsi0.^2);
	iknot0 = iknot;
	jpla0 = jphi(iplasma0);
	jpla0_pred = jpla0;

	% Responses to spline coefficients that scale pprime, ffprim, or change peaking of ffprim
	jphipp = (p1_0(iknot)+p2_0(iknot)*2.*x+p3_0(iknot)*3.*x.^2)*twopi/(psibry-psimag).*R;    
	jphiff = (f1_0(iknot)+f2_0(iknot)*2.*x+f3_0(iknot)*3.*x.^2)*twopi/(psibry-psimag)./R/mu0;    
	jphigg = (g1_0(iknot)+g2_0(iknot)*2.*x+g3_0(iknot)*3.*x.^2)*twopi/(psibry-psimag)./R/mu0;    
	jss = [jphipp jphiff jphigg];

	% Us(ii) d(2*pi*mu0*j)/ds and Us(isides) = dpsizr/ds
	Us = zeros(ngg,3);
	Us(iplasma,:) = -2*pi*mu0*R*ones(1,3).*jss;
	Us(igridedge,:) = xpj(:,iplasma)*jss/1e6;

	% Add djphi to the equation delstar(dpsi)+mu0*R*djphi = 0
	T = T0+sparse(ii, ii, 2*pi*1e6*mu0*rgg(ii).*djphidpsizr(ii), ngg, ngg);

	T = full(T);

	T(igridedge,iplasma) = -ones(ni,1)*djphidpsizr(iplasma)'.*xpj(:,iplasma);

	T(iplasma,iia) = T(iplasma,iia)-2*pi*mu0*R.*(djdpsi.*a+jpla)*wa;
	T(igridedge,iia) = T(igridedge,iia) + xpj(:,iplasma)*(djdpsi.*a+jpla)*wa/1e6;

	T(iplasma,ii0) = T(iplasma,ii0)-2*pi*mu0*R.*(djdpsi.*x-jpla)*w0;
	T(igridedge,ii0) = T(igridedge,ii0) + xpj(:,iplasma)*(djdpsi.*x-jpla)*w0/1e6;

	T(icell,:) = T(icell,:) + 2*pi*1e6*mu0*rgg(icell)'*ones(1,ngg).*djphixedpsizr;
	T(igridedge,:) = T(igridedge,:) - xpj(:,icell)*djphixedpsizr;

	Ua(iplasma) = 2*pi*mu0*R.*djdpsi.*a;
	Ua(ivac) = 0;
	Ua(igridedge) = -xpj(:,iplasma)*(djdpsi.*a)/1e6;

	Ub(iplasma) = 2*pi*mu0*R.*djdpsi.*x;
	Ub(ivac) = 0;
	Ub(igridedge) = -xpj(:,iplasma)*(djdpsi.*x)/1e6;


	% Solve for the plasma response matrix
	Ti = inv(T);
	dpsizrdx = Ti*[mps Ua Ub Ue Us];

	% Reset total changes
	dpsimag_tot = 0;
	dpsibry_tot = 0;
	drbbbs_tot = 0;
	dzbbbs_tot = 0;

	% Remember quantities of the reference equilibrium to be used on response calculations below, the k-loop
	psimag0 = psimag;
	psibry0 = psibry;

	% This is the output for this equilibrium
	yprev = Cmat*[ic; iv; pcurrt(:)];

      else

	ks = nc+nv+3; % Means we will only calculate response to flux error
	dpsizrdx(:,ks) = Ti*Ue;

      end % End of condition for calculating new Ti


      a = 1-x;

      dpsizr = zeros(nz,nr);
      dPkn = zeros(nz,nr); % Pressure change on grid due to spline change of pres
      dPprimekn = zeros(nz,nr); % pprime change on grid due to spline change of pres
      dFFprimkn = zeros(nz,nr); % ffprim change on grid due to spline change of f^2/2

      % RESPONSES FOR TARGETED QUANTITIES
      if new_Ti
	drbbbsdx = zeros(nbbbs,nc+nv+6);
	dzbbbsdx = zeros(nbbbs,nc+nv+6);
      end
      drbbbs = zeros(nbbbs,1);
      dzbbbs = zeros(nbbbs,1);
      dre = zeros(nxe,1);
      dze = zeros(nxe,1);
      dP(ivac) = 0;
      djphixe = zeros(nz,nr);
      djphixekn = zeros(nxe,1);
      djphidx = zeros(ngg,nc+nv+6);
	x = psibarzr(iplasma);

      for k = ks

	dpsizr(:) = dpsizrdx(:,k);
	dpsimag = dpsizr(iia)*wa';
	dpsibry = dpsizr(ii0)*w0';

	if k == nc+nv+4
	  dPkn(iplasma) = p0(iknot)+p1(iknot).*x+p2(iknot).*x.^2+p3(iknot).*x.^3;
	  dPprimekn(iplasma) = (p1_0(iknot)+2*p2_0(iknot).*x+3*p3_0(iknot).*x.^2)*twopi/(psibry-psimag);
	  djphixekn = djphixeds(:,1);
	elseif k == nc+nv+5
	  dPkn(iplasma) = 0;
	  dPprimekn(iplasma) = 0;
	  dFFprimkn(iplasma) = (f1_0(iknot)+2*f2_0(iknot).*x+3*f3_0(iknot).*x.^2)*twopi/(psibry-psimag);
	  djphixekn = djphixeds(:,2);
	elseif k == nc+nv+6
	  dFFprimkn(iplasma) = (g1_0(iknot)+2*g2_0(iknot).*x+3*g3_0(iknot).*x.^2)*twopi/(psibry-psimag);
	  djphixekn = djphixeds(:,3);
	end

	djphi = djphidpsizr.*(dpsizr-dpsibry*psibarzr-dpsimag*(1-psibarzr))-...
        	jphi*(dpsimag-dpsibry)/(psimag-psibry)+(rgg.*dPprimekn+dFFprimkn/mu0./rgg)/1e6;

	djphixe(icell') = djphixedpsizr*dpsizr(:)+djphixekn;

	djphitot = djphi+djphixe;

	dpsipla = dpsizr(:)'*jphitot(:)/sum(sum(jphitot))+psizr(:)'*djphitot(:)/sum(sum(jphitot))-psipla*sum(sum(djphitot))/sum(sum(jphitot));

	dP(iplasma) = Pprime(iplasma).*(dpsizr(iplasma)-dpsibry*x-dpsimag*a)/twopi+dPkn(iplasma); % Note dP(ivac) can not be used

	dcpasma = 1e6*dr*dz*sum(sum(djphitot));
	dpsizr_r = [zeros(nz,1) dpsizr(:,3:end)-dpsizr(:,1:end-2) zeros(nz,1)]/2/dr;
	dpsizr_z = [zeros(1,nr);dpsizr(3:end,:)-dpsizr(1:end-2,:);zeros(1,nr)]/2/dz;
	ca = ashift_Ba*[dpsizr(iia)*war'; dpsizr(iia)*waz'];
	cx1 = x1shift_Bx*[dpsizr(iix1)*wx1r'; dpsizr(iix1)*wx1z'];
	cx2 = x2shift_Bx*[dpsizr(iix2)*wx2r'; dpsizr(iix2)*wx2z'];
	drx = [cx1(1); cx2(1)]; dzx = [cx1(2); cx2(2)];
	if ilimited
	  dlim = -dpsizr(iilim)*wld'/psilimbis*ulim;
	  dr0 = dlim(1); dz0 = dlim(2);
	elseif ix1
	  dr0 = cx1(1); dz0 = cx1(2);
	elseif ix2
	  dr0 = cx2(1); dz0 = cx2(2);
	end	  
	for j = 1:nbbbs
	  if flag_for_r_or_z_edge(j) == 1
            drbbbs(j) = -(wb(j,:)*dpsizr(ib(j)+neighbors)'-dpsibry)/dpsibbbsdr(j);
	  elseif flag_for_r_or_z_edge(j) == 2
            dzbbbs(j) = -(wb(j,:)*dpsizr(ib(j)+neighbors)'-dpsibry)/dpsibbbsdz(j);
	  else
            drbbbs(j) = dr0;
	    dzbbbs(j) = dz0;
	  end
	end
	for j = 1:nxe
          if abs(fe(j)) == 1
	    dre(j) = -(we(j,:)*dpsizr(ie(j)+neighbors)'-dpsibry)/dpsiedr(j);
	  else
	    dze(j) = -(we(j,:)*dpsizr(ie(j)+neighbors)'-dpsibry)/dpsiedz(j);
	  end
	end
	dAe = dAedre1'.*dre([nxe 1:nxe-1])+dAedze1'.*dze([nxe 1:nxe-1])+dAedre2'.*dre+dAedze2'.*dze;
	dDX = diff(drbbbs)';
	dDY = diff(dzbbbs)';
	dCl = sum((DX.*dDX+DY.*dDY)./cpieces); % Change of contour length
	% Change of betap and li
	dBp2V = 2*sum((psizr_r(iplasma).*dpsizr_r(iplasma)+psizr_z(iplasma).*dpsizr_z(iplasma))./rgg(iplasma))/twopi*dr*dz;
	% Store responses to x
	dpsibrydx(1,k) = dpsibry;
	dWth = 3*pi*dr*dz*rgg(iplasma)'*dP(iplasma); % pressure considered 0 in edge elements so ignore contribution from those
	dVtot = 2*pi*rgg(icell)*dAe;
	dWzr = dr*dz*3*pi*rgg.*dP;
	dWzr(icell) = dAe'*3*pi.*rgg(icell).*P(icell) + Ae*3*pi.*rgg(icell).*dP(icell);
	dWth = sum(sum(dWzr));
	dBp2zr = (2*psizr_r.*dpsizr_r+2*psizr_z.*dpsizr_z)/twopi^2./rgg.^2;
	dDX = diff(drbbbs)';
	dDY = diff(dzbbbs)';
	dcpieces = (DX.*dDX+DY.*dDY)./cpieces;
	dCl = sum(dcpieces); % Contour length
	dbp2flx = 2*bp2flx*(dcpasma/cpasma-dCl/Cl);
	dBp2V = sum(sum(Vzr.*dBp2zr))+2*pi*rgg(icell).*Bp2zr(icell)*dAe;
	dli = li*(dBp2V/Bp2V-dVtot/Vtot-dbp2flx/bp2flx);
	dbetap = betap*(dWth/Wth-dVtot/Vtot-dbp2flx/bp2flx);
	
	if new_Ti
          drbbbsdx(:,k) = drbbbs;
          dzbbbsdx(:,k) = dzbbbs;
          dpsimagdx(1,k) = dpsimag;
          dpsibrydx(1,k) = dpsibry;
	end

	dpcurrtdx(:,k) = 1e6*dr*dz*(djphi(:)+djphixe(:));
	dcpasmadx(1,k) = dcpasma;
	dlidx(1,k) = dli;
	dbetapdx(1,k) = dbetap;
	djphidx(:,k) = djphi(:);
	djphitotdx(:,k) = djphitot(:);
	ca = ashift_Ba*[dpsizr(iia)*war'; dpsizr(iia)*waz'];
	drmaxisdx(1,k) = ca(1);
	dzmaxisdx(1,k) = ca(2);
	dr0dx(1,k) = dr0;
	dz0dx(1,k) = dz0;


      end
      
      djpla0dx = djphidx(iplasma0,:);

      % How inputs respond to x-variables
      dudx = [dixdx; dcpasmadx; dlidx; dbetapdx];

      % How x-variables respond to input variables
      dxdu = inv(dudx);

      % How outputs respond to inputs
      dydx = Cmat*[disdx; dpcurrtdx];  

      da = 0; % psimag-psimag0-dpsimag_tot; % if psimag-psimag0 is different from dpsimag_tot due to axis displacement then correct for that
      db = 0; % psibry-psibry0-dpsibry_tot; % if psibry-psibry0 is different from dpsibry_tot due to r0,z0 displacement then correct for it
      de = -1; % Subtract flux error
      sum_djphi = zeros(nz,nr);

    else

      da = 0; % Save time by not finding this correction
      db = 0; % Save time by not finding this correction
      de = 0; % Do not subtract flux error

    end % End of condition for new equilibrium analysis





  %                                                                                                                           
  %     RRRRRRRRRRRRRRRRR                                tttt                                                                 
  %     R::::::::::::::::R                            ttt:::t                                                                 
  %     R::::::RRRRRR:::::R                           t:::::t                                                                 
  %     RR:::::R     R:::::R                          t:::::t                                                                 
  %       R::::R     R:::::R    eeeeeeeeeeee    ttttttt:::::ttttttt    uuuuuu    uuuuuu rrrrr   rrrrrrrrr   nnnn  nnnnnnnn    
  %       R::::R     R:::::R  ee::::::::::::ee  t:::::::::::::::::t    u::::u    u::::u r::::rrr:::::::::r  n:::nn::::::::nn  
  %       R::::RRRRRR:::::R  e::::::eeeee:::::eet:::::::::::::::::t    u::::u    u::::u r:::::::::::::::::r n::::::::::::::nn 
  %       R:::::::::::::RR  e::::::e     e:::::etttttt:::::::tttttt    u::::u    u::::u rr::::::rrrrr::::::rnn:::::::::::::::n
  %       R::::RRRRRR:::::R e:::::::eeeee::::::e      t:::::t          u::::u    u::::u  r:::::r     r:::::r  n:::::nnnn:::::n
  %       R::::R     R:::::Re:::::::::::::::::e       t:::::t          u::::u    u::::u  r:::::r     rrrrrrr  n::::n    n::::n
  %       R::::R     R:::::Re::::::eeeeeeeeeee        t:::::t          u::::u    u::::u  r:::::r              n::::n    n::::n
  %       R::::R     R:::::Re:::::::e                 t:::::t    ttttttu:::::uuuu:::::u  r:::::r              n::::n    n::::n
  %     RR:::::R     R:::::Re::::::::e                t::::::tttt:::::tu:::::::::::::::uur:::::r              n::::n    n::::n
  %     R::::::R     R:::::R e::::::::eeeeeeee        tt::::::::::::::t u:::::::::::::::ur:::::r              n::::n    n::::n
  %     R::::::R     R:::::R  ee:::::::::::::e          tt:::::::::::tt  uu::::::::uu:::ur:::::r              n::::n    n::::n
  %     RRRRRRRR     RRRRRRR    eeeeeeeeeeeeee            ttttttttttt      uuuuuuuu  uuuurrrrrrr              nnnnnn    nnnnnn
  %                                                                                                                           
  %                                                                                                                           
  %                                                                                                                           
  %                                                                                                                           
  %                                                                                                                           
  %                                                                                                                           
  %                                                                                                                           
  %                                                                                                                                        
  %                                                                                                                                        
  %                                               tttt                                                      tttt                           
  %                                            ttt:::t                                                   ttt:::t                           
  %                                            t:::::t                                                   t:::::t                           
  %                                            t:::::t                                                   t:::::t                           
  %        ooooooooooo   uuuuuu    uuuuuuttttttt:::::ttttttt   ppppp   ppppppppp   uuuuuu    uuuuuuttttttt:::::ttttttt        ssssssssss   
  %      oo:::::::::::oo u::::u    u::::ut:::::::::::::::::t   p::::ppp:::::::::p  u::::u    u::::ut:::::::::::::::::t      ss::::::::::s  
  %     o:::::::::::::::ou::::u    u::::ut:::::::::::::::::t   p:::::::::::::::::p u::::u    u::::ut:::::::::::::::::t    ss:::::::::::::s 
  %     o:::::ooooo:::::ou::::u    u::::utttttt:::::::tttttt   pp::::::ppppp::::::pu::::u    u::::utttttt:::::::tttttt    s::::::ssss:::::s
  %     o::::o     o::::ou::::u    u::::u      t:::::t          p:::::p     p:::::pu::::u    u::::u      t:::::t           s:::::s  ssssss 
  %     o::::o     o::::ou::::u    u::::u      t:::::t          p:::::p     p:::::pu::::u    u::::u      t:::::t             s::::::s      
  %     o::::o     o::::ou::::u    u::::u      t:::::t          p:::::p     p:::::pu::::u    u::::u      t:::::t                s::::::s   
  %     o::::o     o::::ou:::::uuuu:::::u      t:::::t    ttttttp:::::p    p::::::pu:::::uuuu:::::u      t:::::t    ttttttssssss   s:::::s 
  %     o:::::ooooo:::::ou:::::::::::::::uu    t::::::tttt:::::tp:::::ppppp:::::::pu:::::::::::::::uu    t::::::tttt:::::ts:::::ssss::::::s
  %     o:::::::::::::::o u:::::::::::::::u    tt::::::::::::::tp::::::::::::::::p  u:::::::::::::::u    tt::::::::::::::ts::::::::::::::s 
  %      oo:::::::::::oo   uu::::::::uu:::u      tt:::::::::::ttp::::::::::::::pp    uu::::::::uu:::u      tt:::::::::::tt s:::::::::::ss  
  %        ooooooooooo       uuuuuuuu  uuuu        ttttttttttt  p::::::pppppppp        uuuuuuuu  uuuu        ttttttttttt    sssssssssss    
  %                                                             p:::::p                                                                    
  %                                                             p:::::p                                                                    
  %                                                            p:::::::p                                                                   
  %                                                            p:::::::p                                                                   
  %                                                            p:::::::p                                                                   
  %                                                            ppppppppp                                                                   
  %                                                                                                                                        
  %     

    % Adjust plasma by linear response

    if ~isfield(u,'cc0t')
      u.cc0t = ic;
    end
    dic = u.cc0t-ic;

    if ~isfield(u,'vc0t')
      u.vc0t = iv;
    end
    div = u.vc0t-iv;

    if ~isfield(u,'ip')
      u.ip = cpasma;
    end
    dip = u.ip-cpasma;

    if ~isfield(u,'li')
      u.li = li;
    end
    dli = u.li-li;

    if ~isfield(u,'betap')
      u.betap = betap;
    end
    dbetap = u.betap-betap;
    
if 0
    dips(end+1) = dip;
    dlis(end+1) = dli;
    dbetaps(end+1) = dbetap;
    derrors(end+1) = sqrt(sum(psizr_err(:).^2));
    Cls(end+1) = Cl;
    Bp2Vs(end+1) = Bp2V;
    Vtots(end+1) = Vtot;
    bp2flxs(end+1) = bp2flx;
    r0s(end+1) = r0;
    Wths(end+1) = Wth;
end
    
if 0
    clf
    subplot(4,1,1)
    plot(dips)
    subplot(4,1,2)
    plot(dlis)
    subplot(4,1,3)
    plot(dbetaps)
    subplot(4,1,4)
    plot(psizr_err(:))
if 0
    clf
    plot(rl,zl,rbbbs,zbbbs)
    hold
    plot(r0,z0,'rx','markers',18,'linew',4)
    drawnow
end
    pause(1)
end
    du = [dic; div; da; db; de; dip; dli; dbetap];

    dx = dxdu*du;
    
    dpsizr = dpsizrdx*dx;
    djpla0 = djpla0dx*dx;    
    dsp    = sp0*dx(nc+nv+4);
    dsf    = sf0*dx(nc+nv+5)+sg0*dx(nc+nv+6);
    
    % Add the full desired change to psizr, jpla0_pred, sp, sf
    psizr(:)   = psizr(:)   + dpsizr;
    jpla0_pred = jpla0_pred + djpla0; % Linear prediction for jphi on grid indices iplasma0
    sp         = sp         + dsp;    % pres versus normalized flux
    sf         = sf         + dsf;    % fpol^2/2 versus normalized flux
    
    
    % Check if the adjustment is too large to be done in one step    
    psimag = wa*psizr(iia');
    psibry = w0*psizr(ii0');        
    p0 = c0*sp; p1 = c1*sp; p2 = c2*sp; p3 = c3*sp;
    f0 = c0*sf; f1 = c1*sf; f2 = c2*sf; f3 = c3*sf;
    
    x = (psizr(iplasma0)-psimag)/(psibry-psimag);
    
    % jpla0_est is an accurate estimate except that rmaxis, zmaxis & r0, z0 haven't been analyzed
    jpla0_est = ((p1(iknot0)+p2(iknot0)*2.*x+p3(iknot0)*3.*x.^2)*twopi/(psibry-psimag).*rgg(iplasma0) + ...
                 (f1(iknot0)+f2(iknot0)*2.*x+f3(iknot0)*3.*x.^2)*twopi/(psibry-psimag)./rgg(iplasma0)/mu0)/1e6;
		 
    % This is a good approximation of the error in the linear prediction of jpla0
    jpla0_err = jpla0_pred - jpla0_est;
    
    rel_nonlin_err = sqrt(sum(jpla0_err.^2)/(1+sum(djpla0.^2)));
    
    damper = 10*rel_nonlin_err;
    
    if damper > 1 & iteration_counter < 9
      f = 1-1/damper; % In this case we need to subtract off a fraction, f of the changes already made
      psizr(:)   = psizr(:)   - f*dpsizr;
      jpla0_pred = jpla0_pred - f*djpla0;     
      sp         = sp         - f*dsp;
      sf         = sf         - f*dsf;
      dx = dx/damper;
      not_done = true;
    else
      not_done = false;
    end
  
    % Update remaining plasma parameters
    ic      = ic      + dx(1:nc);
    iv      = iv      + dx(nc+(1:nv));    
    cpasma  = cpasma  + dcpasmadx*dx;
    li      = li      + dlidx*dx;
    betap   = betap   + dbetapdx*dx;

    % Finally update the outputs
    y = yprev+dydx*dx;
    yprev = y;
    
    iteration_counter = iteration_counter+1;
    
  end % End of while loop with condition not_done 
  
  if nargout > 1
    % Return persistent variables in e
    clear e
    if 01 % Do this when set of persistent variables is under development
      w = which('evolveq');
      f = fopen(w);
      s = char(fread(f))';
      fclose(f);
      k = strfind(s, [10 32 32 'persistent']);
      pvars = '';
      for j = 1:length(k)
	i1 = k(j)+3+length('persistent');
	i2 = strfind(s(i1:end),10);
	pvars = eval(['strvcat(pvars,''' strrep(s(i1+1:i1+i2-2),' ',''',''') ''');']);
      end
      for j = 1:size(pvars,1)
	eval(['e.' pvars(j,:) ' = ' pvars(j,:) ';'])
      end
      if 0 % This will make a file with the matlab code that initializes evolveq4simulink
        f = fopen('load_persistent_from_e.m','w');
        for j = 1:size(pvars,1)
	  fprintf(f,'%s\n',['    ' pvars(j,:) ' = e.' pvars(j,:) ';']);
	end
	fclose(f);
      end
      if 0 % This will make a file with the matlab code that goes after the else below
        f = fopen('store_persistent_to_e.m','w');
        for j = 1:size(pvars,1)
	  fprintf(f,'%s\n',['      e.' pvars(j,:) ' = ' pvars(j,:) ';']);
	end
	fclose(f);
      end
    else % Paste store_persistent_to_e.m here
      e.Cmat             = Cmat            ;
      e.nr               = nr              ;
      e.nz               = nz              ;
      e.ngg              = ngg             ;
      e.rg               = rg              ;
      e.zg               = zg              ;
      e.dr               = dr              ;
      e.dz               = dz              ;
      e.rgg              = rgg             ;
      e.zgg              = zgg             ;
      e.mpc              = mpc             ;
      e.mpv              = mpv             ;
      e.mps              = mps             ;
      e.mpp              = mpp             ;
      e.xpj              = xpj             ;
      e.Rlim             = Rlim            ;
      e.Zlim             = Zlim            ;
      e.nlim             = nlim            ;
      e.rl               = rl              ;
      e.zl               = zl              ;
      e.nl               = nl              ;
      e.wls              = wls             ;
      e.iil              = iil             ;
      e.ilimgg           = ilimgg          ;
      e.klimgg           = klimgg          ;
      e.piccc            = piccc           ;
      e.nc               = nc              ;
      e.nv               = nv              ;
      e.mu0              = mu0             ;
      e.twopi            = twopi           ;
      e.R13              = R13             ;
      e.mx               = mx              ;
      e.neighbors        = neighbors       ;
      e.nkn              = nkn             ;
      e.psikn            = psikn           ;
      e.T0               = T0              ;
      e.T                = T               ;
      e.sp0              = sp0             ;
      e.sf0              = sf0             ;
      e.sg0              = sg0             ;
      e.c0               = c0              ;
      e.c1               = c1              ;
      e.c2               = c2              ;
      e.c3               = c3              ;
      e.iknotg           = iknotg          ;
      e.psibar           = psibar          ;
      e.v1r              = v1r             ;
      e.v2r              = v2r             ;
      e.v3r              = v3r             ;
      e.v4r              = v4r             ;
      e.v1z              = v1z             ;
      e.v2z              = v2z             ;
      e.v3z              = v3z             ;
      e.v4z              = v4z             ;
      e.ii               = ii              ;
      e.ni               = ni              ;
      e.igridedge        = igridedge       ;
      e.dP               = dP              ;
      e.dixdx            = dixdx           ;
      e.disdx            = disdx           ;
      e.Ua               = Ua              ;
      e.Ub               = Ub              ;
      e.p1_0             = p1_0            ;
      e.p2_0             = p2_0            ;
      e.p3_0             = p3_0            ;
      e.p4_0             = p4_0            ;
      e.f1_0             = f1_0            ;
      e.f2_0             = f2_0            ;
      e.f3_0             = f3_0            ;
      e.f4_0             = f4_0            ;
      e.g1_0             = g1_0            ;
      e.g2_0             = g2_0            ;
      e.g3_0             = g3_0            ;
      e.g4_0             = g4_0            ;
      e.psimag0          = psimag0         ;
      e.psibry0          = psibry0         ;
      e.jpla0            = jpla0           ;
      e.iplasma0         = iplasma0        ;
      e.iknot0           = iknot0          ;
      e.djdpsi0          = djdpsi0         ;
      e.sqsum_djdpsi0    = sqsum_djdpsi0   ;
      e.jpla0_pred       = jpla0_pred      ;
      e.jpla0_est        = jpla0_est       ;
      e.cpasma           = cpasma          ;
      e.rmaxis           = rmaxis          ;
      e.zmaxis           = zmaxis          ;
      e.psimag           = psimag          ;
      e.r0               = r0              ;
      e.z0               = z0              ;
      e.psibry           = psibry          ;
      e.x1exists         = x1exists        ;
      e.x1inside         = x1inside        ;
      e.rx1              = rx1             ;
      e.zx1              = zx1             ;
      e.x2exists         = x2exists        ;
      e.x2inside         = x2inside        ;
      e.rx2              = rx2             ;
      e.zx2              = zx2             ;
      e.rzero            = rzero           ;
      e.bzero            = bzero           ;
      e.Pprime           = Pprime          ;
      e.FFprim           = FFprim          ;
      e.P                = P               ;
      e.iplasma          = iplasma         ;
      e.rbbbs            = rbbbs           ;
      e.zbbbs            = zbbbs           ;
      e.zbot             = zbot            ;
      e.ztop             = ztop            ;
      e.psipla           = psipla          ;
      e.li               = li              ;
      e.betap            = betap           ;
      e.jphi             = jphi            ;
      e.pcurrt           = pcurrt          ;
      e.ii0              = ii0             ;
      e.w0               = w0              ;
      e.iia              = iia             ;
      e.wa               = wa              ;
      e.ilimited         = ilimited        ;
      e.ix1              = ix1             ;
      e.ix2              = ix2             ;
      e.nbbbs            = nbbbs           ;
      e.fpol             = fpol            ;
      e.ffprim           = ffprim          ;
      e.pres             = pres            ;
      e.pprime           = pprime          ;
      e.psix1            = psix1           ;
      e.psix2            = psix2           ;
      e.thbbbs           = thbbbs          ;
      e.Wth              = Wth             ;
      e.psizr_err        = psizr_err       ;
      e.Vzr              = Vzr             ;
      e.Bp2zr            = Bp2zr           ;
      e.Cl               = Cl              ;
      e.Bp2V             = Bp2V            ;
      e.Vtot             = Vtot            ;
      e.bp2flx           = bp2flx          ;
      e.ic               = ic              ;
      e.iv               = iv              ;
      e.psizr            = psizr           ;
      e.sp               = sp              ;
      e.sf               = sf              ;
      e.yprev            = yprev           ;
      e.psibarzr         = psibarzr        ;
      e.Ti               = Ti              ;
      e.dpsizrdx         = dpsizrdx        ;
      e.dpsimagdx        = dpsimagdx       ;
      e.dpsibrydx        = dpsibrydx       ;
      e.djphidpsizr      = djphidpsizr     ;
      e.dpcurrtdx        = dpcurrtdx       ;
      e.dcpasmadx        = dcpasmadx       ;
      e.dlidx            = dlidx           ;
      e.dbetapdx         = dbetapdx        ;
      e.dydx             = dydx            ;
      e.dxdu             = dxdu            ;
      e.drbbbsdx         = drbbbsdx        ;
      e.dzbbbsdx         = dzbbbsdx        ;
      e.djphidx          = djphidx         ;
      e.djpla0dx         = djpla0dx        ;
      e.dj1              = dj1             ;
      e.djx              = djx             ;
      e.old_value        = old_value       ;
      e.new_value        = new_value       ;
      e.predicted_change = predicted_change;
      e.dips             = dips            ;
      e.dlis             = dlis            ;
      e.dbetaps          = dbetaps         ;
      e.derrors          = derrors         ;
      e.Cls              = Cls             ;
      e.Bp2Vs            = Bp2Vs           ;
      e.Vtots            = Vtots           ;
      e.bp2flxs          = bp2flxs         ;
      e.r0s              = r0s             ;
      e.Wths             = Wths            ;
    end
  end
  
  
  
  
  
  
  
%                                                                                                                           
%                                                                                                                           
%     RRRRRRRRRRRRRRRRR                                tttt                                                                 
%     R::::::::::::::::R                            ttt:::t                                                                 
%     R::::::RRRRRR:::::R                           t:::::t                                                                 
%     RR:::::R     R:::::R                          t:::::t                                                                 
%       R::::R     R:::::R    eeeeeeeeeeee    ttttttt:::::ttttttt    uuuuuu    uuuuuu rrrrr   rrrrrrrrr   nnnn  nnnnnnnn    
%       R::::R     R:::::R  ee::::::::::::ee  t:::::::::::::::::t    u::::u    u::::u r::::rrr:::::::::r  n:::nn::::::::nn  
%       R::::RRRRRR:::::R  e::::::eeeee:::::eet:::::::::::::::::t    u::::u    u::::u r:::::::::::::::::r n::::::::::::::nn 
%       R:::::::::::::RR  e::::::e     e:::::etttttt:::::::tttttt    u::::u    u::::u rr::::::rrrrr::::::rnn:::::::::::::::n
%       R::::RRRRRR:::::R e:::::::eeeee::::::e      t:::::t          u::::u    u::::u  r:::::r     r:::::r  n:::::nnnn:::::n
%       R::::R     R:::::Re:::::::::::::::::e       t:::::t          u::::u    u::::u  r:::::r     rrrrrrr  n::::n    n::::n
%       R::::R     R:::::Re::::::eeeeeeeeeee        t:::::t          u::::u    u::::u  r:::::r              n::::n    n::::n
%       R::::R     R:::::Re:::::::e                 t:::::t    ttttttu:::::uuuu:::::u  r:::::r              n::::n    n::::n
%     RR:::::R     R:::::Re::::::::e                t::::::tttt:::::tu:::::::::::::::uur:::::r              n::::n    n::::n
%     R::::::R     R:::::R e::::::::eeeeeeee        tt::::::::::::::t u:::::::::::::::ur:::::r              n::::n    n::::n
%     R::::::R     R:::::R  ee:::::::::::::e          tt:::::::::::tt  uu::::::::uu:::ur:::::r              n::::n    n::::n
%     RRRRRRRR     RRRRRRR    eeeeeeeeeeeeee            ttttttttttt      uuuuuuuu  uuuurrrrrrr              nnnnnn    nnnnnn
%                                                                                                                           
%                                                                                                                           
%                                                                                                                           
%                                                                                                                           
%                                                                                                                           
%                                                                                                                           
%                                                                                                                           
%                                                                                                                                                                               
%                                                                                      bbbbbbbb                                                                                 
%                                                                 iiii  lllllll   iiii b::::::b                                 iiii                                            
%                                                                i::::i l:::::l  i::::ib::::::b                                i::::i                                           
%                                                                 iiii  l:::::l   iiii b::::::b                                 iiii                                            
%                                                                       l:::::l         b:::::b                                                                                 
%         eeeeeeeeeeee       qqqqqqqqq   qqqqquuuuuu    uuuuuu  iiiiiii  l::::l iiiiiii b:::::bbbbbbbbb    rrrrr   rrrrrrrrr  iiiiiii uuuuuu    uuuuuu     mmmmmmm    mmmmmmm   
%       ee::::::::::::ee    q:::::::::qqq::::qu::::u    u::::u  i:::::i  l::::l i:::::i b::::::::::::::bb  r::::rrr:::::::::r i:::::i u::::u    u::::u   mm:::::::m  m:::::::mm 
%      e::::::eeeee:::::ee q:::::::::::::::::qu::::u    u::::u   i::::i  l::::l  i::::i b::::::::::::::::b r:::::::::::::::::r i::::i u::::u    u::::u  m::::::::::mm::::::::::m
%     e::::::e     e:::::eq::::::qqqqq::::::qqu::::u    u::::u   i::::i  l::::l  i::::i b:::::bbbbb:::::::brr::::::rrrrr::::::ri::::i u::::u    u::::u  m::::::::::::::::::::::m
%     e:::::::eeeee::::::eq:::::q     q:::::q u::::u    u::::u   i::::i  l::::l  i::::i b:::::b    b::::::b r:::::r     r:::::ri::::i u::::u    u::::u  m:::::mmm::::::mmm:::::m
%     e:::::::::::::::::e q:::::q     q:::::q u::::u    u::::u   i::::i  l::::l  i::::i b:::::b     b:::::b r:::::r     rrrrrrri::::i u::::u    u::::u  m::::m   m::::m   m::::m
%     e::::::eeeeeeeeeee  q:::::q     q:::::q u::::u    u::::u   i::::i  l::::l  i::::i b:::::b     b:::::b r:::::r            i::::i u::::u    u::::u  m::::m   m::::m   m::::m
%     e:::::::e           q::::::q    q:::::q u:::::uuuu:::::u   i::::i  l::::l  i::::i b:::::b     b:::::b r:::::r            i::::i u:::::uuuu:::::u  m::::m   m::::m   m::::m
%     e::::::::e          q:::::::qqqqq:::::q u:::::::::::::::uui::::::il::::::li::::::ib:::::bbbbbb::::::b r:::::r           i::::::iu:::::::::::::::uum::::m   m::::m   m::::m
%      e::::::::eeeeeeee   q::::::::::::::::q  u:::::::::::::::ui::::::il::::::li::::::ib::::::::::::::::b  r:::::r           i::::::i u:::::::::::::::um::::m   m::::m   m::::m
%       ee:::::::::::::e    qq::::::::::::::q   uu::::::::uu:::ui::::::il::::::li::::::ib:::::::::::::::b   r:::::r           i::::::i  uu::::::::uu:::um::::m   m::::m   m::::m
%         eeeeeeeeeeeeee      qqqqqqqq::::::q     uuuuuuuu  uuuuiiiiiiiilllllllliiiiiiiibbbbbbbbbbbbbbbb    rrrrrrr           iiiiiiii    uuuuuuuu  uuuummmmmm   mmmmmm   mmmmmm
%                                     q:::::q                                                                                                                                   
%                                     q:::::q                                                                                                                                   
%                                    q:::::::q                                                                                                                                  
%                                    q:::::::q                                                                                                                                  
%                                    q:::::::q                                                                                                                                  
%                                    qqqqqqqqq                                                                                                                                  
%                                                                                                                                                                               
%      
  
  
  
  
  % output eq, all we really know at this point are: ic iv psizr sp sf 
  
  if nargout > 2 % Create structure with present equilibrium
    
    
    % Analyze the equilibrium
    
    % Remember last analyzed equilibrium
    eqx.cpasma    = cpasma;
    eqx.rmaxis    = rmaxis;
    eqx.zmaxis    = zmaxis;
    eqx.psimag    = psimag;
    eqx.r0        = r0;
    eqx.z0        = z0;
    eqx.psibry    = psibry;
    eqx.x1exists  = x1exists;
    eqx.x1inside  = x1inside;
    eqx.rx1       = rx1;
    eqx.zx1       = zx1;
    eqx.x2exists  = x2exists;
    eqx.x2inside  = x2inside;
    eqx.rx2       = rx2;
    eqx.zx2       = zx2;
    eqx.rzero     = rzero;
    eqx.bzero     = bzero;
    eqx.Pprime    = Pprime;
    eqx.FFprim    = FFprim;
    eqx.P         = P;
    eqx.iplasma   = iplasma;
    eqx.rbbbs     = rbbbs;
    eqx.zbbbs     = zbbbs;
    eqx.zbot      = zbot;
    eqx.ztop      = ztop;
    eqx.psipla    = psipla;
    eqx.li        = li;
    eqx.betap     = betap;
    eqx.jphi      = jphi;
    eqx.pcurrt    = pcurrt;
    eqx.ii0       = ii0;
    eqx.w0        = w0;
    eqx.iia       = iia;
    eqx.wa        = wa;
    eqx.ilimited  = ilimited;
    eqx.ix1       = ix1;
    eqx.ix2       = ix2;
    eqx.nbbbs     = nbbbs;
    eqx.fpol      = fpol;
    eqx.ffprim    = ffprim;
    eqx.pres      = pres;
    eqx.pprime    = pprime;
    eqx.psix1     = psix1;
    eqx.psix2     = psix2;
    eqx.thbbbs    = thbbbs;
    eqx.Wth       = Wth;
    eqx.psizr_err = psizr_err;
    eqx.Vzr       = Vzr;
    eqx.Bp2zr     = Bp2zr;
    eqx.Cl        = Cl;
    eqx.Bp2V      = Bp2V;
    eqx.Vtot      = Vtot;
    eqx.bp2flx    = bp2flx;
    eqx.psibarzr  = psibarzr;
    
    % Store some more diagnostics of the simulation
    eqx.rel_lin_err       = rel_lin_err;
    eqx.rel_nonlin_err    = rel_nonlin_err;
    eqx.rel_totlin_err    = rel_totlin_err;
    eqx.new_Ti            = new_Ti;
    eqx.new_eq_analysis   = new_eq_analysis;
    eqx.rel_nonlin_err    = rel_nonlin_err;
    eqx.iteration_counter = iteration_counter;
    
    % Magnetic axis
    twopirbrzmax = inf; % Maximum of 2*pi*R*Br, 2*pi*R*Bz at calculated null position
    j = 9;
    while j>0 & twopirbrzmax>1e-10
      kr0 = min(nr-3,max(1,floor((rmaxis-rg(1))/dr)));
      kz1 = min(nz-2,max(2,ceil((zmaxis-zg(1))/dz)));
      k = kr0*nz+kz1;
      iia = k+neighbors; % iia indexes 16 points around magnetic axis
      pp = psizr(iia);
      tr = (rmaxis-rgg(k))/dr;
      tz = (zmaxis-zgg(k))/dz;
      wa = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      war = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
      waz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      warr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]/dr^2*mx,1,16);
      wazz = reshape(([0 0 2 6*tz]/dz^2*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      warz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
      ashift_Ba = -inv([pp*warr' pp*warz';pp*warz' pp*wazz']);
      c = ashift_Ba*[pp*war'; pp*waz'];
      rmaxis = rmaxis+c(1);
      zmaxis = zmaxis+c(2);
      if norm(c) > dr
	[dum,k] = min(psizr(iplasma)/(psibry-psimag));
	rmaxis = rgg(iplasma(k));
	zmaxis = zgg(iplasma(k));
      end
      j = j-1;
      twopirbrzmax = max(abs([pp*war' pp*waz']));
    end
    psimag = wa*pp';

    psizr_r = [zeros(nz,1) psizr(:,3:end)-psizr(:,1:end-2) zeros(nz,1)]/2/dr;
    psizr_z = [zeros(1,nr);psizr(3:end,:)-psizr(1:end-2,:);zeros(1,nr)]/2/dz;

    % Check for LOWER NULL
    if ~(x1exists & x1inside)
      j = find(zbbbs < zmaxis);
      [dum, k] = min(interp2(rg,zg,psizr_r.^2+psizr_z.^2,rbbbs(j),zbbbs(j),'spline'));
      rx1 = rbbbs(j(k));
      zx1 = zbbbs(j(k));
    end

    twopirbrzmax = inf; % Maximum of 2*pi*R*Br, 2*pi*R*Bz at calculated null position
    j = 9;
    while j > 0 & twopirbrzmax > 1e-10 % Try zooming in on x-point with Newton Rhapson.
      j = j-1;
      % Find indices and weights for grid points around the x-point
      kr0 = min(nr-3,max(1,floor((rx1-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
      kz1 = min(nz-2,max(2,ceil((zx1-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
      k = kr0*nz+kz1;
      iix1 = k+neighbors; % indexes 16 points around x point
      pp = psizr(iix1)';
      tr = (rx1-rgg(k))/dr;
      tz = (zx1-zgg(k))/dz;
      wx1 = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      wx1r = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
      wx1z = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      wx1rr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]/dr^2*mx,1,16);
      wx1zz = reshape(([0 0 2 6*tz]/dz^2*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      wx1rz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
      x1shift_Bx = -inv([wx1rr*pp wx1rz*pp;wx1rz*pp wx1zz*pp]);
      cx1 = x1shift_Bx*[wx1r*pp; wx1z*pp];
      rx1 = rx1+cx1(1);
      zx1 = zx1+cx1(2);
      twopirbrzmax = max(abs([wx1r*pp wx1z*pp]));
    end
    psix1 = wx1*pp;

    if twopirbrzmax > 1e-10
      x1exists = false; % No x-point found
    else
      x1exists = true; % x-point found
    end
    if x1exists % X-point found. Now check if it is inside the limiter
      if ilimgg(k) == 0 % In this case it is clearly inside
	x1inside = true; % Flag that x1 is a real (and important) x-point	  
      elseif ilimgg(k) == -1 % That means clearly outside
	x1inside = false;
      else
	x1inside = true;
	for j = max(1,ilimgg(k)-5):min(nlim-1,ilimgg(k)+5);
          mlinex = [Rlim(j+1)-Rlim(j) rmaxis-rx1; Zlim(j+1)-Zlim(j) zmaxis-zx1];
	  if rcond(mlinex) > 0
	    kg = inv(mlinex)*[rmaxis-Rlim(j); zmaxis-Zlim(j)];
	    if kg(1)>=0 & kg(1)<=1 & kg(2)>0 & kg(2)<1 % x1 outside limiter
	      x1inside = false;
	    end
	  end
	end
      end
    end

    if x1exists & x1inside
      zbot = zx1; % Min z for where plasma can be
    else
      psix1 = inf/(psibry-psimag);
      zbot = min(zbbbs)-dz;
    end
    psibarx1 = (psix1-psimag)/(psibry-psimag);


    % Check for UPPER NULL
    if ~(x2exists & x2inside)
      j = find(zbbbs > zmaxis);
      [dum, k] = min(interp2(rg,zg,psizr_r.^2+psizr_z.^2,rbbbs(j),zbbbs(j),'spline'));
      rx2 = rbbbs(j(k));
      zx2 = zbbbs(j(k));
    end

    twopirbrzmax = inf; % Maximum of 2*pi*R*Br, 2*pi*R*Bz at calculated null position
    j = 9;
    while j > 0 & twopirbrzmax > 1e-10 % Try zooming in on x-point with Newton Rhapson.
      j = j-1;
      % Find indices and weights for grid points around the x-point
      kr0 = min(nr-3,max(1,floor((rx2-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
      kz1 = min(nz-2,max(2,ceil((zx2-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
      k = kr0*nz+kz1;
      iix2 = k+neighbors; % indexes 16 points around x point
      pp = psizr(iix2)';
      tr = (rx2-rgg(k))/dr;
      tz = (zx2-zgg(k))/dz;
      wx2 = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      wx2r = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
      wx2z = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      wx2rr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]/dr^2*mx,1,16);
      wx2zz = reshape(([0 0 2 6*tz]/dz^2*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      wx2rz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
      x2shift_Bx = -inv([wx2rr*pp wx2rz*pp;wx2rz*pp wx2zz*pp]);
      cx2 = x2shift_Bx*[wx2r*pp; wx2z*pp];
      rx2 = rx2+cx2(1);
      zx2 = zx2+cx2(2);
      twopirbrzmax = max(abs([wx2r*pp wx2z*pp]));
    end
    psix2 = wx2*pp;

    if twopirbrzmax > 1e-10
      x2exists = false; % No x-point found
    else
      x2exists = true; % x-point found
    end
    if x2exists % X-point found. Now check if it is inside the limiter
      if ilimgg(k) == 0 % In this case it is clearly inside
	x2inside = true; % Flag that x2 is a real (and important) x-point	  
      elseif ilimgg(k) == -1 % That means clearly outside
	x2inside = false;
      else
	x2inside = true;
	for j = max(1,ilimgg(k)-5):min(nlim-1,ilimgg(k)+5);
          mlinex = [Rlim(j+1)-Rlim(j) rmaxis-rx2; Zlim(j+1)-Zlim(j) zmaxis-zx2];
	  if rcond(mlinex) > 0
	    kg = inv(mlinex)*[rmaxis-Rlim(j); zmaxis-Zlim(j)];
	    if kg(1)>=0 & kg(1)<=1 & kg(2)>0 & kg(2)<1 % x2 outside limiter
	      x2inside = false;
	    end
	  end
	end
      end
    end

    if x2exists & x2inside
      ztop = zx2; % Max z for where plasma can be
    else
      psix2 = inf/(psibry-psimag);
      ztop = max(zbbbs)+dz;
    end
    psibarx2 = (psix2-psimag)/(psibry-psimag);


    % Find TOUCH POINT candidate
    for j = 1:nl % Calculate for all of limiter to find touch and strike points
      psilim(j) = wls(j,:)*psizr(iil(j,:))';
    end

    k = find(rl>rg(1) & rl<rg(end) & zl>zg(1) & zl<zg(end));
    if x1exists & x1inside
      k = intersect(k,find(zl > zx1));
    end
    if x2exists & x2inside
      k = intersect(k,find(zl < zx2));
    end
    [dum, l] = min(psilim(k)/(psibry-psimag)); k1 = k(l);
    rlim = rl(k1); zlim = zl(k1);
    k2 = k1+1; if k2 > length(rl), k2 = 2; end
    ulim = [rl(k2)-rl(k1) , zl(k1)-zl(k2)]; nulim = norm(ulim);
    ulim = ulim/nulim;
    kr0 = min(nr-3,max(1,floor((rlim-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
    kz1 = min(nz-2,max(2,ceil((zlim-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
    k = kr0*nz+kz1; iilim = k+neighbors;
    tr = (rlim-rgg(k))/dr; tz = (zlim-zgg(k))/dz;
    wl = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
    wld = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]*mx*ulim(1)/dr+...
                   ([0 1 2*tz 3*tz^2]*mx)'*ulim(2)/dz*[1 tr tr^2 tr^3]*mx,1,16);
    wlb = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]*mx*ulim(1)^2/dr^2+...
        	 ([0 0 2 6*tz]*mx)'*ulim(2)^2/dz^2*[1 tr tr^2 tr^3]*mx+...
	       2*([0 1 2*tz 3*tz^2]*mx)'*ulim(2)/dz*[0 1 2*tr 3*tr^2]*mx*ulim(1)/dr,1,16);
    psilimbis = psizr(iilim)*wlb';
    dlim = -psizr(iilim)*wld'/psilimbis*ulim;
    for j = 1:3
      if norm(dlim)<nulim
	rlim = rlim+dlim(1); zlim = zlim+dlim(2);
	kr0 = min(nr-3,max(1,floor((rlim-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
	kz1 = min(nz-2,max(2,ceil((zlim-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
	k = kr0*nz+kz1; iilim = k+neighbors;
	tr = (rlim-rgg(k))/dr; tz = (zlim-zgg(k))/dz;
	wl = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
	wld = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]*mx*ulim(1)/dr+...
                       ([0 1 2*tz 3*tz^2]*mx)'*ulim(2)/dz*[1 tr tr^2 tr^3]*mx,1,16);
	wlb = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]*mx*ulim(1)^2/dr^2+...
                     ([0 0 2 6*tz]*mx)'*ulim(2)^2/dz^2*[1 tr tr^2 tr^3]*mx+...
		   2*([0 1 2*tz 3*tz^2]*mx)'*ulim(2)/dz*[0 1 2*tr 3*tr^2]*mx*ulim(1)/dr,1,16);
	psilimbis = psizr(iilim)*wlb';
	dlim = -psizr(iilim)*wld'/psilimbis*ulim;
      end
    end
    psitouch = psizr(iilim)*wl';
    psibartouch = (psitouch-psimag)/(psibry-psimag);

    % Which point defines the boundary?
    ilimited = 0; ix1 = 0; ix2 = 0; % Flags for what defines boundary
    if psibartouch < min(psibarx1,psibarx2)
      ilimited = 1;
      r0 = rlim;
      z0 = zlim;
      w0 = wl;
      ii0 = iilim;
      psibry = psitouch;
    else
      if psibarx1 < psibarx2
	ix1 = 1;
	r0 = rx1;
	z0 = zx1;
	w0 = wx1;
	ii0 = iix1;
	psibry = psix1;
      else
	ix2 = 1;
	r0 = rx2;
	z0 = zx2;
	w0 = wx2;
	ii0 = iix2;
	psibry = psix2;
      end
    end
    psibarzr = (psizr-psimag)/(psibry-psimag);

    % Trace boundary
    % When tz = 0: psi = [1 tr tr^2 tr^3]*mx*psizr(k+[-nz 0 nz 2*nz]')
    k = find((psizr(klimgg)-psibry).*(psizr(klimgg+nz)-psibry) < 0);
    k = klimgg(k);
    d = psizr(k)-psibry;
    c = (psizr(k+nz)-psizr(k-nz))/2;
    b = (2*psizr(k-nz)-5*psizr(k)+4*psizr(k+nz)-psizr(k+2*nz))/2;
    a = (-psizr(k-nz)+3*psizr(k)-3*psizr(k+nz)+psizr(k+2*nz))/2;
    aa = 2*b.^3-9*a.*b.*c+27*a.^2.*d;
    bb = sqrt(aa.^2-4*(b.^2-3*a.*c).^3);
    xx = d./(psizr(k)-psizr(k+nz));

    q1 = ((aa+bb)/2).^(1/3);
    q2 = ((aa-bb)/2).^(1/3);

    j = find(imag(q1+q2)>1e-6);	q2(j) = R13*q2(j);
    j = find(imag(q1+q2)>1e-6);	q2(j) = R13*q2(j);

    x1 = -(b+real(q1+q2))/3./a;
    x2 = -(b-(1+i*sqrt(3))/2*q1-(1-i*sqrt(3))/2*q2)/3./a;
    x3 = -(b-(1-i*sqrt(3))/2*q1-(1+i*sqrt(3))/2*q2)/3./a;

    for j = 1:length(xx)
      if     abs(imag(x1(j))) < 1e-6 & x1(j) > 0 & x1(j) < 1
	xx(j) = abs(x1(j));
      elseif abs(imag(x2(j))) < 1e-6 & x2(j) > 0 & x2(j) < 1
	xx(j) = abs(x2(j));
      elseif abs(imag(x3(j))) < 1e-6 & x3(j) > 0 & x3(j) < 1
	xx(j) = abs(x3(j));
      end
    end

    rbbbs = rgg(k)+xx*dr;
    zbbbs = zgg(k);
    ib = k;
    flag_for_r_or_z_edge = k*0+1; % 1 for tz=0 and 2 for tr=0

    % When tr = 0: psi = [1 tz tz^2 tz^3]*mx*psizr(k+[-1 0 1 2]')
    k = find((psizr(klimgg)-psibry).*(psizr(klimgg+1)-psibry) < 0);
    k = klimgg(k);
    d = psizr(k)-psibry;
    c = (psizr(k+1)-psizr(k-1))/2;
    b = (2*psizr(k-1)-5*psizr(k)+4*psizr(k+1)-psizr(k+2))/2;
    a = (-psizr(k-1)+3*psizr(k)-3*psizr(k+1)+psizr(k+2))/2;
    aa = 2*b.^3-9*a.*b.*c+27*a.^2.*d;
    bb = sqrt(aa.^2-4*(b.^2-3*a.*c).^3);
    xx = d./(psizr(k)-psizr(k+1));

    q1 = ((aa+bb)/2).^(1/3);
    q2 = ((aa-bb)/2).^(1/3);

    j = find(imag(q1+q2)>1e-6);	q2(j) = R13*q2(j);
    j = find(imag(q1+q2)>1e-6);	q2(j) = R13*q2(j);

    x1 = -(b+real(q1+q2))/3./a;
    x2 = -(b-(1+i*sqrt(3))/2*q1-(1-i*sqrt(3))/2*q2)/3./a;
    x3 = -(b-(1-i*sqrt(3))/2*q1-(1+i*sqrt(3))/2*q2)/3./a;

    for j = 1:length(xx)
      if     abs(imag(x1(j))) < 1e-6 & x1(j) > 0 & x1(j) < 1
	xx(j) = abs(x1(j));
      elseif abs(imag(x2(j))) < 1e-6 & x2(j) > 0 & x2(j) < 1
	xx(j) = abs(x2(j));
      elseif abs(imag(x3(j))) < 1e-6 & x3(j) > 0 & x3(j) < 1
	xx(j) = abs(x3(j));
      end
    end

    rbbbs = [rbbbs; rgg(k)];
    zbbbs = [zbbbs; zgg(k)+xx*dz];
    ib = [ib; k];
    flag_for_r_or_z_edge = [flag_for_r_or_z_edge; k*0+2];

    k = find(zbbbs > zbot & zbbbs < ztop);
    rbbbs = rbbbs(k);
    zbbbs = zbbbs(k);
    ib = ib(k);
    flag_for_r_or_z_edge = flag_for_r_or_z_edge(k);

    rbbbs = [rbbbs; r0];
    zbbbs = [zbbbs; z0];
    ib = [ib; ii0(6)];
    flag_for_r_or_z_edge(end+1) = 0;

    [thbbbs, k] = sort(angle(rbbbs-rmaxis+i*(zbbbs-zmaxis)));
    rbbbs = rbbbs(k);
    zbbbs = zbbbs(k);
    ib = ib(k);
    ibdef = find(k==max(k));
    flag_for_r_or_z_edge = flag_for_r_or_z_edge(k);

    rbbbs(end+1) = rbbbs(1);
    zbbbs(end+1) = zbbbs(1);
    thbbbs(end+1) = thbbbs(1)+2*pi;
    ib(end+1) = ib(1);
    flag_for_r_or_z_edge(end+1) = flag_for_r_or_z_edge(1);
    rhobbbs = sqrt((rbbbs-rmaxis).^2+(zbbbs-zmaxis).^2);
    nbbbs = length(rbbbs);

    for j = 1:nbbbs
      k = ib(j);
      iip = k+neighbors';
      tr = (rbbbs(j)-rgg(k))/dr;
      tz = (zbbbs(j)-zgg(k))/dz;
      if flag_for_r_or_z_edge(j) == 1 % Save a bit of time by only calculating the needed derivative
	dpsibbbsdr(j,1) = reshape(mx'*[1 tz tz^2 tz^3]'*[0 1 2*tr 3*tr^2]/dr*mx,1,16)*psizr(iip);
      else
	dpsibbbsdz(j,1) = reshape(mx'*[0 1 2*tz 3*tz^2]'/dz*[1 tr tr^2 tr^3]*mx,1,16)*psizr(iip);
      end
      wb(j,1:16) = reshape(mx'*[1 tz tz^2 tz^3]'*[1 tr tr^2 tr^3]*mx,1,16);
    end

    [rbbbsmin, irbbbsmin] = min(rbbbs(1:nbbbs-1));
    [rbbbsmax, irbbbsmax] = max(rbbbs(1:nbbbs-1));
    [zbbbsmin, izbbbsmin] = min(zbbbs(1:nbbbs-1));
    [zbbbsmax, izbbbsmax] = max(zbbbs(1:nbbbs-1));
    thgg = angle(rgg-rmaxis+i*(zgg-zmaxis));
    rhbg = spline(thbbbs(1:nbbbs-1),rhobbbs(1:nbbbs-1),thgg);
    rhgg = sqrt((rgg-rmaxis).^2+(zgg-zmaxis).^2);

    % Find grid points in the plasma, and which are new, and which were removed
    iplasma = find(rhgg<=rhbg+dr & psibarzr<=1 & zgg<=max(zbbbs) & zgg>=min(zbbbs));
    ivac = setdiff(1:ngg,iplasma);
    np = length(iplasma);

    x = psibarzr(iplasma);

    iknot = ceil(interp1(psikn,0:nkn,x)); % indices to spline knots
    iknot(x<=0) = 1;
    iknot(x>1) = nkn;

    psigrid = linspace(psimag,psibry,nr)';

    p0 = c0*sp; p1 = c1*sp; p2 = c2*sp; p3 = c3*sp;
    pres = p0(iknotg)+p1(iknotg).*psibar+p2(iknotg).*psibar.^2+p3(iknotg).*psibar.^3; pres(end) = 0;
    pprime = (p1(iknotg)+p2(iknotg)*2.*psibar+p3(iknotg)*3.*psibar.^2)*twopi/(psibry-psimag);

    f0 = c0*sf; f1 = c1*sf; f2 = c2*sf; f3 = c3*sf;
    halffpolsquared = f0(iknotg)+f1(iknotg).*psibar+f2(iknotg).*psibar.^2+f3(iknotg).*psibar.^3+rzero^2*bzero^2/2;
    k = find(halffpolsquared < 0); halffpolsquared(k) = 0;
    fpol = sqrt(2*halffpolsquared);
    ffprim = (f1(iknotg)+f2(iknotg)*2.*psibar+f3(iknotg)*3.*psibar.^2)*twopi/(psibry-psimag);

    Pprime(iplasma) = (p1(iknot)+p2(iknot)*2.*x+p3(iknot)*3.*x.^2)*twopi/(psibry-psimag);
    Pprime(ivac) = 0;

    FFprim(iplasma) = (f1(iknot)+f2(iknot)*2.*x+f3(iknot)*3.*x.^2)*twopi/(psibry-psimag);
    FFprim(ivac) = 0;

    jphi = (rgg.*Pprime+FFprim/mu0./rgg)/1e6;

    % Calculate what fraction of cells are covered by plasma and how that responds to dpsizr
    irbbbs = (rbbbs-rg(1))/dr+1; % rbbbs in grid number units
    izbbbs = (zbbbs-zg(1))/dz+1; % zbbbs in grid number units
    nxe = 0; % Index of last calculated crossing of boundary into new cell
    a = []; b = []; c = []; d = [];
    for j = 1:nbbbs-1
      sz = sign(izbbbs(j+1)-izbbbs(j));
      kz = round(izbbbs(j))+sz*0.5:sz:round(izbbbs(j+1)); % z edges in grid number units
      lkz = length(kz);
      ikz = 0; % Next index in kz to process
      sr = sign(irbbbs(j+1)-irbbbs(j));
      kr = round(irbbbs(j))+sr*0.5:sr:round(irbbbs(j+1)); % r edges in grid number units
      lkr = length(kr);
      ikr = 0; % Next index in kr to process
      while ikr<lkr | ikz<lkz % There are edges to cross
	nxe = nxe+1;
	if ikr < lkr
	  f2r = (kr(ikr+1)-irbbbs(j))/(irbbbs(j+1)-irbbbs(j)); % Fraction of next re,ze that comes from rbbbs(j+1),zbbbs(j+1)
	else
	  f2r = inf;
	end
	if ikz < lkz
	  f2z = (kz(ikz+1)-izbbbs(j))/(izbbbs(j+1)-izbbbs(j)); % Fraction of next re,ze that comes from rbbbs(j+1),zbbbs(j+1)
	else
	  f2z = inf;
	end
	if f2r < f2z % Next crossing is over r edge
	  ikr = ikr+1;
	  f2e(nxe) = f2r;
	  fe(nxe) = 2*sr; % Flag that we crossed over r edge
	else % Next crossing is over z edge
	  ikz = ikz+1;
	  f2e(nxe) = f2z;
	  fe(nxe) = sz; % Flag that we crossed over z edge
	end
	R = (1-f2e(nxe))*rbbbs(j)+f2e(nxe)*rbbbs(j+1);
	Z = (1-f2e(nxe))*zbbbs(j)+f2e(nxe)*zbbbs(j+1);
	% calculate correction of the coordinate R,Z
	kr0 = min(nr-3,max(1,floor((R-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
	kz1 = min(nz-2,max(2,ceil((Z-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
	k = kr0*nz+kz1;
	ie(nxe) = k; % Remember indices for grid points that are inside & below points (re,ze)
	iip = k+neighbors;
	pp = psizr(iip)';
	if f2r < f2z % Next crossing is over r edge, so find exactly right z position
	  a(nxe,1) = v1z*pp;
	  b(nxe,1) = v2z*pp;
	  c(nxe,1) = v3z*pp;
	  d(nxe,1) = v4z*pp-psibry;
	  re(nxe,1) = R;
	else % Next crossing is over z edge, so find exactly right r position
	  a(nxe,1) = v1r*pp;
	  b(nxe,1) = v2r*pp;
	  c(nxe,1) = v3r*pp;
	  d(nxe,1) = v4r*pp-psibry;
	  ze(nxe,1) = Z;
	end
      end
    end
    ie = ie(1:nxe);

    aa = 2*b.^3-9*a.*b.*c+27*a.^2.*d;
    bb = sqrt(aa.^2-4*(b.^2-3*a.*c).^3);

    q1 = ((aa+bb)/2).^(1/3);
    q2 = ((aa-bb)/2).^(1/3);

    j = find(imag(q1+q2)>1e-6);	q2(j) = R13*q2(j);
    j = find(imag(q1+q2)>1e-6);	q2(j) = R13*q2(j);

    x1 = -(b+real(q1+q2))/3./a;
    x2 = -(b-(1+i*sqrt(3))/2*q1-(1-i*sqrt(3))/2*q2)/3./a;
    x3 = -(b-(1-i*sqrt(3))/2*q1-(1+i*sqrt(3))/2*q2)/3./a;

    xx = ones(nxe,1)/2;;
    for j = 1:length(xx)
      if     abs(imag(x1(j))) < 1e-6 & x1(j) > -0.5 & x1(j) < 1.5
	xx(j) = abs(x1(j));
      elseif abs(imag(x2(j))) < 1e-6 & x2(j) > -0.5 & x2(j) < 1.5
	xx(j) = abs(x2(j));
      elseif abs(imag(x3(j))) < 1e-6 & x3(j) > -0.5 & x3(j) < 1.5
	xx(j) = abs(x3(j));
      end
    end

    for j = 1:nxe
      pp = psizr(ie(j)+neighbors');
      if abs(fe(j)) == 1
	tz = 0.5;
	tr = xx(j);
	re(j) = rgg(ie(j))+tr*dr;
      else
        tr = 0.5;
        tz = xx(j);
	ze(j) = zgg(ie(j))+tz*dz;
      end
      dpsiedr(j,1) = reshape(mx'*[1 tz tz^2 tz^3]'*[0 1 2*tr 3*tr^2]/dr*mx,1,16)*pp;
      dpsiedz(j,1) = reshape(mx'*[0 1 2*tz 3*tz^2]'/dz*[1 tr tr^2 tr^3]*mx,1,16)*pp;
      we(j,1:16) = reshape(mx'*[1 tz tz^2 tz^3]'*[1 tr tr^2 tr^3]*mx,1,16);
    end

    % Check for neighboring positions in rec,zec where the order in poloidal angle switched because of the corrections
    dthe = diff(angle(re(1:nxe)-rmaxis+i*(ze(1:nxe)-zmaxis)));
    k = find(abs(dthe) < pi/2 & dthe < 0);
    fec = fe;
    for j = k'
      dum = [re(j) ze(j) fe(j) ie(j) dpsiedr(j) dpsiedz(j) we(j,:)];
      re(j) = re(j+1);
      ze(j) = ze(j+1);
      fe(j) = fe(j+1);
      ie(j) = ie(j+1);
      dpsiedr(j) = dpsiedr(j+1);
      dpsiedz(j) = dpsiedz(j+1);
      we(j,:) = we(j+1,:);
      re(j+1) = dum(1);
      ze(j+1) = dum(2);
      fe(j+1) = dum(3);
      ie(j+1) = dum(4);
      dpsiedr(j+1) = dum(5);
      dpsiedz(j+1) = dum(6);
      we(j+1,:) = dum(7:22);
    end

    % At this point we have: 
    %   re, ze, coordinates where boundary goes into new grid rectangle
    %   fe, flags +1=entered up (+z) into rectange, -1=down (-z) into, +2=outward (+r) into, -2=inward (-r) into grid rectange

    % Calculate index into collapsed grid of cells that are exited at re, ze
      if fe(nxe) == -1
	icell(nxe) = 1+ceil((ze(nxe)-zg(1))/dz)+round((re(nxe)-rg(1))/dr)*nz;
      elseif fe(nxe) == 2
	icell(nxe) = 1+round((ze(nxe)-zg(1))/dz)+floor((re(nxe)-rg(1))/dr)*nz;
      elseif fe(nxe) == +1
	icell(nxe) = ceil((ze(nxe)-zg(1))/dz)+round((re(nxe)-rg(1))/dr)*nz;
      elseif fe(nxe) == -2
	icell(nxe) = 1+round((ze(nxe)-zg(1))/dz)+ceil((re(nxe)-rg(1))/dr)*nz;
      end
    jm = nxe;
    for j = 1:nxe
      if     fe(jm) == +1
	icell(j) = icell(jm)+1;
      elseif fe(jm) == -1
	icell(j) = icell(jm)-1;
      elseif fe(jm) == +2
	icell(j) = icell(jm)+nz;
      elseif fe(jm) == -2
	icell(j) = icell(jm)-nz;
      end
      jm = j;
    end
    icell = icell(1:nxe);

    % Calculate coordinates rc, zc of grid rectangle corners that are inside the plasma and share r or z with re, ze
    for j = 1:nxe
      if     fe(j) == +1
	rc(j) = rg(1+round((re(j)-rg(1))/dr))-dr/2;
	zc(j) = ze(j);
      elseif fe(j) == -1
	rc(j) = rg(1+round((re(j)-rg(1))/dr))+dr/2;
	zc(j) = ze(j);
      elseif fe(j) == +2
	zc(j) = zg(1+round((ze(j)-zg(1))/dz))+dz/2;
	rc(j) = re(j);
      elseif fe(j) == -2
	zc(j) = zg(1+round((ze(j)-zg(1))/dz))-dz/2;
	rc(j) = re(j);
      end
    end

    % Calculate Ae = areas of grid cells that are covered by plasma, and their derivatives w.r.t. re, ze
    % Here dAedre1 is dAe(j) w.r.t. re(j-1) and dAedre2 is dAe(j) w.r.t. dre(j)
    Ae = zeros(1,nxe);
    dAedre1 = zeros(1,nxe);
    dAedre2 = zeros(1,nxe);
    dAedze1 = zeros(1,nxe);
    dAedze2 = zeros(1,nxe);
    for j = 1:nxe
      if j > 1, jm = j-1; else jm = nxe; end
      if abs(fe(j)) == 1 & abs(fe(jm)) == 1 % Coming in and out over a z edge
	Ae(j) = abs(re(j)+re(jm)-rc(j)-rc(jm))*dz/2;
	dAedre1(j) = dz*sign(re(jm)-rc(jm))/2;
	dAedre2(j) = dz*sign(re(j)-rc(j))/2;
       elseif abs(fe(j)) == 2 & abs(fe(jm)) == 2 % Coming in and out over an r edge
	Ae(j) = abs(ze(j)+ze(jm)-zc(j)-zc(jm))*dr/2;
	dAedze1(j) = dr*sign(ze(jm)-zc(jm))/2;
	dAedze2(j) = dr*sign(ze(j)-zc(j))/2;
      else % Coming in over r or z edge and leaving over the other kind
	if norm([rc(j)-rc(jm) zc(j)-zc(jm)]) < 1e-12 % If they share the same corner, we have a triangle of plasma
	  Ae(j) = norm([re(j)-rc(j) ze(j)-zc(j)])*norm([re(jm)-rc(jm) ze(jm)-zc(jm)])/2;
	  if abs(fe(j)) == 1 % re(j), ze(j) is on a z edge, so re(jm), ze(jm) is an r edge
	    dAedze1(j) = abs(re(j)-rc(j))*sign(ze(jm)-zc(jm))/2;
	    dAedre2(j) = abs(ze(jm)-zc(jm))*sign(re(j)-rc(j))/2;
	  else % re(j), ze(j) is on an r edge
	    dAedre1(j) = abs(ze(j)-zc(j))*sign(re(jm)-rc(j))/2;
	    dAedze2(j) = abs(re(jm)-rc(j))*sign(ze(j)-zc(j))/2;
	  end
	else % We have a triangle of vacuum, calculate the common corner rcvac, zcvac
	  if abs(fe(j)) == 1 % Shift corner j in R direction
            rcvac = rc(j)+dr*sign(re(j)-rc(j));
	    zcvac = zc(j);
	    dAedze1(j) = -abs(re(j)-rcvac)*sign(ze(jm)-zcvac)/2;
	    dAedre2(j) = -abs(ze(jm)-zcvac)*sign(re(j)-rcvac)/2;
	  else % Shift corner j in R direction
            rcvac = rc(j);
	    zcvac = zc(j)+dz*sign(ze(j)-zc(j));
	    dAedre1(j) = -abs(ze(j)-zcvac)*sign(re(jm)-rcvac)/2;
	    dAedze2(j) = -abs(re(jm)-rcvac)*sign(ze(j)-zcvac)/2;
	  end
	  Ae(j) = dr*dz-norm([re(j)-rcvac ze(j)-zcvac])*norm([re(jm)-rcvac ze(jm)-zcvac])/2;
	end
      end
    end

    % Current density modification due to partial coverage of grid cell [MA/m2]
    jphixe = zeros(nz,nr);

    % How current is modified within rectangles due to boundary motion due to change of psizr
    djphixedpsizr = zeros(nxe,ngg); % [MA/m2]/Wb
    djphixeds = zeros(nxe,3);

    Ac = Ae; % Will become area correction due to partial cell coverage, so negative for cells where jphi ~= 0

    for j = 1:nxe
      if j > 1, jm = j-1; else jm = nxe; end
      R = rgg(icell(j));
      je = (R*pprime(end)+ffprim(end)/mu0/R);
      if jphi(icell(j)) ~= 0 % Subtract current according to fraction of vacuum within cell
	Ac(j) = Ae(j)-dr*dz;
	jj = (Ae(j)-dr*dz)/dr/dz*je/1e6; % How to modify jphi
        djphixeds(j,1) = (Ae(j)-dr*dz)/dr/dz/1e6*(p1_0(end)+2*p2_0(end)+3*p3_0(end))*twopi/(psibry-psimag)*R;
        djphixeds(j,2) = (Ae(j)-dr*dz)/dr/dz/1e6*(f1_0(end)+2*f2_0(end)+3*f3_0(end))*twopi/(psibry-psimag)/mu0/R;
        djphixeds(j,3) = (Ae(j)-dr*dz)/dr/dz/1e6*(g1_0(end)+2*g2_0(end)+3*g3_0(end))*twopi/(psibry-psimag)/mu0/R;
      else
	jj = Ae(j)/dr/dz*je/1e6; % How to modify jphi
        djphixeds(j,1) = Ae(j)/dr/dz/1e6*(p1_0(end)+2*p2_0(end)+3*p3_0(end))*twopi/(psibry-psimag)*R;
        djphixeds(j,2) = Ae(j)/dr/dz/1e6*(f1_0(end)+2*f2_0(end)+3*f3_0(end))*twopi/(psibry-psimag)/mu0/R;
        djphixeds(j,3) = Ae(j)/dr/dz/1e6*(g1_0(end)+2*g2_0(end)+3*g3_0(end))*twopi/(psibry-psimag)/mu0/R;
      end
      jphixe(icell(j)) = jj;
      iie = ie(j)+neighbors;
      djphixedpsizr(j,iie) = -je*(dAedre2(j)/dpsiedr(j)+dAedze2(j)/dpsiedz(j))*we(j,:)/dr/dz/1e6;
      iie = ie(jm)+neighbors;
      djphixedpsizr(j,iie) = djphixedpsizr(j,iie)-je*(dAedre1(j)/dpsiedr(jm)+dAedze1(j)/dpsiedz(jm))*we(jm,:)/dr/dz/1e6;
      djphixedpsizr(j,ii0) = djphixedpsizr(j,ii0)+jj*w0/(psimag-psibry)+...
	je*(dAedre1(j)/dpsiedr(jm)+dAedze1(j)/dpsiedz(jm)+dAedre2(j)/dpsiedr(j)+dAedze2(j)/dpsiedz(j))*w0/dr/dz/1e6;
      djphixedpsizr(j,iia) = djphixedpsizr(j,iia)-jj*wa/(psimag-psibry);
    end

    jphitot = jphi+jphixe;
    pcurrt = 1e6*dr*dz*jphitot;
    psipla = psizr(:)'*jphitot(:)/sum(sum(jphitot));

    % FLUX ERROR
    psizr_app = reshape(mpc*u.cc0t+mpv*u.vc0t,nz,nr);
    psizr_pla = gscalc(jphitot,rg,zg,mpp);
    psizr_err = psizr-psizr_pla-psizr_app;

    Ue = psizr_err(:);
    Ue(ii) = Ue(ii)*(-2/dr^2-2/dz^2)+Ue(ii-1)/dz^2+Ue(ii+1)/dz^2+...
      Ue(ii-nz).*(1/dr^2+1/2./rgg(ii)/dr)+Ue(ii+nz).*(1/dr^2-1/2./rgg(ii)/dr);

    % Various quantities
    cpasma = 1e6*dr*dz*sum(sum(jphitot));
    Vzr(ivac) = 0;
    Vzr(iplasma) = 2*pi*dr*dz*rgg(iplasma);
    Vzr(icell) = 2*pi*rgg(icell).*Ae(1:nxe);
    Vtot = sum(sum(Vzr));
    P(iplasma) = p0(iknot)+p1(iknot).*x+p2(iknot).*x.^2+p3(iknot).*x.^3;
    P(ivac) = 0;
    Wzr = dr*dz*3*pi*rgg.*P;
    Wzr(icell) = Ae*3*pi.*rgg(icell).*P(icell);
    Wth = sum(sum(Wzr));
    psizr_r = [zeros(nz,1) psizr(:,3:end)-psizr(:,1:end-2) zeros(nz,1)]/2/dr;
    psizr_z = [zeros(1,nr);psizr(3:end,:)-psizr(1:end-2,:);zeros(1,nr)]/2/dz;
    Bp2zr = (psizr_r.^2+psizr_z.^2)/twopi^2./rgg.^2;
    DX = diff(rbbbs)';
    DY = diff(zbbbs)';
    cpieces = sqrt(DX.^2+DY.^2);
    Cl = sum(cpieces); % Contour length
    bp2flx = (mu0*cpasma/Cl)^2;
    Bp2V = sum(sum(Vzr.*Bp2zr));      
    betap = 4/3*mu0*Wth/Vtot/bp2flx;
    li = Bp2V/Vtot/bp2flx; % Volume-averaged Bp^2 / bp2flx


    % End of equilibrium analysis for the sake of debugging

  
  
  
  
  
    clear ds da

    ds.psizr = 'Poloidal flux on grid [Wb]';
    eq.psizr = psizr;

    ds.cc = 'PF coil currents in same units as the template that was used to design this equilibrium';
    eq.cc = piccc*ic;

    ds.ic = 'PF coil currents in Toksys units of terminal Amperes';
    eq.ic = ic;

    ds.iv = 'Vessel currents in Toksys units';
    eq.iv = iv;

    ds.cpasma = 'Plasma current [A]';
    eq.cpasma = cpasma;

    ds.jphi = 'Toroidal current density on grid [MA/m^2]';  
    eq.jphi = jphi;

    ds.pcurrt = 'Toroidal current within grid rectangles including partially covered cells [A]';  
    eq.pcurrt = pcurrt;

    ds.nbbbs = 'Number of boundary points';
    eq.nbbbs = nbbbs;

    ds.rbbbs = 'Radius of boundary points [m]';
    eq.rbbbs = rbbbs(nbbbs:-1:1);

    ds.zbbbs = 'Height of boundary points [m]';
    eq.zbbbs = zbbbs(nbbbs:-1:1);

    for j = 1:nbbbs
      dbbbs(j,1) = sqrt(min((rl-rbbbs(j)).^2+(zl-zbbbs(j)).^2));
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

    ds.r0 = 'Radius of the point that defines the boundary [m]';
    eq.r0 = r0;

    ds.z0 = 'Height of the point that defines the boundary [m]';
    eq.z0 = z0;

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

    ds.betap = 'Normalized plasma pressure';
    eq.betap = betap;
    
    ds.Vzr = 'Volumes of plasma within grid elements';
    eq.Vzr = Vzr;

    ds.Vtot = 'Total plasma volume';
    eq.Vtot = Vtot;

    ds.P = 'Pressure on grid elements';
    eq.P = P;

    ds.Wzr = 'Thermal energy within grid elements';
    eq.Wzr = Wzr;

    ds.Wth = 'Total thermal energy';
    eq.Wth = Wth;

    ds.Cl = 'Length of boundary contour';
    eq.Cl = Cl;
    
    ds.Bp2zr = 'square of poloidal field within grid elements';
    eq.Bp2zr = Bp2zr;

    ds.bp2flx = '(mu0*Ip/cpasma)^2';
    eq.bp2flx = bp2flx;

    ds.Bp2V = 'sum(sum(Vzr.*Bp2zr))';
    eq.Bp2V = Bp2V;
    
    
    
    eq.icell = icell;
    eq.re = re;
    eq.ze = ze;
    eq.Ae = Ae;
    
    
    

    ds.psirz = '-psizr''/2/pi, efitviewer plots contours with this variable';
    eq.psirz = -psizr'/2/pi;

    ds.ssibry = '-psibry/2/pi, efitviewer plots contours with this variable';
    eq.ssibry = -psibry/2/pi;

    ds.ssimag = '-psimag/2/pi, efitviewer plots contours with this variable';
    eq.ssimag = -psimag/2/pi;

    ds.rx1 = 'Radius of closest x point found below plasma [m]';
    eq.rx1 = rx1;

    ds.zx1 = 'Height of closest x point found below plasma [m]';
    eq.zx1 = zx1;

    ds.psix1 = 'Flux at x1 [Wb]';
    eq.psix1 = psix1;

    ds.rx2 = 'Radius of closest x point found above plasma [m]';
    eq.rx2 = rx2;

    ds.zx2 = 'Height of closest x point found above plasma [m]';
    eq.zx2 = zx2;

    ds.psix2 = 'Flux at x2 [Wb]';
    eq.psix2 = psix2;
    
    rsurf = (min(rbbbs)+max(rbbbs))/2;
    ds.rsurf = 'Radius of geometric center of boundary [m]';
    eq.rsurf = rsurf;

    ds.rsurf = 'Height of geometric center of boundary [m]';
    eq.rsurf = (min(zbbbs)+max(zbbbs))/2;

    aminor = (max(rbbbs)-min(rbbbs))/2;
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

    ds.doutl = 'Upper triangularity';
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
    
    th0 = angle(r0-rmaxis + i*(z0-zmaxis));
    da.limloc = 'Shape descriptor, possible values are: IN , OUT, TOP, BOT, DN , SNT, SNB, MAR';
    if ilimited
      if     th0*180/pi > -45 & th0*180/pi <=  45
	ea.limloc = 'OUT';
      elseif th0*180/pi >  45 & th0*180/pi <= 135
	ea.limloc = 'TOP';
      elseif th0*180/pi >  -135 & th0*180/pi <= -45
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

    k = find(thbbbs*180/pi > 135 | thbbbs*180/pi < -135);
    da.oleft = 'Smallest distance to limiter from boundary sector 135<theta<225 [cm]';
    ea.oleft = 100*min(dbbbs(k));

    k = find(thbbbs*180/pi > -45 & thbbbs*180/pi < 45);
    da.oright = 'Smallest distance to limiter from boundary sector -45<theta<45 [cm]';
    ea.oright = 100*min(dbbbs(k));

    k = find(thbbbs*180/pi > 45 & thbbbs*180/pi < 135);
    da.otop = 'Smallest distance to limiter from boundary sector 45<theta<135 [cm]';
    ea.otop = 100*min(dbbbs(k));

    k = find(thbbbs*180/pi > -135 & thbbbs*180/pi < -45);
    da.obott = 'Smallest distance to limiter from boundary sector -135<theta<-45 [cm]';
    ea.obott = 100*min(dbbbs(k));

    da.rseps1 = 'Radial position of lower x-point [cm]';
    ea.rseps1 = eq.rx1*100;

    da.zseps1 = 'Vertical position of lower x-point [cm]';
    ea.zseps1 = eq.zx1*100;

    da.rseps2 = 'Radial position of upper x-point [cm]';
    ea.rseps2 = eq.rx2*100;

    da.zseps2 = 'Vertical position of upper x-point [cm]';
    ea.zseps2 = eq.zx2*100;

    da.betapd = 'Diamagnetic betap';
    ea.betapd = betap;

    da.wplasmd = 'Diamagnetic thermal energy [J]';
    ea.wplasmd = Wth;

    da.wplasm = 'Thermal energy [J]';
    ea.wplasm = Wth;

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
      if ~isempty(eval(s))
	ds = setfield(ds,s,ss);
	eq = setfield(eq,s,eval(s));
      end
    end
    if ~isfield(eq,'ecid')
      eq.ecid = [];
    end

    ds.psizr_err = 'flux error, psizr_err = psizr-psizr_pla-psizr_app;';
    eq.psizr_err = psizr_err;
    
    eq.descriptions = ds;
  
    % Restore last analyzed equilibrium
    cpasma    = eqx.cpasma;
    rmaxis    = eqx.rmaxis;
    zmaxis    = eqx.zmaxis;
    psimag    = eqx.psimag;
    r0        = eqx.r0;
    z0        = eqx.z0;
    psibry    = eqx.psibry;
    x1exists  = eqx.x1exists;
    x1inside  = eqx.x1inside;
    rx1       = eqx.rx1;
    zx1       = eqx.zx1;
    x2exists  = eqx.x2exists;
    x2inside  = eqx.x2inside;
    rx2       = eqx.rx2;
    zx2       = eqx.zx2;
    rzero     = eqx.rzero;
    bzero     = eqx.bzero;
    Pprime    = eqx.Pprime;
    FFprim    = eqx.FFprim;
    P         = eqx.P;
    iplasma   = eqx.iplasma;
    rbbbs     = eqx.rbbbs;
    zbbbs     = eqx.zbbbs;
    zbot      = eqx.zbot;
    ztop      = eqx.ztop;
    psipla    = eqx.psipla;
    li        = eqx.li;
    betap     = eqx.betap;
    jphi      = eqx.jphi;
    pcurrt    = eqx.pcurrt;
    ii0       = eqx.ii0;
    w0        = eqx.w0;
    iia       = eqx.iia;
    wa        = eqx.wa;
    ilimited  = eqx.ilimited;
    ix1       = eqx.ix1;
    ix2       = eqx.ix2;
    nbbbs     = eqx.nbbbs;
    fpol      = eqx.fpol;
    ffprim    = eqx.ffprim;
    pres      = eqx.pres;
    pprime    = eqx.pprime;
    psix1     = eqx.psix1;
    psix2     = eqx.psix2;
    thbbbs    = eqx.thbbbs;
    Wth       = eqx.Wth;
    psizr_err = eqx.psizr_err;
    Vzr       = eqx.Vzr;
    Bp2zr     = eqx.Bp2zr;
    Cl        = eqx.Cl;
    Bp2V      = eqx.Bp2V;
    Vtot      = eqx.Vtot;
    bp2flx    = eqx.bp2flx;
    psibarzr  = eqx.psibarzr;
    
  end
  
  
  
  
  
  
  
  
  
  
  
  
