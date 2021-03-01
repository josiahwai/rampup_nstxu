%  USAGE:   gs_profiles
%
%  PURPOSE: Find how much of quantities are contained within flux-
%           surfaces from axis to boundary, by accumulated sums
%           of quantities within grid elements sorted w.r.t. psibarzr
%
%  INPUTS:   profile_fit_option (default = 2)
%              0 = just fit cumsum of psibarzr-sorted grid cells
%              1 = fit equals the calculated total value at psibar = 1
%              2 = fit also equals calculated derivative at psibar = 0
%            calculate_profile_responses, a flag, default = false
%
%            Remaining inputs are made by gs_eq_analysis:
%            psibarzr, normalized flux on the grid
%            p0, p1, p2, p3, spline function for pressure
%            Wzr, thermal energy within grid cells
%            Wth, total thermal energy
%            Vtot, total plasma volume
%            pprimezr, pprime on the grid
%            ffprimzr, ffprim on the grid
%            cpasma, total plasma current
%            Atot, total plasma area
%            fpolzr, fpol on the grid
%            torflux, total toroidal flux within plasma
%
%  OUTPUTS: SPLINE FUNCTIONS (poly-coefficients for each spline region):
%             v0, v1, v2, v3, spline function for volume, V(psibar)
%             w0, w1, w2, w3, w4, w5, w6, thermal energy, W(psibar)
%             a0, a1, a2, a3, spline function for area, A(psibar)
%             i0, i1, i2, i3, for toroidal current, I(psibar)
%             t0, t1, t2, t3, for toroidal flux, T(psibar)
%           Values at psibar (i.e. nr values from axis to boundary):
%             V, volume contained within flux surfaces
%             W, thermal energy contained within flux surfaces
%             A, area contained within flux surfaces
%             I, toroidal current contained within flux surfaces
%             T, toroidal flux contained within flux surfaces
%             jtav, contour-averaged current density (dI/dA)
%             rhot, square-root of normalized toroidal flux, sqrt(T/T(nr))
%             qpsi, q values (dT/dpsi)
%	
%  METHOD: The grid points inside the plasma are sorted w.r.t. psibarzr.
%          Quantities contained within flux surfaces are found by
%          accumulated sums over the sorted grid points.
%          Spline functions are fitted to these accumulated sums.

%  NOTES:  To do: L consistent with fpol & T, responses to L, V, W,
%                 fits with Lpa, Vpa for L, V
	
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander  ON	4/13/13
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('profile_fit_option','var') | length(profile_fit_option) ~= 1
  profile_fit_option = 2;
end
if ~exist('calculate_profile_responses','var') | length(calculate_profile_responses) ~= 1
  calculate_profile_responses = false;
end

if ~maxis
  sb = zeros(nkn+2,1);
  dsbdx = zeros(nkn+2,nx);
  sh = zeros(nkn+2,1);
  dshdx = zeros(nkn+2,nx);
  return
end


% Variables unique to gs_profiles
ygg = zeros(ngg,1); % Sorted psibarzr for iplasma
kgg  = ones(ngg,1); % indices for sorting psibarzr(iplasma)
Wgg = zeros(ngg,1); % Wzr values sorted using kgg
Igg = zeros(ngg,1); % pcurrt values sorted using kgg
Tgg = zeros(ngg,1); % psitor values sorted using kgg
ikg =  ones(ngg,1); % iknot sorted using kgg
Bgg = zeros(ngg,1);
Cgg = zeros(ngg,1);
Hgg = zeros(ngg,1);

% Needed for all profiles
np = sum(iplasma(:)); % number of grid points that are covered by plasma

% No plasma case
if np == 0
  return
end

[ygg(1:np),kgg(1:np)] = sort(psibarzr(iplasma)); % sorted psibarzr values
dumzr(1:np) = iknot(iplasma);
ikg(1:np) = dumzr(kgg(1:np));
% Use M (which is ngg*ngg in size) as storage space
M(1:np,1:nkn+2) = [             d0(ikg(1:np),:) + ...
    ygg(1:np)   *ones(1,nkn+2).*d1(ikg(1:np),:) + ...
    ygg(1:np).^2*ones(1,nkn+2).*d2(ikg(1:np),:) + ...
    ygg(1:np).^3*ones(1,nkn+2).*d3(ikg(1:np),:)];
MTMi11 = inv(M(1:np,1:nkn+1)'*M(1:np,1:nkn+1));
MTMi12 = inv(M(1:np,1:nkn+2)'*M(1:np,1:nkn+2));
MTMi21 = inv(M(1:np,2:nkn+1)'*M(1:np,2:nkn+1));
% Elliptic surfaces near axis with a, b almost exactly along R, Z
psimagrr = warr*psizr(iia);
psimagzz = wazz*psizr(iia);
% The coordinate system should be a tiny bit rotated for non-zero psimagrz
%psimagrz = warz*psizr(iia); % Not doing this rotation yet
%alpha = asin(psimagrz/psimagrr)/2
%cos(alpha)^2*psimagrr+2*cos(alpha)*sin(alpha)*psimagrz+sin(alpha)^2*psimagzz
%cos(alpha)^2*psimagrr+2*cos(alpha)*sin(alpha)*psimagrz+sin(alpha)^2*psimagzz
if calculate_profile_responses
  dyggdx = zeros(ngg,nx);
  dpsibarzrdx = (dpsizrdx-(1-psibarzr(:))*dpsimagdx-psibarzr(:)*dpsibrydx)/(psibry-psimag);
  dyggdx(1:np,:) = dpsibarzrdx(iplasma,:);
  dyggdx(1:np,:) = dyggdx(kgg(1:np),:);
  ps2(1:np,:) = [                 d1(ikg(1:np),1:nkn+2) + ...
    2*ygg(1:np)   *ones(1,nkn+2).*d2(ikg(1:np),1:nkn+2) + ...
    3*ygg(1:np).^2*ones(1,nkn+2).*d3(ikg(1:np),1:nkn+2)];
  psibarnx = psibar*ones(1,nx);
  tr = (rmaxis-rgg(iia(6)))/dr;
  tz = (zmaxis-zgg(iia(6)))/dz;
  warrr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 0 6]/dr^3*mx,1,16);
  warrz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[0 0 2 6*tr]/dr^2*mx,1,16);
  warzz = reshape(([0 0 2 6*tz]/dz^2*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
  wazzz = reshape(([0 0 0 6]/dz^3*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
  psimagrrr = warrr*psizr(iia);
  psimagrrz = warrz*psizr(iia);
  psimagrzz = warzz*psizr(iia);
  psimagzzz = wazzz*psizr(iia);
  dpsimagrrdx = warr*dpsizrdx(iia,:)+psimagrrr*drmaxisdx+psimagrrz*dzmaxisdx;
  dpsimagzzdx = wazz*dpsizrdx(iia,:)+psimagrzz*drmaxisdx+psimagzzz*dzmaxisdx;
  dTggdx = AR(:)*ones(1,nx).*dfpolzrdx;
  dTggdx(1:np,:) = dTggdx(iplasma,:);
  dTggdx(1:np,:) = cumsum(dTggdx(kgg(1:np),:));
  dIggdx = zeros(ngg,nx);
  dIggdx(1:np,:) = dpcurrtdx(iplasma,:);
  dIggdx(1:np,:) = cumsum(dIggdx(kgg(1:np),:));
end


% Calculate W(psibar) and V(psibar)

dumzr(1:np) = Wzr(iplasma);
Wgg(1:np) = cumsum(dumzr(kgg((1:np)))); % W within fluxes ygg

% dW/dV = W'/V' = 1.5*pres, and pres is given by a spline function.
% Find a V' that produces W' and fits Wgg
% Multiply expression for pres by expression for Vprime and integrate
% Then fit sv to Wgg using pinv

dWggdsv = zeros(ngg,nkn+2); % Response of Wgg(ygg) to sv
dWthdsv = zeros(1,nkn+2); % Response of Wth to sv (sv = "volume function")
wjump = zeros(nkn-1,nkn+2); % Jumps in W at knots
o = ones(ngg,1);
okn = ones(nkn-1,1);
mp = 1.5*[p0 p1 p2 p3];
for jk = 1:6
  o(1:np) = o(1:np).*ygg(1:np);
  okn = okn.*psikn(2:nkn)';
  for j = 1:4
    k = jk-j+1;
    if k == 1
      dWggdsv(1:np,:) = dWggdsv(1:np,:) + k/jk* ...
        mp(ikg(1:np),j).*o(1:np)*ones(1,nkn+2).*d1(ikg(1:np),:);
      dWthdsv(1,:) = dWthdsv(1,:) + k/jk*mp(nkn,j)*d1(nkn,:);
      wjump = wjump + k/jk*(mp(2:nkn,j).*okn*ones(1,nkn+2).*d1(2:nkn,:) - ...
                      mp(1:nkn-1,j).*okn*ones(1,nkn+2).*d1(1:nkn-1,:));
    elseif k == 2
      dWggdsv(1:np,:) = dWggdsv(1:np,:) + k/jk* ...
        mp(ikg(1:np),j).*o(1:np)*ones(1,nkn+2).*d2(ikg(1:np),:);
      dWthdsv(1,:) = dWthdsv(1,:) + k/jk*mp(nkn,j)*d2(nkn,:);
      wjump = wjump + k/jk*(mp(2:nkn,j).*okn*ones(1,nkn+2).*d2(2:nkn,:) - ...
                      mp(1:nkn-1,j).*okn*ones(1,nkn+2).*d2(1:nkn-1,:));
    elseif k == 3
      dWggdsv(1:np,:) = dWggdsv(1:np,:) + k/jk* ...
        mp(ikg(1:np),j).*o(1:np)*ones(1,nkn+2).*d3(ikg(1:np),:);
      dWthdsv(1,:) = dWthdsv(1,:) + k/jk*mp(nkn,j)*d3(nkn,:);
      wjump = wjump + k/jk*(mp(2:nkn,j).*okn*ones(1,nkn+2).*d3(2:nkn,:) - ...
                      mp(1:nkn-1,j).*okn*ones(1,nkn+2).*d3(1:nkn-1,:));
    end  
  end
end

dw0dsv = -cumsum([zeros(1,nkn+2); wjump; zeros(1,nkn+2)]);

dWggdsv(1:np,:) = dWggdsv(1:np,:)+dw0dsv(ikg(1:np),:);

dWthdsv = dWthdsv+dw0dsv(nkn,:);

if profile_fit_option == 0
  % Simple fit
  sv = dWggdsv(1:np,:)\Wgg(1:np);

  % Fit that requires V(1) = Vtot (from gs_eq_analysis)
  sv = [zeros(nkn+1,1); Vtot];
  sv(1:nkn+1) = dWggdsv(1:np,1:nkn+1)\[Wgg(1:np)-dWggdsv(1:np,nkn+2)*Vtot];
  sv(nkn+1) = (Wth-dWthdsv(1:nkn)*sv(1:nkn)-dWthdsv(nkn+2)*Vtot)/dWthdsv(nkn+1);
else
  % Fit that also requires W(1) = Wth (from gs_eq_analysis)
  sv = [zeros(nkn+1,1); Vtot];
  dum = (Wth-dWthdsv(nkn+2)*Vtot)/dWthdsv(nkn+1);
  sv(1:nkn) = (dWggdsv(1:np,1:nkn)-dWggdsv(1:np,nkn+1)*[dWthdsv(1:nkn)/dWthdsv(nkn+1)])\...
              [Wgg(1:np)-dWggdsv(1:np,nkn+2)*Vtot-dWggdsv(1:np,nkn+1)*dum];
  sv(nkn+1) = (Wth-dWthdsv(1:nkn)*sv(1:nkn)-dWthdsv(nkn+2)*Vtot)/dWthdsv(nkn+1);
end

% Spline function for volume, V
v0 = d0*sv;
v1 = d1*sv;
v2 = d2*sv;
v3 = d3*sv;
V = v0(iknotg) + v1(iknotg).*psibar + v2(iknotg).*psibar.^2 + v3(iknotg).*psibar.^3;

% Spline function for thermal energy, W
w0 = dw0dsv*sv; % The constants in the spline function for W
w1 = 1.5*(p0.*v1                                        );
w2 = 1.5*(p0.*v2 + p1.*v1/2                             );
w3 = 1.5*(p0.*v3 + p1.*v2/1.5  + p2.*v1/3               );
w4 = 1.5*(         p1.*v3*0.75 + p2.*v2/2   + p3.*v1/4  );
w5 = 1.5*(                       p2.*v3*0.6 + p3.*v2*0.4);
w6 = 1.5*(                                    p3.*v3/2  );
W = w0(iknotg)            + ...
    w1(iknotg).*psibar    + ...
    w2(iknotg).*psibar.^2 + ...
    w3(iknotg).*psibar.^3 + ...
    w4(iknotg).*psibar.^4 + ...
    w5(iknotg).*psibar.^5 + ...
    w6(iknotg).*psibar.^6;


% Calculate current within flux surfaces, I(psibar)

dumzr(1:np) = pcurrt(iplasma);
Igg(1:np) = cumsum(dumzr(kgg(1:np)));

if profile_fit_option == 0
  % Simple fit
  si = M(1:np,1:nkn+2)\Igg(1:np);
  if calculate_profile_responses
    dsidx = zeros(nkn+2,nx);
    for j = 1:nx
      ps3(1:np,:) = dyggdx(1:np,j)*ones(1,nkn+2).*ps2(1:np,:);
      MTM12 = M(1:np,1:nkn+2)'*ps3(1:np,1:nkn+2)+ps3(1:np,1:nkn+2)'*M(1:np,1:nkn+2);
      dsidx(:,j) = MTMi12*(ps3(1:np,:)'*Igg(1:np) - MTM12*si + ...
        M(1:np,1:nkn+2)'*dIggdx(1:np,j));
    end
  end
elseif profile_fit_option == 1
  % Fit that requires total current to be cpasma (from gs_eq_analysis)
  Cgg(1:np) = Igg(1:np)-M(1:np,nkn+2)*cpasma;
  si = [M(1:np,1:nkn+1)\Cgg(1:np); cpasma];
  if calculate_profile_responses
    dsidx = [zeros(nkn+1,nx); dcpasmadx];
    for j = 1:nx
      ps3(1:np,:) = dyggdx(1:np,j)*ones(1,nkn+2).*ps2(1:np,:);
      MTM11 = M(1:np,1:nkn+1)'*ps3(1:np,1:nkn+1)+ps3(1:np,1:nkn+1)'*M(1:np,1:nkn+1);
      dsidx(1:nkn+1,j) = MTMi11*(ps3(1:np,1:nkn+1)'*Cgg(1:np) - ...
        MTM11*si(1:nkn+1) + M(1:np,1:nkn+1)'*dIggdx(1:np,j) - ...
	M(1:np,1:nkn+1)'*(ps3(1:np,nkn+2)*cpasma + M(1:np,nkn+2)*dcpasmadx(j)));
    end
  end
else
  % Fit that also locks Iprime on axis
  jmaxis = rmaxis*pprime(1) + ffprim(1)/mu0/rmaxis;
  Ipa = jmaxis*2*pi*(psimag-psibry)/sqrt(psimagrr*psimagzz);
  Cgg(1:np) = Igg(1:np) - M(1:np,1)*Ipa - M(1:np,nkn+2)*cpasma;
  si = [Ipa; M(1:np,2:nkn+1)\Cgg(1:np); cpasma];
  if calculate_profile_responses
    djmaxisdx = -(dpsibrydx-dpsimagdx)/(psibry-psimag)*jmaxis + ...
      drmaxisdx*pprime(1)-ffprim(1)/mu0/rmaxis^2*drmaxisdx;
    djmaxisdx(indsp) = djmaxisdx(indsp) + c1(1,:)*rmaxis;
    djmaxisdx(indsf) = djmaxisdx(indsf) + c1(1,:)/mu0/rmaxis;
    dIpadx = Ipa*(djmaxisdx/jmaxis + (dpsimagdx-dpsibrydx)/(psimag-psibry)-...
      (dpsimagrrdx*psimagzz+psimagrr*dpsimagzzdx)/(psimagrr*psimagzz)/2);
    dsidx = [dIpadx; zeros(nkn,nx); dcpasmadx];
    for j = 1:nx
      ps3(1:np,:) = dyggdx(1:np,j)*ones(1,nkn+2).*ps2(1:np,:);
      MTM21 = M(1:np,2:nkn+1)'*ps3(1:np,2:nkn+1)+ps3(1:np,2:nkn+1)'*M(1:np,2:nkn+1);
      dsidx(2:nkn+1,j) = MTMi21*(ps3(1:np,2:nkn+1)'*Cgg(1:np) - ...
        MTM21*si(2:nkn+1) + M(1:np,2:nkn+1)'*dIggdx(1:np,j) - ...
	M(1:np,2:nkn+1)'*(ps3(1:np,nkn+2)*cpasma + M(1:np,nkn+2)*dcpasmadx(j) + ...
	ps3(1:np,1)*Ipa + M(1:np,1)*dIpadx(j)));
    end
  end
end
% The function current within flux surfaces, I(psibar)
i0 = d0*si;
i1 = d1*si;
i2 = d2*si;
i3 = d3*si;
I = i0(iknotg) + i1(iknotg).*psibar + i2(iknotg).*psibar.^2 + i3(iknotg).*psibar.^3;
Iprime = 2*pi/(psibry-psimag)*(i1(iknotg) + 2*i2(iknotg).*psibar + 3*i3(iknotg).*psibar.^2);
if calculate_profile_responses
  di0dx = d0*dsidx;
  di1dx = d1*dsidx;
  di2dx = d2*dsidx;
  di3dx = d3*dsidx;
  dIdx = di0dx(iknotg,:)              + di1dx(iknotg,:).*psibarnx + ...
         di2dx(iknotg,:).*psibarnx.^2 + di3dx(iknotg,:).*psibarnx.^3;
  dIprimedx = 2*pi/(psibry-psimag)*(...
      di1dx(iknotg,:)           + ...
    2*di2dx(iknotg,:).*psibarnx + ...
    3*di3dx(iknotg,:).*psibarnx.^2) - ...
    Iprime/(psibry-psimag)*(dpsibrydx-dpsimagdx);
end


% Calculate area within flux surfaces, A(psibar)

if profile_fit_option == 0
  % Simple fit
  Cgg(1:np) = Ag*(1:np)';
  sa = M(1:np,1:nkn+2)\Cgg(1:np);
  if calculate_profile_responses
    dsadx = zeros(nkn+2,nx);
    for j = 1:nx
      ps3(1:np,:) = dyggdx(1:np,j)*ones(1,nkn+2).*ps2(1:np,:);
      MTM12 = M(1:np,1:nkn+2)'*ps3(1:np,1:nkn+2)+ps3(1:np,1:nkn+2)'*M(1:np,1:nkn+2);
      dsadx(:,j) = MTMi12*(ps3(1:np,:)'*Cgg(1:np) - MTM12*sa);
    end
  end
elseif profile_fit_option == 1
  % Fit that requires total area to be Atot (gs_eq_analysis)
  Cgg(1:np) = Ag*(1:np)' - M(1:np,nkn+2)*Atot;
  sa = [M(1:np,1:nkn+1)\Cgg(1:np); Atot];
  if calculate_profile_responses
    dAtotdx = sum(dAcelldpsizr*dpsizrdx);
    dsadx = [zeros(nkn+1,nx); dAtotdx];
    for j = 1:nx
      ps3(1:np,:) = dyggdx(1:np,j)*ones(1,nkn+2).*ps2(1:np,:);
      MTM11 = M(1:np,1:nkn+1)'*ps3(1:np,1:nkn+1)+ps3(1:np,1:nkn+1)'*M(1:np,1:nkn+1);
      dsadx(1:nkn+1,j) = MTMi11*(ps3(1:np,1:nkn+1)'*Cgg(1:np) - MTM11*sa(1:nkn+1) - ...
	M(1:np,1:nkn+1)'*(ps3(1:np,nkn+2)*Atot + M(1:np,nkn+2)*dAtotdx(j)));
    end
  end
else
  % Fit that also locks Aprime on axis
  Apa = 2*pi*(psimag-psibry)/sqrt(psimagrr*psimagzz);
  Cgg(1:np) = Ag*(1:np)' - M(1:np,1)*Apa - M(1:np,nkn+2)*Atot;
  sa = [Apa; M(1:np,2:nkn+1)\Cgg(1:np); Atot];
  if calculate_profile_responses
    dApadx = Apa*((dpsimagdx-dpsibrydx)/(psimag-psibry)-...
      (dpsimagrrdx*psimagzz+psimagrr*dpsimagzzdx)/(psimagrr*psimagzz)/2);
    dAtotdx = sum(dAcelldpsizr*dpsizrdx);
    dsadx = [dApadx; zeros(nkn,nx); dAtotdx];
    for j = 1:nx
      ps3(1:np,:) = dyggdx(1:np,j)*ones(1,nkn+2).*ps2(1:np,:);
      MTM21 = M(1:np,2:nkn+1)'*ps3(1:np,2:nkn+1)+ps3(1:np,2:nkn+1)'*M(1:np,2:nkn+1);
      dsadx(2:nkn+1,j) = MTMi21*(ps3(1:np,2:nkn+1)'*Cgg(1:np) - MTM21*sa(2:nkn+1) - ...
	M(1:np,2:nkn+1)'*(ps3(1:np,nkn+2)*Atot + M(1:np,nkn+2)*dAtotdx(j) + ...
	ps3(1:np,1)*Apa + M(1:np,1)*dApadx(j)));
    end
  end
end
% The function area within flux surfaces, A(psibar)
a0 = d0*sa;
a1 = d1*sa;
a2 = d2*sa;
a3 = d3*sa;
A = a0(iknotg) + a1(iknotg).*psibar + a2(iknotg).*psibar.^2 + a3(iknotg).*psibar.^3;
Aprime = 2*pi/(psibry-psimag)*(a1(iknotg) + 2*a2(iknotg).*psibar + 3*a3(iknotg).*psibar.^2);
jtav = Iprime./Aprime;
if calculate_profile_responses
  da0dx = d0*dsadx;
  da1dx = d1*dsadx;
  da2dx = d2*dsadx;
  da3dx = d3*dsadx;
  dAdx = da0dx(iknotg,:)              + da1dx(iknotg,:).*psibarnx + ...
         da2dx(iknotg,:).*psibarnx.^2 + da3dx(iknotg,:).*psibarnx.^3;
  dAprimedx = 2*pi/(psibry-psimag)*(...
      da1dx(iknotg,:)           + ...
    2*da2dx(iknotg,:).*psibarnx + ...
    3*da3dx(iknotg,:).*psibarnx.^2) - ...
    Aprime/(psibry-psimag)*(dpsibrydx-dpsimagdx);
  djtavdx = dIprimedx./(Aprime*ones(1,nx))-jtav./Aprime*ones(1,nx).*dAprimedx;
end


% Calculate toroidal flux within flux surfaces and its derivative q

dumzr(1:np) = AR(iplasma).*fpolzr(iplasma);
Tgg(1:np) = cumsum(dumzr(kgg(1:np)));
if profile_fit_option == 0
  % Simple fit
  st = M(1:np,1:nkn+2)\Tgg(1:np);
  if calculate_profile_responses
    dstdx = zeros(nkn+2,nx);
    for j = 1:nx
      ps3(1:np,:) = dyggdx(1:np,j)*ones(1,nkn+2).*ps2(1:np,:);
      MTM12 = M(1:np,1:nkn+2)'*ps3(1:np,1:nkn+2)+ps3(1:np,1:nkn+2)'*M(1:np,1:nkn+2);
      dstdx(:,j) = MTMi12*(ps3(1:np,:)'*Tgg(1:np) - MTM12*st + ...
        M(1:np,1:nkn+2)'*dTggdx(1:np,j));
    end
  end
elseif profile_fit_option == 1
  % Fit that requires total toroidal flux to be torflux (gs_eq_analysis)
  Cgg(1:np) = Tgg(1:np)-M(1:np,nkn+2)*torflux;
  st = [M(1:np,1:nkn+1)\Cgg(1:np); torflux];
  if calculate_profile_responses
    dstdx = [zeros(nkn+1,nx); dtorfluxdx];
    for j = 1:nx
      ps3(1:np,:) = dyggdx(1:np,j)*ones(1,nkn+2).*ps2(1:np,:);
      MTM11 = M(1:np,1:nkn+1)'*ps3(1:np,1:nkn+1)+ps3(1:np,1:nkn+1)'*M(1:np,1:nkn+1);
      dstdx(1:nkn+1,j) = MTMi11*(ps3(1:np,1:nkn+1)'*Cgg(1:np) - ...
        MTM11*st(1:nkn+1) + M(1:np,1:nkn+1)'*dTggdx(1:np,j) - ...
	M(1:np,1:nkn+1)'*(ps3(1:np,nkn+2)*torflux + M(1:np,nkn+2)*dtorfluxdx(j)));
    end
  end
else
  % Fit that also requires derivative on axis to be given by value of Tpa below
  % The ellipse has b/a = sqrt(psimagrr/psimagzz)
  % Toroidal flux = pi*a*b*Bt, Poloidal flux = psimagrr*a^2/2, q=2*pi*kappa/psimagrr
  fpol1 = fpol(1);
  Tpa = 2*pi*(psimag-psibry)*fpol(1)/rmaxis/sqrt(psimagrr*psimagzz);  
  Cgg(1:np) = Tgg(1:np)-M(1:np,1)*Tpa-M(1:np,nkn+2)*torflux;
  st = [Tpa; M(1:np,2:nkn+1)\Cgg(1:np); torflux];
  if calculate_profile_responses
    dfpol1dx = zeros(1,nx);
    dfpol1dx(:,indsf) = dfpol1dx(:,indsf)+c0(1,:)/fpol1;  
    dTpadx = Tpa*((dpsimagdx-dpsibrydx)/(psimag-psibry) + ...
      dfpol1dx/fpol1 - drmaxisdx/rmaxis - ...
      (dpsimagrrdx*psimagzz+psimagrr*dpsimagzzdx)/(psimagrr*psimagzz)/2);
    dstdx = [dTpadx; zeros(nkn,nx); dtorfluxdx];
    for j = 1:nx
      ps3(1:np,:) = dyggdx(1:np,j)*ones(1,nkn+2).*ps2(1:np,:);
      MTM21 = M(1:np,2:nkn+1)'*ps3(1:np,2:nkn+1)+ps3(1:np,2:nkn+1)'*M(1:np,2:nkn+1);
      dstdx(2:nkn+1,j) = MTMi21*(ps3(1:np,2:nkn+1)'*Cgg(1:np) - ...
        MTM21*st(2:nkn+1) + M(1:np,2:nkn+1)'*dTggdx(1:np,j) - ...
	M(1:np,2:nkn+1)'*(ps3(1:np,nkn+2)*torflux + M(1:np,nkn+2)*dtorfluxdx(j) + ...
	ps3(1:np,1)*Tpa + M(1:np,1)*dTpadx(j)));
    end
  end
end
% Toroidal flux function
t0 = d0*st;
t1 = d1*st;
t2 = d2*st;
t3 = d3*st;
T = t0(iknotg) + t1(iknotg).*psibar + t2(iknotg).*psibar.^2 + t3(iknotg).*psibar.^3;
rhot = sqrt(T/T(nr));
% The derivative of T is q, evaluate for psibar values
qpsi = sign(rzero*bzero)/(psimag-psibry)*...
  (t1(iknotg) + 2*t2(iknotg).*psibar + 3*t3(iknotg).*psibar.^2);
if calculate_profile_responses
  dt0dx = d0*dstdx;
  dt1dx = d1*dstdx;
  dt2dx = d2*dstdx;
  dt3dx = d3*dstdx;
  dTdx = dt0dx(iknotg,:)              + dt1dx(iknotg,:).*psibarnx + ...
         dt2dx(iknotg,:).*psibarnx.^2 + dt3dx(iknotg,:).*psibarnx.^3;
  drhotdx = 1/2/T(nr)./rhot*ones(1,nx).*(dTdx-T/T(nr)*dTdx(nr,:));
  drhotdx(1,:) = 0;
  drhotdx(nr,:) = 0;
end



if exist('eta_vs_rhot') & ~isempty(eta_vs_rhot)
% Calculate area-integral of resistive voltage within flux surfaces

Tzr = t0(iknot) + ...
  t1(iknot).*psibarzr + ...
  t2(iknot).*psibarzr.^2 + ...
  t3(iknot).*psibarzr.^3;
rhotzr = sqrt(Tzr/T(nr));

% Calculate eta versus psibar
rhoe = linspace(0,1,length(eta_vs_rhot))';
eta = spline(rhoe, eta_vs_rhot, rhot);
etazr = nan(nz,nr);
etazr(iplasma) = spline(rhoe, eta_vs_rhot, rhotzr(iplasma));
dumzr(1:np) = 2*pi*rgg(iplasma).*etazr(iplasma).*pcurrt(iplasma);
Bgg = zeros(ngg,1);
Bgg(1:np) = cumsum(dumzr(kgg(1:np)));
Btot = Bgg(np);
pprimea = p1(nkn)*twopi/(psibry-psimag);
ffprima = f1(nkn)*twopi/(psibry-psimag);
pprimeb = (p1(nkn) + 2*p2(nkn) + 3*p3(nkn))*twopi/(psibry-psimag);
ffprimb = (f1(nkn) + 2*f2(nkn) + 3*f3(nkn))*twopi/(psibry-psimag);
for j = 1:ngg
  if ncell(j) > 0
    if iplasma(j)
      dum = (RAcell(j)-RA(j))*pprime(nr)+(ARcell(j)-AR(j))*ffprim(nr)/mu0;
    else
      dum = RAcell(j)*pprime(nr)+ARcell(j)*ffprim(nr)/mu0;
    end
    Btot = Btot + 2*pi*rgg(j)*eta(nr)*dum;
  end
end

if calculate_profile_responses
  detadrhot = (spline(rhoe, eta_vs_rhot, rhot+1e-6) - ...
               spline(rhoe, eta_vs_rhot, rhot-1e-6))/2e-6;
  detadx = detadrhot*ones(1,nx).*drhotdx;
  detadpsibar = (spline(psibar, eta, psibar+1e-6) - ...
                 spline(psibar, eta, psibar-1e-6))/2e-6;

  detadrhotzr = spline(psibar, detadrhot, psibarzr);
  detadpsibarzr = spline(psibar, detadpsibar, psibarzr);
  dTzrdpsibarzr = t1(iknot) + 2*t2(iknot).*psibarzr + 3*t3(iknot).*psibarzr.^2;
  dTzrdx = dTzrdpsibarzr(:)*ones(1,nx).*dpsibarzrdx + ...
                                      dt0dx(iknot,:) + ...
           psibarzr(:)   *ones(1,nx).*dt1dx(iknot,:) + ...
           psibarzr(:).^2*ones(1,nx).*dt2dx(iknot,:) + ...
	   psibarzr(:).^3*ones(1,nx).*dt3dx(iknot,:);
  drhotzrdx = 1/2/T(nr)./rhotzr(:)*ones(1,nx).*(dTzrdx-Tzr(:)/T(nr)*dTdx(nr,:));
  detazrdrhotzr = (...
    spline(rhoe, eta_vs_rhot, rhotzr+1e-6)-...
    spline(rhoe, eta_vs_rhot, rhotzr-1e-6))/2e-6;
  detazrdx = detazrdrhotzr(:)*ones(1,nx).*drhotzrdx;
  dBggdx = zeros(ngg,nx);
  dBggdx(1:np,:) = ...
    2*pi*rgg(iplasma).*etazr(iplasma)*ones(1,nx).*dpcurrtdx(iplasma,:) + ...
    2*pi*rgg(iplasma).*pcurrt(iplasma)*ones(1,nx).*detazrdx(iplasma,:);
  dBggdx(1:np,:) = cumsum(dBggdx(kgg(1:np),:));
  dBtotdx = dBggdx(np,:);
  dpprimeadx = -pprimea/(psibry-psimag)*(dpsibrydx-dpsimagdx);
  dpprimeadx(1,indsp) = dpprimeadx(1,indsp) + c1(nkn,:)*twopi/(psibry-psimag);
  dffprimadx = -ffprima/(psibry-psimag)*(dpsibrydx-dpsimagdx);
  dffprimadx(1,indsf) = dffprimadx(1,indsf) + c1(nkn,:)*twopi/(psibry-psimag);
  dpprimebdx = -pprimeb/(psibry-psimag)*(dpsibrydx-dpsimagdx);
  dpprimebdx(1,indsp) = dpprimebdx(1,indsp) + ...
    (c1(nkn,:) + 2*c2(nkn,:) + 3*c3(nkn,:))*twopi/(psibry-psimag);
  dffprimbdx = -ffprimb/(psibry-psimag)*(dpsibrydx-dpsimagdx);
  dffprimbdx(1,indsf) = dffprimbdx(1,indsf) + ...
    (c1(nkn,:) + 2*c2(nkn,:) + 3*c3(nkn,:))*twopi/(psibry-psimag);
  for j = 1:ngg
    if ncell(j) > 0
      if iplasma(j)
	dumxr = dRAcelldx(j,:)*pprime(nr) + ...
	         dARcelldx(j,:)*ffprim(nr)/mu0 + ...
		 (RAcell(j)-RA(j))*dpprimebdx + ...
		 (ARcell(j)-AR(j))/mu0*dffprimbdx;
      else
	dumxr = dRAcelldx(j,:)*pprime(nr) + ...
	         dARcelldx(j,:)*ffprim(nr)/mu0 + ...
	         RAcell(j)*dpprimebdx + ...
		 ARcell(j)/mu0*dffprimbdx;
      end
      dBtotdx = dBtotdx + 2*pi*rgg(j)*eta(nr)*dumxr;
    end
  end
end
if profile_fit_option == 0
  % Simple fit
  sb = M(1:np,1:nkn+2)\Bgg(1:np);
  if calculate_profile_responses
    dsbdx = zeros(nkn+2,nx);
    for j = 1:nx
      ps3(1:np,:) = dyggdx(1:np,j)*ones(1,nkn+2).*ps2(1:np,:);
      MTM12 = M(1:np,1:nkn+2)'*ps3(1:np,1:nkn+2)+ps3(1:np,1:nkn+2)'*M(1:np,1:nkn+2);
      dsbdx(:,j) = MTMi12*(ps3(1:np,:)'*Bgg(1:np) - MTM12*sb + ...
        M(1:np,1:nkn+2)'*dBggdx(1:np,j));
    end
  end
elseif profile_fit_option == 1
  % Fit that requires total area-integral of resistive voltage to be Btot
  Cgg(1:np) = Bgg(1:np)-M(1:np,nkn+2)*Btot;
  sb = [M(1:np,1:nkn+1)\Cgg(1:np); Btot];
  if calculate_profile_responses
    dsbdx = [zeros(nkn+1,nx); dBtotdx];
    for j = 1:nx
      ps3(1:np,:) = dyggdx(1:np,j)*ones(1,nkn+2).*ps2(1:np,:);
      MTM11 = M(1:np,1:nkn+1)'*ps3(1:np,1:nkn+1)+ps3(1:np,1:nkn+1)'*M(1:np,1:nkn+1);
      dsbdx(1:nkn+1,j) = MTMi11*(ps3(1:np,1:nkn+1)'*Cgg(1:np) - ...
        MTM11*sb(1:nkn+1) + M(1:np,1:nkn+1)'*dBggdx(1:np,j) - ...
	M(1:np,1:nkn+1)'*(ps3(1:np,nkn+2)*Btot + M(1:np,nkn+2)*dBtotdx(j)));
    end
  end
else
  % Fit that also requires derivative on axis to be consistent with eta(1)
  % The ellipse has b/a = sqrt(psimagrr/psimagzz)
  % Area-integral of Vres near axis = pi*a*b*eta(1)*j_maxis, 
  Bpa = 2*pi*eta(1)*Apa*(rmaxis^2*pprimea+ffprima/mu0);  
  Cgg(1:np) = Bgg(1:np)-M(1:np,1)*Bpa-M(1:np,nkn+2)*Btot;
  sb = [Bpa; M(1:np,2:nkn+1)\Cgg(1:np); Btot];
  if calculate_profile_responses
    dBpadx = 2*pi*eta(1)*dApadx*(rmaxis^2*pprimea+ffprima/mu0) + ...
      4*pi*eta(1)*Apa*rmaxis*drmaxisdx*pprimea + ...
      2*pi*eta(1)*Apa*(rmaxis^2*dpprimeadx+dffprimadx/mu0);
    dsbdx = [dBpadx; zeros(nkn,nx); dBtotdx];
    for j = 1:nx
      ps3(1:np,:) = dyggdx(1:np,j)*ones(1,nkn+2).*ps2(1:np,:);
      MTM21 = M(1:np,2:nkn+1)'*ps3(1:np,2:nkn+1)+ps3(1:np,2:nkn+1)'*M(1:np,2:nkn+1);
      dsbdx(2:nkn+1,j) = MTMi21*(ps3(1:np,2:nkn+1)'*Cgg(1:np) - ...
        MTM21*sb(2:nkn+1) + M(1:np,2:nkn+1)'*dBggdx(1:np,j) - ...
	M(1:np,2:nkn+1)'*(ps3(1:np,nkn+2)*Btot + M(1:np,nkn+2)*dBtotdx(j) + ...
	ps3(1:np,1)*Bpa + M(1:np,1)*dBpadx(j)));
    end
  end
end
% area-integral of resistive voltage function
b0 = d0*sb;
b1 = d1*sb;
b2 = d2*sb;
b3 = d3*sb;
B = b0(iknotg) + b1(iknotg).*psibar + b2(iknotg).*psibar.^2 + b3(iknotg).*psibar.^3;
Bprime = 2*pi/(psibry-psimag)*(b1(iknotg) + 2*b2(iknotg).*psibar + 3*b3(iknotg).*psibar.^2);
Vres = Bprime./Aprime;
if calculate_profile_responses
  db0dx = d0*dsbdx;
  db1dx = d1*dsbdx;
  db2dx = d2*dsbdx;
  db3dx = d3*dsbdx;
  dBdx = db0dx(iknotg,:)              + db1dx(iknotg,:).*psibarnx + ...
         db2dx(iknotg,:).*psibarnx.^2 + db3dx(iknotg,:).*psibarnx.^3;
  dBprimedx = 2*pi/(psibry-psimag)*(...
      db1dx(iknotg,:)           + ...
    2*db2dx(iknotg,:).*psibarnx + ...
    3*db3dx(iknotg,:).*psibarnx.^2) - ...
    Bprime/(psibry-psimag)*(dpsibrydx-dpsimagdx);
  dVresdx = dBprimedx./(Aprime*ones(1,nx))-Vres./Aprime*ones(1,nx).*dAprimedx;
end

% Calculate area-integral of poloidal flux within flux surfaces

dumzr(1:np) = psizr(iplasma);
Hgg(1:np) = Ag*cumsum(dumzr(kgg(1:np)));
Htot = Hgg(np);
for j = 1:ngg
  if ncell(j) > 0
    if iplasma(j)
      dum = (Acell(j)-Ag)*psibry;
    else
      dum = Acell(j)*psibry;
    end
    Htot = Htot + dum;
  end
end

if calculate_profile_responses
  dHggdx = zeros(ngg,nx);
  dHggdx(1:np,:) = dpsizrdx(iplasma,:);
  dHggdx(1:np,:) = Ag*cumsum(dHggdx(kgg(1:np),:));
  dHtotdx = dHggdx(np,:);
  for j = 1:ngg
    if ncell(j) > 0
      if iplasma(j)
	dumxr = (Acell(j)-Ag)*dpsibrydx + dAcelldx(j,:)*psibry;
      else
	dumxr = dAcelldx(j,:)*psibry + Acell(j)*dpsibrydx;
      end
      dHtotdx = dHtotdx + dumxr;
    end
  end
end

if profile_fit_option == 0
  % Simple fit
  sh = M(1:np,1:nkn+2)\Hgg(1:np);
  if calculate_profile_responses
    dshdx = zeros(nkn+2,nx);
    for j = 1:nx
      ps3(1:np,:) = dyggdx(1:np,j)*ones(1,nkn+2).*ps2(1:np,:);
      MTM12 = M(1:np,1:nkn+2)'*ps3(1:np,1:nkn+2)+ps3(1:np,1:nkn+2)'*M(1:np,1:nkn+2);
      dshdx(:,j) = MTMi12*(ps3(1:np,:)'*Hgg(1:np) - MTM12*sh + ...
        M(1:np,1:nkn+2)'*dHggdx(1:np,j));
    end
  end
elseif profile_fit_option == 1
  % Fit that requires total enclosed poloidal flux to be Htot
  Cgg(1:np) = Hgg(1:np)-M(1:np,nkn+2)*Htot;
  sh = [M(1:np,1:nkn+1)\Cgg(1:np); Htot];
  if calculate_profile_responses
    dshdx = [zeros(nkn+1,nx); dHtotdx];
    for j = 1:nx
      ps3(1:np,:) = dyggdx(1:np,j)*ones(1,nkn+2).*ps2(1:np,:);
      MTM11 = M(1:np,1:nkn+1)'*ps3(1:np,1:nkn+1)+ps3(1:np,1:nkn+1)'*M(1:np,1:nkn+1);
      dshdx(1:nkn+1,j) = MTMi11*(ps3(1:np,1:nkn+1)'*Cgg(1:np) - ...
        MTM11*sh(1:nkn+1) + M(1:np,1:nkn+1)'*dHggdx(1:np,j) - ...
	M(1:np,1:nkn+1)'*(ps3(1:np,nkn+2)*Htot + M(1:np,nkn+2)*dHtotdx(j)));
    end
  end
else
  % Fit that also requires derivative on axis to be psimag*Apa
  Hpa = psimag*Apa;  
  Cgg(1:np) = Hgg(1:np)-M(1:np,1)*Hpa-M(1:np,nkn+2)*Htot;
  sh = [Hpa; M(1:np,2:nkn+1)\Cgg(1:np); Htot];
  if calculate_profile_responses
    dfpol1dx = zeros(1,nx);
    dfpol1dx(:,indsf) = dfpol1dx(:,indsf)+c0(1,:)/fpol1;  
    dHpadx = psimag*dApadx+Apa*dpsimagdx;
    dshdx = [dHpadx; zeros(nkn,nx); dHtotdx];
    for j = 1:nx
      ps3(1:np,:) = dyggdx(1:np,j)*ones(1,nkn+2).*ps2(1:np,:);
      MTM21 = M(1:np,2:nkn+1)'*ps3(1:np,2:nkn+1)+ps3(1:np,2:nkn+1)'*M(1:np,2:nkn+1);
      dshdx(2:nkn+1,j) = MTMi21*(ps3(1:np,2:nkn+1)'*Cgg(1:np) - ...
        MTM21*sh(2:nkn+1) + M(1:np,2:nkn+1)'*dHggdx(1:np,j) - ...
	M(1:np,2:nkn+1)'*(ps3(1:np,nkn+2)*Htot + M(1:np,nkn+2)*dHtotdx(j) + ...
	ps3(1:np,1)*Hpa + M(1:np,1)*dHpadx(j)));
    end
  end
end
% (area-integral of poloidal flux) function
h0 = d0*sh;
h1 = d1*sh;
h2 = d2*sh;
h3 = d3*sh;
H = h0(iknotg) + h1(iknotg).*psibar + h2(iknotg).*psibar.^2 + h3(iknotg).*psibar.^3;
if calculate_profile_responses
  dh0dx = d0*dshdx;
  dh1dx = d1*dshdx;
  dh2dx = d2*dshdx;
  dh3dx = d3*dshdx;
  dHdx = dh0dx(iknotg,:)              + dh1dx(iknotg,:).*psibarnx + ...
         dh2dx(iknotg,:).*psibarnx.^2 + dh3dx(iknotg,:).*psibarnx.^3;
end
end

