%  USAGE:   gs_profile_response
%
%  PURPOSE: Calculate response of profile functions to x
%
%  INPUTS: output from gs_contour_profiles
%
%  OUTPUTS: dVcdx, d(volume within contours)/dx
%           dAcdx, d(area within contours)/dx
%           dLcdx, d(area-integral of 1/R within contours)/dx
%	
%  METHOD: 
	
%  NOTES:  
	
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander  ON	5/13/14
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~plasma
  return
end


dpprimecdpsibarc = ...
  (2*p2(iknotc) + 6*p3(iknotc).*psibarc(:))*twopi/(psibry-psimag);
dpprimecdpsimag = pprimec/(psibry-psimag);
dpprimecdpsibry = dpprimecdpsimag;

dhalffpolsquaredcdpsibarc = ...
  f1(iknotc) + 2*f2(iknotc).*psibarc(:)+3*f3(iknotc).*psibarc(:).^2;

dfpolcdpsibarc = sign(rzero*bzero)./sqrt(2*halffpolsquaredc).*...
  dhalffpolsquaredcdpsibarc;

dffprimcdpsibarc = ...
  (2*f2(iknotc) + 6*f3(iknotc).*psibarc(:))*twopi/(psibry-psimag);
dffprimcdpsimag = ffprimc/(psibry-psimag);
dffprimcdpsibry = dffprimcdpsimag;

dfprimcdpsibarc = dffprimcdpsibarc./fpolc-fprimc./fpolc.*dfpolcdpsibarc;


dprescdx = zeros(ncont,nx);
dprescdx(:,indsp) = c0(iknotc,:) + ...
  c1(iknotc,:).*(psibarc(:)*ones(1,nkn+2)) + ...
  c2(iknotc,:).*(psibarc(:).^2*ones(1,nkn+2)) + ...
  c3(iknotc,:).*(psibarc(:).^3*ones(1,nkn+2));
dprescdx(psibarc==1,:) = 0;

dpprimecdx = -pprimec/(psibry-psimag)*(dpsibrydx-dpsimagdx);
dpprimecdx(:,indsp) = dpprimecdx(:,indsp) + (...
  c1(iknotc,:) + ...
  c2(iknotc,:)*2.*(psibarc(:)*ones(1,nkn+2)) + ...
  c3(iknotc,:)*3.*(psibarc(:).^2*ones(1,nkn+2)))*...
  twopi/(psibry-psimag);

dhalffpolsquaredcdsf = c0(iknotc,:) + ...
  c1(iknotc,:).*(psibarc(:)*ones(1,nkn+2)) + ...
  c2(iknotc,:).*(psibarc(:).^2*ones(1,nkn+2)) +...
  c3(iknotc,:).*(psibarc(:).^3*ones(1,nkn+2));

dfpolcdx = zeros(ncont,nx);
dfpolcdx(:,indsf) = sign(rzero*bzero)./sqrt(2*halffpolsquaredc)*...
  ones(1,nkn+2).*dhalffpolsquaredcdsf;

dffprimcdx = -ffprimc/(psibry-psimag)*(dpsibrydx-dpsimagdx);
dffprimcdx(:,indsf) = dffprimcdx(:,indsf) + (...
  c1(iknotc,:) + ...
  c2(iknotc,:)*2.*(psibarc(:)*ones(1,nkn+2)) + ...
  c3(iknotc,:)*3.*(psibarc(:).^2*ones(1,nkn+2)))*...
  twopi/(psibry-psimag);
  
dfprimcdx = dffprimcdx./(fpolc*ones(1,nx)) - ...
  fprimc./fpolc*ones(1,nx).*dfpolcdx;

dVcdx = zeros(ncont,nx);
dAcdx = zeros(ncont,nx);
dLcdx = zeros(ncont,nx);
for i = 1:ncont

  i1 = (i-1)*npola+1;
  i2 = i*npola;
  
  r_x = drcontdx(i1+1:i2,:)+drcontdx(i1:i2-1,:);
  z_x = dzcontdx(i1+1:i2,:)-dzcontdx(i1:i2-1,:);
  rnx = (rcont(i1+1:i2)+rcont(i1:i2-1))'*ones(1,nx);
  znx = (zcont(i1+1:i2)-zcont(i1:i2-1))'*ones(1,nx);
  lnr = log(rcont(i1+1:i2)+rcont(i1:i2-1))'*ones(1,nx);
  
  dVcdx(i,:) = pi/4*(2*sum(r_x.*rnx.*znx) + sum(rnx.*rnx.*z_x));
  dAcdx(i,:) = (sum(r_x.*znx) + sum(rnx.*z_x))/2;
  dLcdx(i,:) = (sum(r_x./rnx.*znx) + sum(lnr.*z_x));

end

dTcdx = 0.50*[zeros(1,nx); ...
  cumsum((dfpolcdx(2:end,:)+dfpolcdx(1:end-1,:)).*(diff(Lc)*ones(1,nx)) + ...
  (fpolc(2:end)+fpolc(1:end-1))*ones(1,nx).*diff(dLcdx))];

if isfield(index_in_y,'T')
  dydx(index_in_y.T,:) = dTcdx;
end

