%  USAGE:   gs_contour_profiles
%
%  PURPOSE: Calculate profile functions using contours rcont, zcont
%
%  INPUTS: output from gs_trace_contours
%
%  OUTPUTS: Vc, volume within contours
%           Ac, area within contours
%           Lc, area-integral of 1/R within contours
%           All contour-derived profiles of size [ncont,1] have suffix c
%	
%  METHOD: 
	
%  NOTES:  
	
%
%  WRITTEN BY:  Anders Welander  ON	5/12/14
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lcont = sum(sqrt(diff(rcont).^2+diff(zcont).^2))';

if ~plasma
  return
end

iknotc = ones(ncont,1);
for j = 2:nkn
  iknotc(psibarc > psikn(j)) = j;
end

% Coefficients for the nkn third degree polynomials for pres
p0 = c0*sp;
p1 = c1*sp;
p2 = c2*sp;
p3 = c3*sp;

% Coefficients for the nkn third degree polynomials for fpol^2/2-fpol(nr)^2/2
f0 = c0*sf;
f1 = c1*sf;
f2 = c2*sf;
f3 = c3*sf;

presc = p0(iknotc) + ...
        p1(iknotc).*psibarc(:) + ...
        p2(iknotc).*psibarc(:).^2 + ...
        p3(iknotc).*psibarc(:).^3;
presc(psibarc==1) = 0;

pprimec = (p1(iknotc) + ...
           p2(iknotc)*2.*psibarc(:) + ...
	   p3(iknotc)*3.*psibarc(:).^2)*twopi/(psibry-psimag);

halffpolsquaredc = f0(iknotc) + ...
                   f1(iknotc).*psibarc(:) + ...
		   f2(iknotc).*psibarc(:).^2+...
                   f3(iknotc).*psibarc(:).^3 + ...
		   rzero^2*bzero^2/2;

halffpolsquaredc(halffpolsquaredc < 0) = 0;

fpolc = sign(rzero*bzero)*sqrt(2*halffpolsquaredc);

ffprimc = (f1(iknotc) + ...
           f2(iknotc)*2.*psibarc(:) + ...
	   f3(iknotc)*3.*psibarc(:).^2)*twopi/(psibry-psimag);
fprimc = ffprimc./fpolc;

rbrcont = -1/2/pi*reshape(sum(wcont_z'.*psizr(icont)')',npola,ncont);
rbzcont = +1/2/pi*reshape(sum(wcont_r'.*psizr(icont)')',npola,ncont);
rbpcont = sqrt(rbrcont.^2+rbzcont.^2);
rbtcont = ones(npola,1)*fpolc';
bpcont = rbpcont./rcont;
btcont = rbtcont./rcont;
brcont = rbrcont./rcont;
bzcont = rbzcont./rcont;

drcont = diff(rcont);
dzcont = diff(zcont);
dpolc = sqrt(drcont.^2+dzcont.^2);

qc = 1/2/pi*sign(rzero*bzero)*sum(dpolc.*...
 (btcont(1:npola-1,:)+btcont(2:npola,:))./...
 (rbpcont(1:npola-1,:)+rbpcont(2:npola,:)))';

for i = 1:ncont
  if psibarc(i) == 0
    emaxis = sqrt(abs(warr*psizr(iia)/(wazz*psizr(iia))));
    j_axis = rmaxis*pprimec(i)+ffprimc(i)/rmaxis/mu0;
    qc(i) = abs(fpolc(i)/j_axis)/mu0 * ...
            (1+emaxis*emaxis)/(emaxis*rmaxis^2);
  end
end
lae.psibarc = psibarc;
lae.qc = qc;

% primes, i.e. derivatives w.r.t. (poloidal flux per radian)
% Suffix for prime is "p" for quantities evaluated at psibarc
% Suffix for prime is "prime" for quantities evaluated at psibar

Vc = pi/4*sum((rcont(2:npola,:)+rcont(1:npola-1,:)).^2.*...
         (zcont(2:npola,:)-zcont(1:npola-1,:)))';

Ac = sum((rcont(2:npola,:)+rcont(1:npola-1,:)).*...
         (zcont(2:npola,:)-zcont(1:npola-1,:)))'/2;

Lc = sum(log(rcont(2:npola,:)+rcont(1:npola-1,:)).*...
	    (zcont(2:npola,:)-zcont(1:npola-1,:)))';

Acp = pi/(psibry-psimag)*(...
      sum((coscont(2:npola  ,:)./psibarcont_rho(2:npola  ,:) + ...
           coscont(1:npola-1,:)./psibarcont_rho(1:npola-1,:)).*...
          (zcont(2:npola,:)-zcont(1:npola-1,:)))' + ...
      sum((rcont(2:npola,:)+rcont(1:npola-1,:)).*... 	   
          (sincont(2:npola  ,:)./psibarcont_rho(2:npola  ,:) - ...
	   sincont(1:npola-1,:)./psibarcont_rho(1:npola-1,:)))');

Tc = 0.50*[0; cumsum((fpolc(2:end)+fpolc(1:end-1)).*diff(Lc))];
Wc = 0.75*[0; cumsum((presc(2:end)+presc(1:end-1)).*diff(Vc))];
Ic = [0; cumsum((pprimec(2:end)+pprimec(1:end-1)).*diff(Vc))/(4*pi)] + ...
     [0; cumsum((ffprimc(2:end)+ffprimc(1:end-1)).*diff(Lc))/(2*mu0)];

rhot = sqrt(Tc/Tc(end));

BP2FLX = (mu0*Ic./lcont).^2;
Bc = 4/3*mu0*Wc./Vc./BP2FLX;

if isfield(index_in_y,'T')
  lae.y(index_in_y.T) = Tc;
end

if isfield(index_in_y,'qpsi')
  lae.y(index_in_y.qpsi) = qc;
end

lae.Tc = Tc;
lae.rhot = rhot;
