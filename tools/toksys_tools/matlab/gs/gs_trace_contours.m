%  USAGE:   gs_trace_contours
%
%  PURPOSE: Trace contours of constant flux in the plasma
%           This version iterates to find the contours
%
%  INPUTS: psibarzr, psibarzr = (psizr-psimag)/(psibry-psimag);
%          psibarc, normalized fluxes at contours (default psibar)
%          npola, number of poloidal angles (default 2*nr-1)
%          rbdef, zbdef, position that defines the boundary
%          anglesoption, for poloidal angles of contour points:
%            1. start at bdef point, equally spaced angles (default)
%            2. Like 1 for boundary but along gradient toward axis
%          rmaxis, zmaxis, position of axis
%          rg, zg, dr, dz, nr, nz (grid variables)
%          For cubic interpolation on the grid:
%            mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2
%            neighbors = reshape((-1:2)'*ones(1,4)+ones(4,1)*(-1:2)*nz,1,16)
%
%  OUTPUTS: rcont, zcont, ncont, = R, Z, number of contours
%	
%  METHOD: Newton-Rhapson method used on all contour-points in parallel
	
%  NOTES:  To-do: deduce rhocont min, max after each iteration
%          to safe-guard against points shooting off?
%          Right now a speed limit of dr/3 per iteration is imposed
	
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander  ON	4/21/13
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~plasma
  return
end

if ~exist('psibarzr','var')
  psibarzr = (psizr-psimag)/(psibry-psimag);
end

if ~exist('psibarc','var')
  psibarc = psibar;
end
ncont = length(psibarc);

if ~exist('npola','var') | isempty(npola)
  npola = 2*nr-1;
end

if ~exist('anglesoption','var')
  anglesoption = 1;
end

rcont = zeros(npola,ncont);
zcont = zeros(npola,ncont);
coscont = zeros(npola,ncont);
sincont = zeros(npola,ncont);
psibarcont_target = ones(npola,1)*psibarc(:)';

thab = angle((rbdef-rmaxis)+1i*(zbdef-zmaxis));
dths = linspace(0,2*pi,npola);

n = 0; % Number of turns around bbbs
k = 1; % points to a bbbs point
th_previous = thbbbs(nbbbs-1)-2*pi; % previous theta
rh_previous = rhobbbs(nbbbs-1); % previous rho
th = thbbbs(k)+n*2*pi;
rh = rhobbbs(k);
for j = 1:npola
  costh = cos(thab+dths(j));
  sinth = sin(thab+dths(j));
  while thab+dths(j) > th
    if k == nbbbs
      k = 1;
      n = n+1;
    else
      k = k+1;
    end
    th_previous = th; % previous theta
    rh_previous = rh; % previous rho
    th = thbbbs(k)+n*2*pi;
    rh = rhobbbs(k);
  end
  dth = th - th_previous;
  x = (th-thab-dths(j))/dth;
  rho = x*rh_previous + (1-x)*rh;
  rcont(j,:) = rmaxis+rho*costh*psibarc;
  zcont(j,:) = zmaxis+rho*sinth*psibarc;
  coscont(j,:) = costh;
  sincont(j,:) = sinth;
end
rcont([1 npola],psibarc==1) = rbdef;
zcont([1 npola],psibarc==1) = zbdef;
rcont(:,psibarc==0) = rmaxis;
zcont(:,psibarc==0) = zmaxis;

contour_iterations = 0;
drhomax = dr;

while drhomax > dr/1e9 & contour_iterations < 50

contour_iterations = contour_iterations+1;

kcont = floor((rcont-rg(1))/dr)*nz+floor((zcont-zg(1))/dz)+1;

xcont = (rcont-rgg(kcont))/dr;
wrc = [ones(npola*ncont,1) xcont(:) xcont(:).^2 xcont(:).^3]*mx;
wrc_r = [zeros(npola*ncont,1) ones(npola*ncont,1) 2*xcont(:) 3*xcont(:).^2]*mx/dr;

xcont = (zcont-zgg(kcont))/dz;
wzc = [ones(npola*ncont,1) xcont(:) xcont(:).^2 xcont(:).^3]*mx;
wzc_z = [zeros(npola*ncont,1) ones(npola*ncont,1) 2*xcont(:) 3*xcont(:).^2]*mx/dz;

icont = [kcont(:)-(nz+1)   kcont(:)-nz     kcont(:)-(nz-1)   kcont(:)-(nz-2) ...
         kcont(:)-1        kcont(:)        kcont(:)+1        kcont(:)+2      ...
	 kcont(:)+(nz-1)   kcont(:)+nz     kcont(:)+(nz+1)   kcont(:)+(nz+2) ...
	 kcont(:)+(2*nz-1) kcont(:)+(2*nz) kcont(:)+(2*nz+1) kcont(:)+(2*nz+2)];

wcont = [wzc(:,1).*wrc(:,1) wzc(:,2).*wrc(:,1) wzc(:,3).*wrc(:,1) wzc(:,4).*wrc(:,1) ...
	 wzc(:,1).*wrc(:,2) wzc(:,2).*wrc(:,2) wzc(:,3).*wrc(:,2) wzc(:,4).*wrc(:,2) ...
	 wzc(:,1).*wrc(:,3) wzc(:,2).*wrc(:,3) wzc(:,3).*wrc(:,3) wzc(:,4).*wrc(:,3) ...
	 wzc(:,1).*wrc(:,4) wzc(:,2).*wrc(:,4) wzc(:,3).*wrc(:,4) wzc(:,4).*wrc(:,4)];

wcont_r = [wzc(:,1).*wrc_r(:,1) wzc(:,2).*wrc_r(:,1) wzc(:,3).*wrc_r(:,1) wzc(:,4).*wrc_r(:,1) ...
	   wzc(:,1).*wrc_r(:,2) wzc(:,2).*wrc_r(:,2) wzc(:,3).*wrc_r(:,2) wzc(:,4).*wrc_r(:,2) ...
	   wzc(:,1).*wrc_r(:,3) wzc(:,2).*wrc_r(:,3) wzc(:,3).*wrc_r(:,3) wzc(:,4).*wrc_r(:,3) ...
	   wzc(:,1).*wrc_r(:,4) wzc(:,2).*wrc_r(:,4) wzc(:,3).*wrc_r(:,4) wzc(:,4).*wrc_r(:,4)];

wcont_z = [wzc_z(:,1).*wrc(:,1) wzc_z(:,2).*wrc(:,1) wzc_z(:,3).*wrc(:,1) wzc_z(:,4).*wrc(:,1) ...
	   wzc_z(:,1).*wrc(:,2) wzc_z(:,2).*wrc(:,2) wzc_z(:,3).*wrc(:,2) wzc_z(:,4).*wrc(:,2) ...
	   wzc_z(:,1).*wrc(:,3) wzc_z(:,2).*wrc(:,3) wzc_z(:,3).*wrc(:,3) wzc_z(:,4).*wrc(:,3) ...
	   wzc_z(:,1).*wrc(:,4) wzc_z(:,2).*wrc(:,4) wzc_z(:,3).*wrc(:,4) wzc_z(:,4).*wrc(:,4)];

psibarcont = reshape(sum(wcont'.*psibarzr(icont)')',npola,ncont);
psibarcont_r = reshape(sum(wcont_r'.*psibarzr(icont)')',npola,ncont);
psibarcont_z = reshape(sum(wcont_z'.*psibarzr(icont)')',npola,ncont);
psibarcont_rho = coscont.*psibarcont_r + sincont.*psibarcont_z;

drhocont = (psibarcont_target-psibarcont) ./ psibarcont_rho;
%drhocont([1 npola],psibarc==1) = 0;
%drhocont(:,psibarc==0) = 0;
drhocont(drhocont > dr/3) = dr/3;
drhocont(drhocont < -dr/3) = -dr/3;
rcont = rcont+drhocont.*coscont;
zcont = zcont+drhocont.*sincont;

drhomax = max(abs(drhocont(psibarcont_target>0)));

end

for j = 1:ncont
  if psibarc(j) == 0
    rcont(:,j) = rmaxis;
    zcont(:,j) = zmaxis;
  elseif psibarc(j) == 1
    rcont(1,j) = rbdef;
    zcont(1,j) = zbdef;
    rcont(npola,j) = rbdef;
    zcont(npola,j) = zbdef;
  end
end

xcont = (rcont-rgg(kcont))/dr;
wrc_rr = [zeros(npola*ncont,2) 2*ones(npola*ncont,1) 6*xcont(:)]*mx/dr^2;
xcont = (zcont-zgg(kcont))/dz;
wzc_zz = [zeros(npola*ncont,2) 2*ones(npola*ncont,1) 6*xcont(:)]*mx/dz^2;

wcont_rr = [wzc(:,1).*wrc_rr(:,1) wzc(:,2).*wrc_rr(:,1) ...
            wzc(:,3).*wrc_rr(:,1) wzc(:,4).*wrc_rr(:,1) ...
	    wzc(:,1).*wrc_rr(:,2) wzc(:,2).*wrc_rr(:,2) ...
	    wzc(:,3).*wrc_rr(:,2) wzc(:,4).*wrc_rr(:,2) ...
	    wzc(:,1).*wrc_rr(:,3) wzc(:,2).*wrc_rr(:,3) ...
	    wzc(:,3).*wrc_rr(:,3) wzc(:,4).*wrc_rr(:,3) ...
	    wzc(:,1).*wrc_rr(:,4) wzc(:,2).*wrc_rr(:,4) ...
	    wzc(:,3).*wrc_rr(:,4) wzc(:,4).*wrc_rr(:,4)];

wcont_rz = [wzc_z(:,1).*wrc_r(:,1) wzc_z(:,2).*wrc_r(:,1) ...
            wzc_z(:,3).*wrc_r(:,1) wzc_z(:,4).*wrc_r(:,1) ...
	    wzc_z(:,1).*wrc_r(:,2) wzc_z(:,2).*wrc_r(:,2) ...
	    wzc_z(:,3).*wrc_r(:,2) wzc_z(:,4).*wrc_r(:,2) ...
	    wzc_z(:,1).*wrc_r(:,3) wzc_z(:,2).*wrc_r(:,3) ...
	    wzc_z(:,3).*wrc_r(:,3) wzc_z(:,4).*wrc_r(:,3) ...
	    wzc_z(:,1).*wrc_r(:,4) wzc_z(:,2).*wrc_r(:,4) ...
	    wzc_z(:,3).*wrc_r(:,4) wzc_z(:,4).*wrc_r(:,4)];

wcont_zz = [wzc_zz(:,1).*wrc(:,1) wzc_zz(:,2).*wrc(:,1) ...
            wzc_zz(:,3).*wrc(:,1) wzc_zz(:,4).*wrc(:,1) ...
	    wzc_zz(:,1).*wrc(:,2) wzc_zz(:,2).*wrc(:,2) ...
	    wzc_zz(:,3).*wrc(:,2) wzc_zz(:,4).*wrc(:,2) ...
	    wzc_zz(:,1).*wrc(:,3) wzc_zz(:,2).*wrc(:,3) ...
	    wzc_zz(:,3).*wrc(:,3) wzc_zz(:,4).*wrc(:,3) ...
	    wzc_zz(:,1).*wrc(:,4) wzc_zz(:,2).*wrc(:,4) ...
	    wzc_zz(:,3).*wrc(:,4) wzc_zz(:,4).*wrc(:,4)];
psibarcont_rr = reshape(sum(wcont_rr'.*psibarzr(icont)')',npola,ncont);
psibarcont_rz = reshape(sum(wcont_rz'.*psibarzr(icont)')',npola,ncont);
psibarcont_zz = reshape(sum(wcont_zz'.*psibarzr(icont)')',npola,ncont);

lae.rcont = rcont;
lae.zcont = zcont;
