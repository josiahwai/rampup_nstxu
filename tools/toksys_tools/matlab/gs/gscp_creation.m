%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gscp_creation
%
%  PURPOSE: Create initial circular plasma with aspect ratio = 50
%
%  INPUTS: run gs_breakdown to create:
%          rbdef, zbdef, boundary-defining point
%          lim, index into rl, zl where bdef is in range lim:lim+1
%
%  OUTPUTS: Initialized equilibrium

%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
%  METHOD: 
	
%  NOTES:  
	
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander ON 2015-03-03
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%gscp_configure

% Initial a0 is smallest that converges well with circular model
if exist('smallest_aminor','var') & length(smallest_aminor) == 1
  a0 = smallest_aminor;
else
  a0 = rbdef/50;
end

% lim is index into rl, zl with bdef within lim:lim+1
vlim = [turnin*[drl(lim);dzl(lim)]]'/dl(lim);
r0 = rbdef + vlim(1)*a0;
z0 = zbdef + vlim(2)*a0;

% Some extra numbers for the boundary
rmin = r0 - a0;
rmax = r0 + a0;
zmin = z0 - a0;
zmax = z0 + a0;

% The moving coordinates where flux is measured
rh = linspace(rmin,rmax,nrh);
drh = (rmax-rmin)/(nrh-1);

% Create initial flux at rh
n1 = ceil(nrh/2)+1;
n2 = nrh - n1;
psih = [linspace(1,0,n1) linspace(0,1,n2)].^2;
sf = sf0;
sp = sp0;
plasma = 1;
gscp_analysis_response

% The change of r, z, a that is needed for vertical force balance
drza = cpMi(np-2:np,np-1)*perr(np-1);

% The factor for psihpla to make the plasma radius = a0
f = -(psihapp(end)-psihapp(1))/(psihpla(end)-psihpla(1));

% Flux for the created plasma
psih = psihapp + f*psihpla;

% Current for the created plasma
cpasma1 = f*cpasma;

% Calculate what the current is with sp, sf
gscp_analysis_response

% Scale sp, sf to get cpasma1
sp = cpasma1/cpasma*sp;
sf = cpasma1/cpasma*sf;

% Run once more
gscp_analysis_response

if exist('open_field_line_volume','var') & ~isempty(open_field_line_volume)
  if isfield(index_in_y,'Vtot')
    lae.y(index_in_y.Vtot) = open_field_line_volume;
  end
end
if exist('open_field_line_rmaxis','var') & ~isempty(open_field_line_rmaxis)
  if isfield(index_in_y,'rmaxis')
    lae.y(index_in_y.rmaxis) = open_field_line_rmaxis;
  end
end

disp(['Switch to GS evolution at Ip > ' num2str(round(cpasma1)) ...
  ' A (for the vertical field at plasma initiation)'])

return
plot(rh,psiherr-psiherr(1),rh,psih-psih(1))

% THIS CODE BELOW NEEDS SAME LOGIC AS IN gscp_converge
% MAKE IT AN OPTION TO CONVERGE WITHOUT CHANGING a0

% Construct an equation for converging *without changing a0*
% The variables in this M version are: total flux change at rh points,
% shifts of r0, z0, and a factor * [sp;sf], that replaces a0
% The equations are: changes in (total - plasma = applied flux), and:
% 1) flux difference psih(end)-psih(1), 
% 2) psitopapp-psibotapp
% 3) displacement into the wall
% The right-hand-side to 1 & 2 can be non-zero when fixing errors
dpsihappdp(:,end) = 0; % Replacing a0 with a factor f * [sp;sf]
dpsihpladp(:,end) = psihpla(:);
dpsibotappdp(end) = 0;
dpsitopappdp(end) = 0;
dnbdefdp(end) = 0;
cpM = ...
  [eye(nrh,np)-dpsihappdp-dpsihpladp; ... % tot- pla = app flux
  -1 zeros(1,nrh-2) 1 0 0 0 ; ...   % psih(end) - psih(1)
  dpsitopappdp-dpsibotappdp; ... % psitopapp-psibotapp
  dnbdefdp];                        % = 0 to not go through wall
perr = [psiherr(:);psih(end)-psih(1);psitopapp-psibotapp;0];
cpMi = inv(cpM);
dp = cpMi*perr;
psih(:) = psih(:) - dp(1:nrh);
r0 = r0 - dp(nrh+1);
z0 = z0 - dp(nrh+2);
%sp = (1-dp(np))*sp;
%sf = (1-dp(np))*sf;
gscp_analysis_response
plot(rh,dp(1:nrh)-dp(1),rh,psih-psih(1))

cpasma

% The internal states depend on the flag constraints
states.sp = sp;
states.sf = sf;
states.cpasma = cpasma;
states.li = li;
states.betap = betap;
states.Wth = Wth;





return

cpM2 = ...
  [eye(nrh,np-1)-dpsihappdp(:,1:end-1)-dpsihpladp(:,1:end-1); ...
  dpsitopappdp(1:end-1)-dpsibotappdp(1:end-1); ...
  dnbdefdp(1:end-1)];
p2err = [psiherr(:);psitopapp-psibotapp;0];
cpM2i = inv(cpM2);
dp2 = cpM2i*p2err;
psih(:) = psih(:) - dp2(1:nrh);
r0 = r0 - dp2(nrh+1);
z0 = z0 - dp2(nrh+2);
gscp_analysis_response
dp = cpMi*perr;
plot(rh,psiherr-psiherr(1),rh,psih-psih(1),...
  rh,dp(1:nrh)-dp(1))
