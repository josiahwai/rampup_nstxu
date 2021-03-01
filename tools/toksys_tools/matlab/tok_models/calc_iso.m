function [isoresp,iso] = calc_iso(rziso,tok_data_struct,resp,eq)
%
%  USAGE: [isoresp,iso] = calc_iso(rziso,tok_data_struct,resp,eq)
%
%  PURPOSE: Calculate flux, field and responses thereof at isoflux points
%
%  INPUTS: rziso: R, Z of isoflux points arranged as [R(:) Z(:)] (unit: meters)
%          tok_data_struct: toksys description of tokamak
%          resp (optional): the output from a response model (gspert or rzrig)
%          eq (optional): either equilibrium with field psizr or psizr itself
%
%  OUTPUTS: isoresp: response structure with fields:
%                    dpsidis: flux response to conductors [Wb/A]
%                    dpsidip,*dli,*dbetap, etc: exogenous responses
%                    dbrdis, dbzdis: Br, Bz responses to conductors [T/A]
%                    dbrdip,dbzdip, etc: exogenous Br, Bz responses
%               iso: equilibrium structure with fields:
%                    psi: flux at rziso points
%                    br, bz: Br, Bz at rziso points
%
%  RESTRICTIONS: rziso must be within the grid

%
%  METHOD: 
%	
%  VERSION @(#)calc_iso.m	1.4 08/30/11
%
%  WRITTEN BY:  Anders Welander  ON	5/17/11
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  struct_to_ws(tok_data_struct); dr = mean(diff(rg)); dz = mean(diff(zg));
  nss = nc+nv;
  ngg = nr*nz;
  % Find weights in grid points to calculate value at a point using cubic Hermite spline (ref. wikipedia)
  mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2;
  neighbors = reshape([-1-nz -1 -1+nz -1+2*nz;-nz 0 nz 2*nz;1-nz 1 1+nz 1+2*nz;2-nz 2 2+nz 2+2*nz],1,16);

  % Retrieve psizr to calculate fluxes and fields at rziso for the equilibrium
  if exist('eq','var') & ~isempty(eq)
    if isstruct(eq)
      psizr = eq.psizr;
    else
      psizr = eq;
    end
    % Now calculate fluxes and fields
    for j = 1:size(rziso,1)
      kr0 = min(nr-3,max(1,floor((rziso(j,1)-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
      kz1 = min(nz-2,max(2,ceil((rziso(j,2)-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
      k = kr0*nz+kz1;
      ii = k+neighbors; % ii indexes 16 points around rziso(j,:)
      tr = (rziso(j,1)-rgg(k))/dr; tz = (rziso(j,2)-zgg(k))/dz;
      w = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      wr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
      wz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      iso.psi(j,1) = psizr(ii)*w';
      iso.br(j,1) = -psizr(ii)*wz'/rziso(j,1);
      iso.bz(j,1) = +psizr(ii)*wr'/rziso(j,1);
    end
  end
  % Always make this field so that iso exists even if eq was not supplied in call
  iso.descriptions.flux = 'flux at rziso points for supplied equilibrium [Wb]';
  iso.descriptions.br = 'Radial field at rziso points for supplied equilibrium [T]';
  iso.descriptions.bz = 'Vertical field at rziso points for supplied equilibrium [T]';
  
  if nargin<2
    isoresp = [];
    help calc_iso
    disp('**********************')
    disp('Supply more arguments!')
    disp('**********************')
    return
  end
  
  if ~exist('resp','var') | isempty(resp)
    resp.dcphidis = zeros(ngg,nss);
  end
  
  % Parse the resp structure
  if isfield(resp,'dcphidis') % true for gspert
  elseif isfield(resp,'drdis') & isfield(resp,'dcdr') % true for rzrig
    resp.dcphidis = resp.dcdr(:)*resp.drdis+resp.dcdz(:)*resp.dzdis;
    resp.dcphidip = resp.dcdr(:)*resp.drdip+resp.dcdz(:)*resp.dzdip;
    resp.dcphidli = resp.dcdr(:)*resp.drdli;
    resp.dcphidbetap = resp.dcdr(:)*resp.drdbetap+resp.dcdz(:)*resp.dzdbetap;
  else
    wait('error calc_iso: can not interpret structure resp')
    return
  end
  % Since exogenous variables may vary, use field names to extract all that apply
  f = fieldnames(resp); % f holds names of fields in resp
  ind_resp_in_f = []; % will become indices in f containing dcphid*
  for j = 1:length(f)
    if strfind(char(f(j)),'dcphid') == 1
      ind_resp_in_f(end+1) = j;
    end
  end
  
  % Now calculate responses using dcphid* for plasma response and tok_data_struct for mutuals
  for j = 1:size(rziso,1)
    kr0 = min(nr-3,max(1,floor((rziso(j,1)-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
    kz1 = min(nz-2,max(2,ceil((rziso(j,2)-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
    k = kr0*nz+kz1;
    ii = k+neighbors; % ii indexes 16 points around rziso(j,:)
    m = zeros(16,nz*nr); % will hold mutuals to the 16 points around an rziso point
    for jj = 1:16
      z_idx = mod(ii(jj)-1,nz)+1; % z index of output element ii(jj)
      r_idx = ceil(ii(jj)/nz);    % r index of output element ii(jj)
      for jr=1:nr
        k1 = (r_idx-1)*nz;
	k2 = (jr-1)*nz;
	for jz=1:nz
	  m(jj,k2+jz) = mpp(k1+abs(z_idx-jz)+1,jr);
	end
      end
    end
    % We now have dpsizr(ii)/dis = m*resp.dcphidis+[mpc(ii,:) mpv(ii,:)], etc
    tr = (rziso(j,1)-rgg(k))/dr; tz = (rziso(j,2)-zgg(k))/dz;
    w = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
    wr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
    wz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
    % We now have dpsi(rziso)/dis = w*(m*dpsizr(ii))/dis, etc
    for jj = 1:length(ind_resp_in_f)
      respobj = char(f(ind_resp_in_f(jj)));
      if strmatch(respobj,'dcphidis')
        isoresp.dpsidis(j,:) = w*(m*resp.dcphidis + [mpc(ii,:) mpv(ii,:)]);
        isoresp.dbrdis(j,:) = -wz*(m*resp.dcphidis + [mpc(ii,:) mpv(ii,:)])/rziso(j,1)/2/pi;
        isoresp.dbzdis(j,:) = +wr*(m*resp.dcphidis + [mpc(ii,:) mpv(ii,:)])/rziso(j,1)/2/pi;
      elseif eval(['prod(size(resp.' respobj '))']) == nr*nz % respobj may be on form (nz,nr) so collapse
        eval(['isoresp.dpsid' respobj(7:end) '(j,:) = w*(m*resp.' respobj '(:));'])
        eval(['isoresp.dbrd' respobj(7:end) '(j,:) = -wz*(m*resp.' respobj '(:))/rziso(j,1)/2/pi;'])
        eval(['isoresp.dbzd' respobj(7:end) '(j,:) = +wr*(m*resp.' respobj '(:))/rziso(j,1)/2/pi;'])
      else
        eval(['isoresp.dpsid' respobj(7:end) '(j,:) = w*(m*resp.' respobj ');'])
        eval(['isoresp.dbrd' respobj(7:end) '(j,:) = -wz*(m*resp.' respobj ')/rziso(j,1)/2/pi;'])
        eval(['isoresp.dbzd' respobj(7:end) '(j,:) = +wr*(m*resp.' respobj ')/rziso(j,1)/2/pi;'])
      end
    end
  end
  for jj = 1:length(ind_resp_in_f)
    respobj = char(f(ind_resp_in_f(jj)));
    eval(['isoresp.descriptions.dpsid' respobj(7:end) ' = ''flux response to ' respobj(7:end) ''';'])
    eval(['isoresp.descriptions.dbrd' respobj(7:end) ' = ''Br response to ' respobj(7:end) ''';'])
    eval(['isoresp.descriptions.dbzd' respobj(7:end) ' = ''Bz response to ' respobj(7:end) ''';'])
  end
  isoresp.descriptions.dpsidis = 'flux response to is (conductor currents) [Wb/A]';
  isoresp.descriptions.dbrdis = 'Br response to is (conductor currents) [T/A]';
  isoresp.descriptions.dbzdis = 'Bz response to is (conductor currents) [T/A]';
  
  
