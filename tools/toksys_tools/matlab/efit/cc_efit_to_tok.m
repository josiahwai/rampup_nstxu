function equil_I = cc_efit_to_tok(vac_objs,equil_data)
 %
%  SYNTAX: equil_I = cc_efit_to_tok(vac_objs,equil_data)
%
%  PURPOSE: Convert equilibrium conductor currents in internal EFIT 
%	representation (contained in equil_data) into a representation 
%	compatible with toksys environment.
%
%  INPUT:
%    vac_objs = vacuum model data object structure
%    equil_data = equilibrium data structure
%    idx_efit_to_tok = optional map of efit pf coil indices to toksys indices, s.t. if
%	                    I=currents in efit order, then
%`                     I(idx_efit_to_tok)=currents in toksys order
%    idxvv_efit_to_tok = " vv indices " " "
%
%  OUTPUT: (units defined by imks, iterminal in vac_objs)
%    equil_I = structure containing:
%             cc0t = coil current equilibrium vector, toksys conductors
%             cc0e = coil current equilibrium vector, EFIT conductors
%             vc0t = vessel current equilibrium vector, toksys conductors
%             vc0e = vessel current equilibrium vector, EFIT conductors
%             cprojet = mapping: cc0t = cprojet*cc0e
%             cprojte = mapping: cc0e = cprojte*cc0t
%             vprojte = mapping: vc0e = vprojte*vc0t
%    (Equilibrium vessel currents are usually all zero. In those cases, vc0e
%	and vprojte are not produced.)
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	8/31/09 (from in-line code in rzrig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Derived values:
   if vac_objs.imks, iscale=1e6; else iscale=1; end
   iterminal = vac_objs.iterminal;
   ncc = size(vac_objs.mcc,1);

% Extract Needed EFIT data from equil_data:

   if(~isfield(vac_objs,'ecnturn') | isempty(vac_objs.ecnturn))
      ecnturn = 0;
   else
      ecnturn = vac_objs.ecnturn;
   end
   if(ecnturn~=0)
      ccnturn = [ecnturn; vac_objs.fcnturn];
      if(~iterminal)
         ccnturn = ones(size(ccnturn));
      end
   else
      ccnturn = vac_objs.fcnturn;
      if(~iterminal)
         ccnturn = ones(size(ccnturn));
      end
   end
% make all "turn" objects into columns, consistent with coil current vectors
   [s1,s2]=size(ccnturn); if(s1==1), ccnturn = ccnturn'; end

% Define desired currents for calculating equilibrium field.
% Assumes cc came out of read_gfile always in MA-turns.

if(~isfield(equil_data,'cc') | isempty(equil_data.cc))

   wait('ERROR cc_efit_to_tok: no valid cc data - returning empty matrices')
   equil_I = struct( ...
   'cc0t', [], ...
   'cc0e', [], ...
   'cprojte', []);

else
   ccx = iscale*equil_data.cc(:);      % ccx now in proper units (MA or Amps)

% Complicated logic required here because either turnfc or fcturn can hold
% the fraction of current values - depends on how used internally in EFIT.

   if(any(equil_data.fcturn<1))
      Ifrac = equil_data.fcturn;
      nturn = equil_data.turnfc;
   elseif(any(equil_data.turnfc<1))
      Ifrac = equil_data.turnfc;
      nturn = equil_data.fcturn;
   elseif(any(equil_data.turnfc>1))
      Ifrac = equil_data.fcturn;
      nturn = equil_data.turnfc;
   elseif(any(equil_data.fcturn>1))
      Ifrac = equil_data.turnfc;
      nturn = equil_data.fcturn;
   else                                 % all must be = 1
      Ifrac = equil_data.turnfc;
      nturn = equil_data.fcturn;
   end
% if any ecoils, insert before fcoils:

   idxf = abs(equil_data.fcid);	% use abs in case anti-series connection
   isign = sign(equil_data.fcid);
   if(~isempty(equil_data.ecid))
      nee = length(unique(equil_data.ecid));
      Ifrac = [ones(1,nee) Ifrac];
      for i=nee:-1:1
         idx = find(equil_data.ecid==i);
         neturn = sum(equil_data.ecturn(idx));
         nturn = [neturn nturn];
      end
      idx = [[1:nee] idxf+nee];
      isign = [ones(1,nee) isign];
   else
      idx = idxf;
   end

   [mm,i1] = min(size(Ifrac)); [mm,i2] = min(size(ccx));
   if(i1~=i2) Ifrac = Ifrac'; end;
   cc0 = isign'.*ccx(idx).*Ifrac;

   nn = length(cc0);
   if(isfield(equil_data,'idx_efit_to_tok'))
      idx_efit_to_tok = equil_data.idx_efit_to_tok;
   else
      idx_efit_to_tok = 1:length(nturn);
   end

   if(iterminal)
      [mm,i1] = min(size(nturn)); [mm,i2] = min(size(ccnturn));
      if(i1~=i2) nturn = nturn'; end;
% (length(ccnturn) > length(nturn) is allowed. Extra currents are 0.)
      if(norm(ccnturn(1:length(nturn)) - nturn(idx_efit_to_tok)))
        fprintf('WARNING cc_efit_to_tok:\n');
        fprintf('  #turns equil. data=');
        fprintf('%d ',nturn); fprintf('\n');
        fprintf('  #turns vacuum objects=');
        fprintf('%d ',ccnturn(1:length(nturn))); fprintf('\n');
%        wait
      end
      cc0 = cc0./nturn;
   end

% The following accounts for coils assumed to carry no steady state current
% NOTE that these are always assumed to be listed last.
   if(length(cc0)<ncc)
      cc0 = [cc0; zeros(ncc-length(cc0),1)];
   end

% At this point, cc0 is in correct units, lumped/terminal mode, and represents
% equilibrium currents assigned to all coils represented by vacuum objects.

cprojet = zeros(length(cc0),length(equil_data.cc));
cprojte = zeros(length(equil_data.cc),length(cc0));
if(isfield(vac_objs,'ecnturn'))
   nec = length(vac_objs.ecnturn);
   for k=1:nec
      cprojet(k,k)=1;        % ohmic coil
      cprojte(k,k)=1;        % ohmic coil
   end
else
   nec = 0;
end
for k=nec+1:size(cprojte,1)
   idx = find(equil_data.fcid==k-nec);
   if(~isempty(idx))
      cprojet(idx+nec,k)=1;
      cprojte(k,idx(1)+nec)=1;
   end
end

if(isfield(equil_data,'idx_efit_to_tok'))
   nn = length(cc0);
   ni = max(idx_efit_to_tok);
   temp = zeros(nn);
   for k=1:ni
      temp(idx_efit_to_tok(k),k)=1;
   end
   cc0t = zeros(size(cc0));
   cc0t(1:ni) = cc0(idx_efit_to_tok);
   cprojte = cprojte*temp;
   cprojet = temp'*cprojet;
else
   cc0t = cc0;
end

cc0e = cprojte*cc0t;

equil_I = struct( ...
'cc0t', cc0t, ...
'cc0e', cc0e, ...
'cprojet', cprojet, ...
'cprojte', cprojte);

description = struct( ...
'cc0t', 'coil currents corresponding to toksys vacuum objects', ...
'cc0e', 'coil currents corresponding to EFIT fitted coils', ...
'cprojet', 'mapping: cc0t = cprojet*cc0e', ...
'cprojte', 'mapping: cc0e = cprojte*cc0t');

end

% If vessel currents in equilibrium, these must be used in rzrig calculation.
if(isfield(equil_data,'vc') & ~isempty(equil_data.vc))
   vcx = iscale*equil_data.vc;       % vcx now in proper units (MA or Amps)
   vprojte= proj_turn(equil_data.vvid,equil_data.vvfrac);
   vc0= vprojte'*vcx;  
   

%  construct a projection without inverting
   vprojte= proj_turn(equil_data.vvid,1./equil_data.vvfrac);
   ii = find(vprojte~=0);
   vprojtej= zeros(size(vprojte));
   vprojtej(ii)= 1;
   num= sum(vprojtej')';
   vprojtej= (1./num)*ones(1,size(vprojte,2));
   vprojte= vprojte.*vprojtej;
   
   % Translate conductors if toksys & efit have different indices or #
   if(isfield(equil_data,'idxvv_efit_to_tok'))
      nn = length(vc0);
      ni = max(equil_data.idxvv_efit_to_tok);
      temp = zeros(nn,ni);
      for k=1:ni
         temp(equil_data.idxvv_efit_to_tok(k),k)=1;
      end
      vc0t = zeros(vac_objs.nv,1);
      vc0t(1:ni) = vc0(equil_data.idxvv_efit_to_tok);
      vprojte = vprojte*temp;
   else
      vc0t = vc0;
   end
   
   vc0e = vprojte*vc0t;
   equil_I.vc0t = vc0t(1:vac_objs.nv);
   equil_I.vc0e = vc0e;
   equil_I.vprojte = vprojte(1:vac_objs.nv);

   description.vc0t='vessel currents corresponding to toksys vacuum objects from equil_data';
   description.cc0e='vessel currents corresponding to toksys.def_connect fitted conductors';
   description.vprojte= 'mapping: vc0e = vprojte*vc0t';
else
   equil_I.vc0t = zeros(vac_objs.nv,1);
   description.vc0t='vessel currents corresponding to toksys vacuum objects';
end

equil_I.description = description;
