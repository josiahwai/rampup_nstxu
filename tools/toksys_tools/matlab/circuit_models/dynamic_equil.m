function [moving_equil,dbg_out,ier] = dynamic_equil(inputs)
 %
%  SYNTAX: moving_equil = dynamic_equil(inputs)
%
%  PURPOSE: Calculate "moving equilibrium" so that delta coordinates can
%   be defined. Get objects to fit model
%   where we take d/dt i_PF as input disturbance or forcing function.
%
%  INPUT: (in structure "inputs")
%    shotnum = shot number
%    ip0   = approximate plasma current (Amps)
%    vmat = mapping of input voltages to rows in state equation, after reduction
%    vdata = approximate equil. power supply voltages (V) corresponding to cols in vmat
%    good_vdata = vector of flags indicating which vdata values are "good" (0 or 1)
%    tok_system = system model data structure for device being modeled
%    Rp = plasma resistance
%    ohstates_idx = indices of PF current states that participate in Ip control
%   Optional (first 3 must be included or omitted together):
%    LHSconstraint = defines extra LHS * (current states) = RHS constraints
%    RHSconstraint = defines extra LHS * (current states) = RHS constraints
%    constraint_str = string description of these constraints
%    weight_coils = scalar or length ncx vector wtg factor for fitting coil currents
%			 (default=100)
%    weight_volts = scalar or vector wtg factor (length(vdata)) for fitting power
%			 supply voltages (default=10)
%
%  OUTPUT: (in structure moving_equil)
%       equil_pt   = nominal equilibrium to be used for linearization:
% 
%  RESTRICTIONS: If a particular voltage measurement is of poor or unknown quality, set
%	the good_vdata flag = 0  and define a low fitting weight. 
%
%  METHOD:  Performs a weighted least squares fit of free variables (see text
%	variable "var" in output structure moving_equil) to shot data and
%	constraint equations containing those variables defined by model 
%	circuit equations and physical constraints.  "Bad" voltages are replaced by 0
%	during the fitting, but plots show comparison of input and fitted values.
 
%  WRITTEN BY:  Mike Walker     ON      7/1/09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)dynamic_equil.m	1.6 04/29/11

%    coil_states = approximate equilibrium PF-coil current states (Amps):

show_equil_plots = 1;

moving_equil=[]; dbg_out=[];

struct_to_ws(inputs);
struct_to_ws(tok_system);
coil_states = tok_system.Ieq(1:ncx);

if ~exist('weight_coils')
   weight_coils = 100*ones(1,length(coil_states));;
else
   if(length(weight_coils)~=length(coil_states))
      wait('ERROR dynamic_equil: coil_states and weight_coils must be same length')
      return;
   end
end
if ~exist('weight_volts')
   weight_volts = 10*ones(1,length(vdata));
else
   if(length(weight_volts)~=length(vdata))
      wait('ERROR dynamic_equil: vdata and weight_volts must be same length')
      return;
   end
end
if(exist('LHSconstraint')~=exist('RHSconstraint'))
   wait('ERROR dynamic_equil: one of LHS and RHS constraints not specified')
   return;
end
if(exist('LHSconstraint')~=1)
   LHSconstraint = [];
   RHSconstraint = [];
   constraint_strs = '';
else
   if(size(LHSconstraint,1) ~= length(RHSconstraint))
      wait('ERROR dynamic_equil: LHSconstraint and RHSconstraint are incompatible sizes')
      return;
   end
   if(size(LHSconstraint,1) ~= size(constraint_str,1))
      wait('ERROR dynamic_equil: LHSconstraint and constraint_str are incompatible sizes')
      return;
   end
   constraint_strs = constraint_str;
end
if(size(vmat,1) ~= size(tok_system.bmat,1))
   wait('ERROR dynamic_equil: size of vmat inconsistent with system model')
   return;
end

% blockR = resistance matrix for model, after state reduction
% blockM = mutual inductance matrix for model, after state reduction
if(isfield(tok_system,'netlist_model'))
   blockR = netlist_model.Rhat;
   blockM = netlist_model.Mhat;
else
   blockR = Pxx' * (rxx + Rextra) * Pxx;
   blockM = Pxx' * (mxx + Lextra) * Pxx;
end

% Variable to solve for is:
%   x = [IOHdot0, IPF0, Iv0, v0, Ip0]

if(size(coil_states,1)==1), coil_states = coil_states'; end;
if(size(vdata,1)==1), vdata = vdata'; end;
lenv = length(vdata);
ncx = tok_system.ncx;
nvx = tok_system.nvx; %    number of vessel segment states
noh = length(ohstates_idx);
nvar = noh + ncx + nvx + lenv + (ip0>0);
v_states = ncx+1:ncx+nvx; % vessel current indices in overall set of states

if(length(coil_states)~=ncx)
   wait('ERROR dynamic_equil: length of input coil_states must match ncx in system model');
   return;
end

var = ['IOHdot0(' int2str(noh) '),IPFstates0(' int2str(ncx) ...
	'), Iv0(' int2str(nvx) '),v0(' int2str(lenv) ')'];
if(ip0 > 0)
   var = [var ',Ip0'];
end

LHS=[]; RHS=[];
constraint_index = [];
if(~exist('constraint_strs'))
   constraint_strs = [];
end

[mM,nM]=size(blockM);

%------------------------------------------------------------------
%		CONSTRAINT EQUATIONS
%------------------------------------------------------------------

% Impose specified added constraints:
if(~isempty(LHSconstraint))
   [LHS,RHS,ier] = update_constraint_eqns(LHS,RHS,nvar,LHSconstraint,RHSconstraint);
   if(ier)
      wait(['ERROR dynamic_equil: update_constraint_eqns returned ier = ' int2str(ier)]);
      return;
   end
   [m,n]=size(LHS);

   constraint_index = [constraint_index; size(LHS,1)];
   fprintf('constraints %d to %d are %s\n',1, ...
	constraint_index(end),constraint_strs(end,:));
else
   constraint_index = 0;
end
%------------------------------------------------------------------

% Ic0 should be approximately equal to measured coil currents at equilibrium:
clear addLHS addRHS
Iwt = diag(weight_coils);
addLHS(1:ncx,noh+[1:ncx]) = Iwt;
addRHS = [Iwt*coil_states];

[LHS,RHS,ier] = update_constraint_eqns(LHS,RHS,nvar,addLHS,addRHS);
if(ier)
   wait(['ERROR dynamic_equil: update_constraint_eqns returned ier = ' int2str(ier)]);
   return;
end
[m,n]=size(LHS);

constraint_index = [constraint_index; size(LHS,1)];
constraint_strs = strvcat(constraint_strs, 'Ic0=Ic,eq');
fprintf('constraints %d to %d are %s\n',constraint_index(end-1)+1, ...
	constraint_index(end),constraint_strs(end,:));
%------------------------------------------------------------------

if(ip0 > 0)
% plasma current Ip0 should be approximately Ip,eq
   clear addLHS addRHS
   addLHS = 1*[zeros(1,nvar-1) 1];
   addRHS = 1*[ip0];

   [LHS,RHS,ier] = update_constraint_eqns(LHS,RHS,nvar,addLHS,addRHS);
   if(ier)
      wait(['ERROR dynamic_equil: update_constraint_eqns returned ier = ' int2str(ier)]);
      return;
   end

   constraint_index = [constraint_index; size(LHS,1)];
   constraint_strs = strvcat(constraint_strs, 'Ip0=Ip,eq');
   fprintf('constraints %d to %d are %s\n',constraint_index(end-1)+1, ...
	constraint_index(end),constraint_strs(end,:));
end
%------------------------------------------------------------------

% vessel current Iv0 should satisfy linearized equation:
clear addLHS addRHS
MvIdot = blockM(v_states,ohstates_idx);
Rv  = blockR(v_states,v_states);
addLHS = 0.1*[MvIdot zeros(nvx,ncx) Rv zeros(nvx,lenv)];
if(ip0 > 0)
   addLHS = [addLHS  zeros(nvx,1)];
end
addRHS = zeros(nvx,1);

[LHS,RHS,ier] = update_constraint_eqns(LHS,RHS,nvar,addLHS,addRHS);
if(ier)
   wait(['ERROR dynamic_equil: update_constraint_eqns returned ier = ' int2str(ier)]);
   return;
end

constraint_index = [constraint_index; size(LHS,1)];
constraint_strs = strvcat(constraint_strs, 'Iv0 equation');
fprintf('constraints %d to %d are %s\n',constraint_index(end-1)+1, ...
	constraint_index(end),constraint_strs(end,:));
%------------------------------------------------------------------

% 0th order terms in conductor circuit equation:
clear addLHS addRHS
addLHS = 10000*[blockM(1:ncx,ohstates_idx),blockR(1:ncx,1:ncx),...
					zeros(ncx,nvx),-vmat(1:ncx,:)];
if(ip0>0)
   addLHS = [addLHS zeros(ncx,1)];
end
addRHS = [zeros(ncx,1)];

[LHS,RHS,ier] = update_constraint_eqns(LHS,RHS,nvar,addLHS,addRHS);
if(ier)
   wait(['ERROR dynamic_equil: update_constraint_eqns returned ier = ' int2str(ier)]);
   return;
end
[m,n]=size(LHS);

constraint_index = [constraint_index; size(LHS,1)];
constraint_strs = strvcat(constraint_strs, 'conductor circuit equation');
fprintf('constraints %d to %d are %s\n',constraint_index(end-1)+1, ...
	constraint_index(end),constraint_strs(end,:));
%------------------------------------------------------------------

if(ip0 > 0) 	% must satisfy 0th order plasma circuit equation:
   clear addLHS addRHS
   addLHS(1,1:noh) = blockM(mM,ohstates_idx);
   addLHS(1,nvar) = Rp;
   addRHS(1) = 0;
   addLHS = 10000*addLHS;
   addRHS = 10000*addRHS;
   %mM=mM-1;
   %nM=nM-1;

   [LHS,RHS,ier] = update_constraint_eqns(LHS,RHS,nvar,addLHS,addRHS);
   if(ier)
      wait(['ERROR dynamic_equil: update_constraint_eqns returned ier = ' int2str(ier)]);
      return;
   end

   [m,n]=size(LHS);

   constraint_index = [constraint_index; size(LHS,1)];
   constraint_strs = strvcat(constraint_strs, 'plasma circuit equation');
   fprintf('constraints %d to %d are %s\n',constraint_index(end-1)+1, ...
	constraint_index(end),constraint_strs(end,:));
end
%------------------------------------------------------------------

% voltages should be approximately equal to data
clear addLHS addRHS
temp = vdata;
for k=1:length(vdata)
   if(good_vdata(k))
      temp(k) = vdata(k);
   else
      temp(k) = 0;
   end
end
Vwt = diag(weight_volts);
addLHS = [zeros(lenv,noh+ncx+nvx) Vwt*eye(lenv)];
if(ip0 > 0)
   addLHS = [addLHS zeros(lenv,1)];
end
addRHS = Vwt*[temp];

[LHS,RHS,ier] = update_constraint_eqns(LHS,RHS,nvar,addLHS,addRHS);
if(ier)
   wait(['ERROR dynamic_equil: update_constraint_eqns returned ier = ' int2str(ier)]);
   return;
end

constraint_index = [constraint_index; size(LHS,1)];
constraint_strs = strvcat(constraint_strs, 'Vc0=Vc,eq');
fprintf('constraints %d to %d are %s\n',constraint_index(end-1)+1, ...
	constraint_index(end),constraint_strs(end,:));

equil_pt = LHS\RHS;		

fit_error = LHS*equil_pt-RHS;
fprintf('normalized error = %f\n',norm(fit_error)/norm(RHS));

descriptions = struct( ...
'var','list of variables in dynamic equilibrium (equil_pt)', ...
'equil_pt','nominal equilibrium to be used for linearization', ...
'constraint_index','number of equations in each group in constraint_strs', ...
'constraint_strs','descriptions of groups of constraint equations');

% If needed, remove non-existent index for "added constraints".
if(constraint_index(1)==0)
   constraint_index = constraint_index(2:end);
end

% variable = [IOHdot0, IPF0, Iv0, v0, Ip0]
leni = length(ohstates_idx);
m_Ioh = equil_pt(1:leni);				% current slope(s) of OH coil(s)
IPF_states = equil_pt(leni+1:leni+ncx);			% nominal PF currents
IPF0 = Pcc*IPF_states;
Iv0 = equil_pt(leni+ncx+1:leni+ncx+nvx);		% nominal vessel currents
V0 = equil_pt(leni+ncx+nvx+1:leni+ncx+nvx+lenv);	% nominal coil voltages
Ip0 = equil_pt(nvar);				% nominal plasma current

moving_equil = struct( ...
'var',var, ...
'equil_pt', equil_pt, ...
'm_Ioh', m_Ioh, ... 
'IPF0', IPF0, ...
'Iv0', Iv0, ... 
'V0', V0, ... 
'Ip0', Ip0, ... 
'fit_error', fit_error, ...
'constraint_index',constraint_index, ...
'constraint_strs',constraint_strs, ...
'LHS', LHS, ...
'RHS', RHS, ...
'descriptions',descriptions);

dbg_out = struct( ...
'blockM',blockM, ...
'blockR',blockR, ...
'LHS', LHS, ...
'RHS', RHS);

%save test_new
%wait('stop')

if show_equil_plots
   figure_num=100;

   next_figure
   plot(coil_states,'x-')
   hold on
   plot(IPF_states,'ro-')
   hold off
   xlabel('coil number')
   ylabel('current (A)')
   title('dynamic equilibrium currents(o) vs. measured(x)')

   next_figure
   plot(vdata,'x-')
   hold on
   plot(V0,'ro-')
   plot(temp,'c*-')
   hold off
   xlabel('input voltage number')
   ylabel('voltage (V)')
   title('dynamic equilibrium voltages(o) vs. data(x), constraints(*)')

if 0
   next_figure
   tsec = equil_time*1e-3;
   [d,t]=getptd(shotnum,'ecoila',tsec-.5,tsec+.5);
   plot(t,d)
   grid on
   hold on
   [mm,mi]=min(abs(t-tsec));
   dd = d(mi);
   plot([t(1);tsec;t(length(t))],[dd-0.5*m_Ioh(1);dd;dd+0.5*m_Ioh(1)],'r');
   hold off
   title('dynamic equilibrium ramp vs. measured (ecoila)')
   ylabel('current (A)')
   xlabel('time (sec)')
end

   next_figure
   plot(fit_error)
   title('dynamic equilibrium fit errors')
   xlabel('constraint number')
   ylabel('error')
end
