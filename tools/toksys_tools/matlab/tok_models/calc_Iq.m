function [Iqresp,Iq,eqx] = calc_Iq(q,tok_data_struct,eq,options)
%
%  USAGE: [Iqresp,Iq,eqx] = calc_Iq(q,tok_data_struct,eq,options)
%
%  PURPOSE: Calculate amount of current within q-surfaces and responses thereof
%
%  INPUTS: q: q-values for surfaces
%          tok_data_struct: toksys description of tokamak
%          eq: equilibrium structure
%          options (optional), do help gspert for description
%
%  OUTPUTS: Iqresp: response structure with fields:
%                   dIqdis: response of current inside q to conductor currents
%                   dIqdip,*dli,*dbphi, etc: exogenous responses
%               Iq: current within q-contours
%              eqx: extra equilibrium information, help gspert for description

%
%  RESTRICTIONS: 

%
%  METHOD: 
%	
%  VERSION @(#)calc_Iq.m	1.1 05/18/11
%
%  WRITTEN BY:  Anders Welander  ON	5/17/11
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  
  if nargin<3
    Iqresp = [];
    help calc_Iq
    disp('**********************')
    disp('Supply more arguments!')
    disp('**********************')
    return
  end
  
  struct_to_ws(tok_data_struct);
  
  % Call gspert to get all information about currents with flux contours
  options.q = q;
  if ~isfield(options,'iconstraints')
    options.iconstraints = 1;
  end
  [response,eqx] = gspert(eq,tok_data_struct,options);
  struct_to_ws(eqx);
  psigrid = linspace(eq.psimag,eq.psibry,nr);
  
  % I is current with surfaces of flux psigrid, psiq is flux at specified q
  Iq = spline(psigrid,I,psiq)';
  
  % dIdpsi is how (current within flux surface) varies with the flux of the surface
  dIdpsi = (spline(psigrid,I,psigrid+1e-6)-spline(psigrid,I,psigrid-1e-6))/2e-6;
  dIqdpsi = spline(psigrid,dIdpsi,psiq);
  
  v = ones(1,size(response.dqpsidis,2));
  % Below the first term is how current within surfaces of flux psiq responds
  % The second term is due to the change of the flux at q-surfaces
  Iqresp.dIqdis = spline(psigrid,response.dIdis',psiq)'-2*pi*(dIqdpsi./qprime)'*v.*spline(psigrid,response.dqpsidis',psiq)';
  
  % Now we calculate the response to exogenous variables
  f = fieldnames(response); % f holds names of fields in response
  for j = 1:length(f)
    respobj = char(f(j));
    if strfind(respobj,'dcphid') == 1
      if strmatch(respobj,'dcphidis')
      else
        x = respobj(7:end);
        v = ones(1,size(getfield(response,['dqpsid' x]),2));
        eval(['Iqresp.dIqd' x ' = spline(psigrid,response.dId' x ''',psiq)''- ' ...
          '2*pi*(dIqdpsi./qprime)''*v.*spline(psigrid,response.dqpsid' x ''',psiq)'';'])
      end
    end
  end
  
  Iqresp.dIqdbphi = -2*pi*(dIqdpsi./qprime)'.*spline(psigrid,response.dqpsidbphi',psiq)';
