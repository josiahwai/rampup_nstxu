   function [resvv,vals,vecs] = calc_vv_resvals(tok_data_struct)
 %
%  SYNTAX:  [resvv,vals,vecs] = calc_vv_resvals(tok_data_struct)
%
%  PURPOSE:  Calculate one-turn resistance of vacuum vessel from 
%		tok_data_struct objects, and also calculate eigenvalues/
%		vectors from the vessel-only passive system. Prints
%		various values to standard output and passes to output
%		variables.
%
%  INPUTS:
%	tok_data_struct = standard TokSys passive system data structure
%		which must include mvv, resv
%
%  OUTPUTS:
%	resvv = one-turn resistance of VV [Ohms]
%	vals = eigenvalues of vessel system amatvv=-inv(mvv)*diag(resv);
%	vecs = eigenvectors of vessel system amatvv=-inv(mvv)*diag(resv);
%
%  RESTRICTIONS:
%
%  METHOD:  

%  WRITTEN BY:  Dave Humphreys 	ON	6/6/08
%
%  MODIFICATION HISTORY:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prelims:
   mu0 = 0.4*pi;

% Derived Values:


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Caculate 1-Turn VV Resistance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  resvv = 1/sum(1./tok_data_struct.resv);
  disp(' ')
  disp(['One-turn VV resistance = ',num2str(resvv),' Ohms'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct eigensystem and calc eigenvalues/vecs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  rvv = diag(tok_data_struct.resv);
  mvv = tok_data_struct.mvv;

  amatvv = -inv(mvv)*rvv;
  [vecs,vals] = eigsort(amatvv);

  idx = min(length(vals),10);

  disp(' ')
  disp(['Top ',int2str(idx),' eigenvalues of VV [rad/sec  sec]:'])
 
  tmp = [vals(1:idx) 1./vals(1:idx)];
  disp(tmp)


