function [dpsizr, dpsimag, dpsibry] = gspert_dpsizr(response,perturbation,tok_data_struct,eq0,idoplot,options)
%
%  USAGE: [dpsizr, dpsimag, dpsibry] = gspert_dpsizr(response,perturbation,tok_data_struct,eq0,idoplot,options)
%
%  PURPOSE: Return the perturbed flux on the grid, dpsizr,
%           given perturbation of set of independent variables defined in gspert
%
%  INPUTS:
%    response = the plasma response calculated by gspert
%    perturbation = row vector of perturbations of the independent variables in gspert
%      May contain several columns with different perturbations 
%      If the gspert response was calculated with iconstraints = 1
%        the columns should contain [dis; dbetap; dli; dip]
%      If the gspert response was calculated with iconstraints = 2
%        the columns should contain [dis; dw; dli; dip]
%      If the gspert response was calculated with iconstraints = 3
%        the columns should contain [dis; dW; dI]
%      Here, 'is' is conductor currents (coils and vessel)
%        Should be compatible with tok_data_struct so that dpsizr_app = [mpc mpv]*dis
%      w is total thermal energy
%      W is thermal energy in a few profile points (help gspert for more info)
%      I is toroidal current in a few profile points (help gspert for more info)
%    tok_data_struct is standard toksys description of tokamak, including mutual inductances
%      should be the same that was used in gspert to obtain response
%    eq0 is the (unperturbed) equilibrium that response was calculated from
%    idoplots = (optional) if scalar, show contours of eq.psizr+dpsizr for idoplot seconds: 
%       default is 0 = don't plot
%       idoplots can be an array to plot only select time samples
%    options = optional options:
%      options.trace: if ~0, boundary is traced to make small correction of predicted dpsibry.
%        This feature is only used in the contour plot
%
%  OUTPUTS:
%    dpsizr = the perturbed flux
%      If perturbation is only one time sample then dpsizr is a (nz,nr) matrix
%      If perturbation is matrix with nt time samples, dpsizr is a (nr*nz,nt) matrix
%
%  RESTRICTIONS:
%
%  METHOD:
	
%  VERSION @(#)gspert_dpsizr.m	1.2 05/09/11
%
%  WRITTEN BY:  Anders Welander  ON	1/12/11
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  ns = size(response.dcphidis,2);
  
  if ~exist('options','var')
    options.trace = 0;
  end
  if ~isfield(options,'trace')
    options.trace = 0;
  end
  
  if isfield(response,'dcphidbetap')
    dcphidpert = [response.dcphidis response.dcphidbetap(:) response.dcphidli(:) response.dcphidip(:)];
    dpsimagdpert = [response.dpsimagdis response.dpsimagdbetap response.dpsimagdli response.dpsimagdip];
    dpsibrydpert = [response.dpsibrydis response.dpsibrydbetap response.dpsibrydli response.dpsibrydip];
  elseif isfield(response,'dcphidw')
    dcphidpert = [response.dcphidis response.dcphidw(:) response.dcphidli(:) response.dcphidip(:)];
    dpsimagdpert = [response.dpsimagdis response.dpsimagdw response.dpsimagdli response.dpsimagdip];
    dpsibrydpert = [response.dpsibrydis response.dpsibrydw response.dpsibrydli response.dpsibrydip];
  elseif isfield(response,'dcphidW')
    dcphidpert = [response.dcphidis response.dcphidW response.dcphidI];
    dpsimagdpert = [response.dpsimagdis response.dpsimagdW response.dpsimagdI];
    dpsibrydpert = [response.dpsibrydis response.dpsibrydW response.dpsibrydI];
  else
    disp('Que? Something is wrong with response.')
  end
  
  struct_to_ws(tok_data_struct);
  Mgg = zeros(nr*nz,nr*nz);
  for j = 1:nr, for k = 1:nz, for l=1:nr, for m=1:nz
    Mgg(k+(j-1)*nz,m+(l-1)*nz) = mpp(1+abs(m-k)+(j-1)*nz,l);
  end, end, end, end
  
  nt = size(perturbation,2); % Number of time samples
  if exist('idoplot','var') & length(idoplot)==1, idoplot = idoplot*ones(nt,1); end
  idoplot(nt+1) = 0;
  
  for j = 1:nt
    dpsizr(:,j) = [mpc mpv]*perturbation(1:ns,j) + Mgg*(dcphidpert*perturbation(:,j));
    dpsimag(j) = dpsimagdpert*perturbation(:,j);
    dpsibry(j) = dpsibrydpert*perturbation(:,j);
    if idoplot(j)
      clf, hold on
      plot(limdata(2,:),limdata(1,:),'k','linew',2)
      if options.trace
        [rb,zb,rx,zx,ra,za,r0,z0,ilimited,psimag,psibry]=trace_boundary(eq0.psizr+reshape(dpsizr(:,j),nz,nr),rg,zg,limdata);
	disp(['gspert-predicted dpsibry = ' num2str(dpsibry) ', displacement-corrected value = ' num2str(psibry-eq0.psibry)])
        contour(rg,zg,eq0.psizr+reshape(dpsizr(:,j),nz,nr),linspace(psibry,psimag,8))
      else
        contour(rg,zg,eq0.psizr+reshape(dpsizr(:,j),nz,nr),linspace(eq0.psibry+dpsibry(j),eq0.psimag+dpsimag(j),8))
      end
      title(['Perturbation # ' num2str(j)])
      axis('image')
      if idoplot(j)>0
        pause(idoplot(j))
      else
        disp('Press return to continue plotting. ');
	pause
      end
    end
  end
  
  if nt == 1
    dpsizr = reshape(dpsizr,nz,nr);
  end
