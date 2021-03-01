%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gscp_converge
%
%  PURPOSE: Converge circular plasma model to solution
%
%  INPUTS: The work space for gscp codes and:
%          max_iterations (default = 99)
%          max_psibarerr (default 1e-6)
%
	
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander  ON	1/22/14
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('max_iterations','var')
  max_iterations = 99;
end

gscp_analysis_response

iteration_counter = 0;
done = 0;
while ~done
  iteration_counter = iteration_counter+1;
  gs_find_dx
  dp = dpdx*dx;
  fc = abs(psimag-psibry);
  f = min(1,fc/max(abs(dp(1:nrh)))/2);
  f = min(f,dz/abs(dp(nrh+2))/2);
  f = min(f,a0/abs(dp(nrh+3))/5);
  % Stay inside grid
  if r0 + dp(nrh+1)*f + a0 + dp(nrh+3)*f > rg(nr-1)
    f = (rg(nr-1) - r0 - a0)/(dp(nrh+1) + dp(nrh+3));
  end
  if r0 + dp(nrh+1)*f - a0 - dp(nrh+3)*f < rg(2)
    f = (rg(2) - r0 + a0)/(dp(nrh+1) - dp(nrh+3));
  end
  if max(abs(psiherr(:))) < fc
    psih = psih + dp(1:nrh)'*f;
    r0 = r0 + dp(nrh+1)*f;
    z0 = z0 + dp(nrh+2)*f;
    a0 = a0 + dp(nrh+3)*f;
    sp = sp + dx(indsp)*f;
    sf = sf + dx(indsf)*f;
    ic    = ic    + dx(indic)*f;
    iv    = iv    + dx(indiv)*f;
  else
    psih = (psih + psihpla + psihapp)/2;
  end
  if a0 < r0/1e2
    dum = r0/1e2 - a0;
    r0 = r0 + dum;
    a0 = a0 + dum;
  end
  gscp_analysis_response
  if plotit > 1
    gs_plot_progress
  end
  
  costfun = max(abs(psiherr(:)/(psibry-psimag)));
  if evolve_option == 0 & exist('xc','var')
    if (constraints == 1 | constraints == 2) & isfield(xc,'cpasma')
      if xc.cpasma ~= 0
        costfun = max(costfun,abs(xc.cpasma-cpasma)/abs(xc.cpasma));
      else
        costfun = max(costfun,abs(xc.cpasma-cpasma));
      end
    end
    if (constraints == 1 | constraints == 2) & ...
       isfield(xc,'li') & xc.li > 0
      costfun = max(costfun,abs(xc.li-li)/abs(xc.li));
    end
    if (constraints == 1 | constraints == 2) & ...
       isfield(xc,'betap') & xc.betap > 0
      costfun = max(costfun,abs(xc.betap-betap)/abs(xc.betap));
    end
  elseif evolve_option == 9 & exist('xs','var')
    costfun = max(costfun,abs(xs(end-2)-cpasma));
  end
  
  done = costfun < 1e-6 | iteration_counter >= max_iterations;
end
