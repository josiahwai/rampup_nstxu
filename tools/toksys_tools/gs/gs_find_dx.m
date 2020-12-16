%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gs_find_dx
%
%  PURPOSE: Find change of x (dx) since last analyzed equilibrium (lae)
%
%  INPUTS:  One of either xc or xs
%           xc, array or structure, may contain:
%               ic, iv, sp, sf, cpasma, li, betap
%               (content depending on constraints)
%           xs, state vector (content controlled by evolve_option)
%
%  OUTPUTS: dx, the change of the state vector, x (=[ic;iv;sp;sf;er])
%           dcpasma, dli, dbetap, dWth, changes of state plasma quantities
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander  ON	1/7/14
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if evolve_option == 0 % No evolving. The new x is specified with xc
  if exist('xc','var') & isstruct(xc)
    if isfield(xc,'ic')
      dic = xc.ic - lae.ic;
    elseif isfield(xc,'cc0t')
      dic = xc.cc0t - lae.ic;
    else
      dic = zeros(nc,1);
    end
    if isfield(xc,'iv')
      div = xc.iv - lae.iv;
    elseif isfield(xc,'vc0t')
      div = xc.vc0t - lae.iv;
    else
      div = zeros(nv,1);
    end
    if constraints == 0
      if isfield(xc,'sp')
	dsp = xc.sp - lae.sp;
      else
	dsp = zeros(nkn+2,1);
      end
      if isfield(xc,'sf')
	dsf = xc.sf - lae.sf;
      else
	dsf = zeros(nkn+2,1);
      end
      dx = [dic; div; dsp; dsf; -er];
      plasma = any(lae.sf+dsf) | any(lae.sp+dsp);
    elseif constraints == 1
      if isfield(xc,'ip')
	dip = xc.ip - lae.cpasma;
      elseif isfield(xc,'cpasma')
	dip = xc.cpasma - lae.cpasma;
      else
	dip = 0;
      end
      if isfield(xc,'li')
	dli = xc.li - lae.li;
      else
	dli = 0;
      end
      if isfield(xc,'betap')
	dbetap = xc.betap - lae.betap;
      else
	dbetap = 0;
      end
      dxc = [dic; div; dip; dli; dbetap; -er];
      dx = dxdxc*dxc;
      plasma = lae.cpasma+dip ~= 0;
      if ~lae.plasma
        cpasma = dip;
      end
    elseif constraints == 2
      if isfield(xc,'ip')
	dip = xc.ip - lae.cpasma;
      elseif isfield(xc,'cpasma')
	dip = xc.cpasma - lae.cpasma;
      else
	dip = 0;
      end
      if isfield(xc,'li')
	dli = xc.li - lae.li;
      else
	dli = 0;
      end
      if isfield(xc,'Wth')
	dWth = xc.Wth - lae.Wth;
      else
	dWth = 0;
      end
      dxc = [dic; div; dip; dli; dWth; -er];
      dx = dxdxc*dxc;
      plasma = lae.cpasma+dip ~= 0;
      if ~lae.plasma
        cpasma = dip;
      end
    end
  elseif exist('xc','var') % xc is a vector
    if length(xc) >= nc
      dic = xc(indic) - lae.ic;
    else
      dic = zeros(nc,1);
    end
    if length(xc) >= nc+nv
      div = xc(indiv) - lae.iv;
    else
      div = zeros(nv,1);
    end
    if constraints == 0
      if length(xc) >= indsp(end)
	dsp = xc(indsp) - lae.sp;
      else
	dsp = zeros(nkn+2,1);
      end
      if length(xc) >= indsf(end)
	dsf = xc(indsf) - lae.sf;
      else
	dsf = zeros(nkn+2,1);
      end
      dx = [dic; div; dsp; dsf; -er];
      plasma = any(lae.sf+dsf) | any(lae.sp+dsp);
    elseif constraints == 1
      if length(xc) >= nc+nv+1
	dip = xc(nc+nv+1) - lae.cpasma;
      else
	dip = 0;
      end
      if length(xc) >= nc+nv+2
	dli = xc(nc+nv+2) - lae.li;
      else
	dli = 0;
      end
      if length(xc) >= nc+nv+3
	dbetap = xc(nc+nv+3) - lae.betap;
      else
	dbetap = 0;
      end
      dxc = [dic; div; dip; dli; dbetap; -er];
      dx = dxdxc*dxc;
      plasma = lae.cpasma+dip ~= 0;
    elseif constraints == 2
      if length(xc) >= nc+nv+1
	dip = xc(nc+nv+1) - lae.cpasma;
      else
	dip = 0;
      end
      if length(xc) >= nc+nv+2
	dli = xc(nc+nv+2) - lae.li;
      else
	dli = 0;
      end
      if length(xc) >= nc+nv+3
	dWth = xc(nc+nv+3) - lae.Wth;
      else
	dWth = 0;
      end
      dxc = [dic; div; dip; dli; dWth; -er];
      dx = dxdxc*dxc;
      plasma = lae.cpasma+dip ~= 0;
      if ~lae.plasma
        cpasma = dip;
      end
    end    
  else % xc doesn't exist  
    dx = [zeros(nx-1,1); -er];
    plasma = any(lae.sp) | any(lae.sf);  
  end  
else % evolve_option > 0, so xs is used to specify new x
  if exist('xs','var') & ~isempty(xs)
    if constraints == 0
      dx = [Pxx*xs; 0]-[lae.ic; lae.iv; lae.sp; lae.sf; er];
      if ~lae.plasma
        if any(dx([indsp indsf])) & cpasma == 0;
	  cpasma = 1e3;
	end
      end
      plasma = any([lae.sp+dx(indsp); lae.sf+dx(indsf)]);
    elseif constraints == 1
      dx = dxdxc*([Pxx*xs; 0]-[lae.ic; lae.iv; lae.cpasma; lae.li; lae.betap; er]);
      if ~lae.plasma
        cpasma = xs(end-2);
      end
      plasma = xs(end-2) ~= 0;
    elseif constraints == 2
      dx = dxdxc*([Pxx*xs; 0]-[lae.ic; lae.iv; lae.cpasma; lae.li; lae.Wth; er]);
      if ~lae.plasma
        cpasma = xs(end-2);
      end
      plasma = xs(end-2) ~= 0;
    end
  else  
    dx = [zeros(nx-1,1); -er];
  end
end
