%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gs_limiter_radiation
%
%  PURPOSE: Calculate power per unit distance on limiter
%           from 1 W of isotropic radiation at grid elements
%
%  INPUTS:  rl, zl, ilimgg, rg, zg, rgg, zgg
%
%  OUTPUTS: Plg, radiated power to limiter from grid
%
%  METHOD:   
	
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander  ON	6/30/14
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Radiation in W/m2 at limiter points from radiating axisymmetric filaments
Plg = zeros(nl,ngg);

% Super simple approximation
for i = 1:nl
  Plg(i,:) = 1./((rgg(:)-rl(i)).^2+(zgg(:)-zl(i)).^2);
end

return
% Remove shadowed limiter points
for i = 1:ngg
  if ilimgg(i) == 0
    alg = angle(rgg(i)-rl + 1i*zgg(i)-1i*zl);
    dlg = (rgg(i)-rl).^2 + (zgg(i)-zl).^2;
    for j = 1:nl-1
      for k = 1:nl-2
        alg1 = alg(k);
	alg2 = alg(k+1);
	if alg1 > alg2+pi
	  alg1 = alg1-2*pi;
	end
	if alg1 < alg2-pi
	  alg1 = alg1+2*pi;
	end
	if alg1 > alg(j)+pi
	  alg1 = alg1-2*pi;
	  alg2 = alg2-2*pi;
	end
	if alg1 < alg(j)-pi
	  alg1 = alg1+2*pi;
	  alg2 = alg2+2*pi;
	end
        if (alg2-alg(j))*(alg1-alg(j)) <= 0
	  if j ~= k & j ~= k+1
	    if dlg(k) < dlg(j) | dlg(k+1) < dlg(j)
	      % limiter point j shadowed from grid i by limiter k:k+1
	      Plg(j,i) = 0;
	    end
	  end
	end
      end
    end
  end
end
