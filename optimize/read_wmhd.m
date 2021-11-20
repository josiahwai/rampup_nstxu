
function wmhd = read_wmhd(eqs, tok_data_struct)

struct_to_ws(tok_data_struct);
dr = mean(diff(rg));
dz = mean(diff(zg));
dA = dr*dz;

if ~isfield(eqs, 'gdata')
  eqs.gdata = eqs;
end

neq = length(eqs.gdata);

for ieq = 1:neq

  eq = eqs.gdata(ieq);
      
  in = inpolygon(rgg, zgg, eq.rbbbs, eq.zbbbs);
  psin = linspace(0, 1, length(rg));
  psin_grid = (eq.psizr - eq.psimag) / (eq.psibry - eq.psimag);  
  pres_grid = interp1(psin, eq.pres, psin_grid(:), 'linear');
  pres_grid(~in(:)) = 0;

  wmhd(ieq,:) = nansum(3*pi*rgg(:).*pres_grid * dA);
end




