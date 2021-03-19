

Mvc = []; 
Mvv = [];
Mvp = [];
Rv = [];
Lv = [];

for i = 1:40
  load(['fit_vessel_' num2str(i) '.mat'])
  
  Mvc = [Mvc; vess_fit.Mvc];
  Mvv = [Mvv; vess_fit.Mvv];
  Mvp = [Mvp; vess_fit.Mvp];
  Rv = [Rv; vess_fit.Rv];
  Lv = [Lv; vess_fit.Lv];
end

Mvv = Mvv - diag(Mvv) + diag(Lv);

vessel_sysid_fit = variables2struct(Rv, Mvc, Mvv, Mvp);
save('/Users/jwai/Research/rampup_nstxu/sysid/vessel_sysid/fits/vessel_sysid_fit', 'vessel_sysid_fit')


% figure
% contourf( (Mvv - Mvv') / median(abs(Mvv(:))))
% colorbar
% 
% 
% % Mvv_sym = (Mvv + Mvv') / 2;


  












































