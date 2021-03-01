%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gsdesign_iss
%
%  PURPOSE: Design ITER-Similiar-Shape for generic tokamak
%           Demonstrates gsdesign capabilities by designing
%           the ITER shape for ITER and 9 additional tokamaks
%           No coil connections are specified
%           No coil current limits specified but gsdesign has
%           defaults for some tokamaks
%
%  INPUTS:  none
%
%  OUTPUTS: eq, equilibrium for specified tokamak
	
%
%  WRITTEN BY:  Anders Welander ON 1/30/15
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('iter','var')
  gsdesign_demo_iter
  set(gcf,'Name',['Designed shape for ' upper(config.tokamak)])
  spec.targets.rbdef = spec.targets.rx;
  spec.targets.zbdef = spec.targets.zx;
  iter.eq = eq;
  iter.spec = spec;
  iter.rsep = spec.targets.rsep;
  iter.zsep = spec.targets.zsep;
  pause(5)
  kbbbs = eq.zbbbs ~= min(eq.zbbbs);
  [~,gaps,gapspec] = calc_gaps([eq.rbbbs(kbbbs) eq.zbbbs(kbbbs)],config,[],eq);
end
tokamaks = char(...
  'pegasus/make/pegasus_obj_mks_struct.mat', ...        % 1. Pegasus
  'sst/make/sst_obj_mks_struct.mat', ...                % 2. SST
  'fdf/make/standard/fdf_obj_mks_struct.mat', ...       % 3. FDF
  'd3d/make/after_ADP/d3d_obj_mks_struct_3333.mat', ... % 4. D3D
  'east/make/east_obj_2014_3333.mat', ...               % 5. EAST
  'kstar/make/kstar_obj_mks_struct_2014_6565.mat', ...  % 6. KSTAR
  'nstx-u/make/15feb/nstxu_obj_15feb_6565.mat', ...     % 7. NSTX-U
  'cfetr/make/cfetr_obj_2014_3333.mat', ...             % 8. CFETR
  'hl2m/make/hl2m_obj_2014_3333.mat');                  % 9. HL2M
for itok = 1:9
  figure_handle = figure;
  if isa(figure_handle,'matlab.ui.Figure')
    figure_handle = figure_handle.Number;
  end
  fpos = get(figure_handle,'position');
  fpos(4) = fpos(3)/2;
  fpos(3:4) = [1000 500];
  set(figure_handle,'position',fpos)
  eq = iter.eq;
  load(['/m/GAtools/tokamaks/' deblank(tokamaks(itok,:))])

  config = regrid(65,65,tok_data_struct);

  % Plot ITER with designed equilibrium
  ix = 0; % index of x-point that defines boundary
  for i = 1:eq.nulls.count
    if eq.nulls.r(i) == eq.rbdef & eq.nulls.z(i) == eq.zbdef
      ix = i;
    end
  end
  clf
  set(zoom(gca),'ActionPostCallback','');
  hold on
  plot(eq.xlim,eq.ylim,'k','linew',6)
  patch(eq.xlim,eq.ylim,[0.8 0.8 0.8])
  patch(eq.rbbbs(1:eq.nbbbs),eq.zbbbs(1:eq.nbbbs),[1 1 0.9])
  psia = eq.psimag;
  psib = eq.psibry;
  nhires = 10;
  nzhr = 1+nhires*(eq.nh-1);
  nrhr = 1+nhires*(eq.nw-1);
  rghr = linspace(eq.rg(1),eq.rg(eq.nw),nrhr);
  zghr = linspace(eq.zg(1),eq.zg(eq.nh),nzhr)';
  px = gs_interp2(eq.rg,eq.zg,eq.psizr,ones(nzhr,1)*rghr, zghr*ones(1,nrhr));
  %contour(rghr,zghr,1-(px-psia)/(psib-psia),[0:7]/8,'linew',3)
  for i = eq.nulls.count:-1:1
    if i ~= ix
      plot(eq.nulls.r(i),eq.nulls.z(i),'bx','markersize',12,'linew',2)
    end
  end
  plot(iter.spec.targets.rx,iter.spec.targets.zx,'yx','markersize',18,'linew',3)
  plot(iter.rsep,iter.zsep,'mo','markersize',9,'linew',2)
  plot(eq.rbbbs(1:eq.nbbbs),eq.zbbbs(1:eq.nbbbs),'linew',3,'color',[0,0, 0.5625])
  plot(eq.rbdef,eq.zbdef,'yx','markers',24,'linew',5)
  plot(eq.rbdef,eq.zbdef,'rx','markers',21,'linew',3)
  axis('image')
  axis([0 eq.rg(end) eq.zg(1) eq.zg(end)])
  title(['ISS for ' upper(config.tokamak)]);

  % Plot generic tokamak
  if size(config.limdata,2) > 2
    config.limdata = config.limdata';
  end
  if min(config.limdata(:,2)) < min(config.limdata(:,1))
    config.limdata = config.limdata(:,[2 1]);
  end
  config.Rlim = config.limdata([1:end 1],2);
  config.Zlim = config.limdata([1:end 1],1);
  plot(config.Rlim,config.Zlim,'k','linew',6)
  patch(config.Rlim,config.Zlim,[0.8 0.8 0.8])
  h = patch(eq.rbbbs(1:eq.nbbbs),eq.zbbbs(1:eq.nbbbs),[1 1 0.9]);
  %hc = contour(rghr,zghr,1-(px-psia)/(psib-psia),[0:7]/8,'linew',3);
  h(end+1) = plot(iter.spec.targets.rx,iter.spec.targets.zx,'yx','markersize',18,'linew',3);
  h(end+1) = plot(iter.rsep,iter.zsep,'mo','markersize',2,'linew',2);
  h(end+1) = plot(eq.rbbbs(1:eq.nbbbs),eq.zbbbs(1:eq.nbbbs),'linew',3,'color',[0,0,0.5625]);
  h(end+1) = plot(eq.rbdef,eq.zbdef,'yx','markers',24,'linew',5);
  h(end+1) = plot(eq.rbdef,eq.zbdef,'rx','markers',21,'linew',3);

  scale = 1;
  r_offset = 0;
  z_offset = 0;
  r = eq.rbbbs;
  z = eq.zbbbs;

  rmin = min(r);
  rmax = max(r);
  zmin = min(z);
  zmax = max(z);
  r2min = min(config.Rlim);
  r2max = max(config.Rlim);
  z2min = min(config.Zlim);
  z2max = max(config.Zlim);

  scale = min((r2max-r2min)/(rmax-rmin),(z2max-z2min)/(zmax-zmin));
  r_offset = (r2min-scale*rmin+r2max-scale*rmax)/2;
  z_offset = (z2min-scale*zmin+z2max-scale*zmax)/2;

  n = 50;
  ss = linspace(1,scale,n);
  rs = linspace(0,r_offset,n);
  zs = linspace(0,z_offset,n);
  for i = 1:n
    scale = ss(i);
    r_offset = rs(i);
    z_offset = zs(i);
    set(h(1),'XData',r_offset + scale*eq.rbbbs)
    set(h(1),'YData',z_offset + scale*eq.zbbbs)
    set(h(3),'XData',r_offset + scale*iter.rsep)
    set(h(3),'YData',z_offset + scale*iter.zsep)
    ymin = min(min(z_offset + scale*iter.zsep),min(config.Zlim));
    ymax = max(max(z_offset + scale*iter.zsep),max(config.Zlim));
    pos = get(gcf,'Position');
    y = max(0,(eq.rg(end)*pos(4)/pos(3)-(ymax-ymin))/2);
    %axis([0 eq.rg(end) ymin-y ymax+y])
    pause(0.04)
  end
  ax = axis;
  x1 = linspace(ax(1),config.rg(1),n);
  x2 = linspace(ax(2),config.rg(end),n);
  y1 = linspace(ax(3),config.zg(1),n);
  y2 = linspace(ax(4),config.zg(end),n);
  for i = 1:n
    axis([x1(i) x2(i) y1(i) y2(i)])
    pause(0.04)
  end
  eq2 = eq;

  for i = 1:22
    eq2.rbbbs = r_offset+scale*eq.rbbbs;
    eq2.zbbbs = z_offset+scale*eq.zbbbs;
    eq2.rg = r_offset+scale*eq.rg;
    eq2.zg = z_offset+scale*eq.zg;
    eq2.dr = scale*eq.dr;
    eq2.dz = scale*eq.dz;
    [~,gaps2] = calc_gaps([eq2.rbbbs(kbbbs) eq2.zbbbs(kbbbs)],config,[],eq2);

    scale = scale+1e-3;
    eq2.rbbbs = r_offset+scale*eq.rbbbs;
    eq2.zbbbs = z_offset+scale*eq.zbbbs;
    eq2.rg = r_offset+scale*eq.rg;
    eq2.zg = z_offset+scale*eq.zg;
    eq2.dr = scale*eq.dr;
    eq2.dz = scale*eq.dz;
    [~,gaps2b] = calc_gaps([eq2.rbbbs(kbbbs) eq2.zbbbs(kbbbs)],config,[],eq2);
    scale = scale-1e-3;
    dgapsdscale = (gaps2b-gaps2)/1e-3;

    r_offset = r_offset+1e-3;
    eq2.rbbbs = r_offset+scale*eq.rbbbs;
    eq2.zbbbs = z_offset+scale*eq.zbbbs;
    eq2.rg = r_offset+scale*eq.rg;
    eq2.zg = z_offset+scale*eq.zg;
    eq2.dr = scale*eq.dr;
    eq2.dz = scale*eq.dz;
    [~,gaps2b] = calc_gaps([eq2.rbbbs(kbbbs) eq2.zbbbs(kbbbs)],config,[],eq2);
    r_offset = r_offset-1e-3;
    dgapsdr_offset = (gaps2b-gaps2)/1e-3;

    z_offset = z_offset+1e-3;
    eq2.rbbbs = r_offset+scale*eq.rbbbs;
    eq2.zbbbs = z_offset+scale*eq.zbbbs;
    eq2.rg = r_offset+scale*eq.rg;
    eq2.zg = z_offset+scale*eq.zg;
    eq2.dr = scale*eq.dr;
    eq2.dz = scale*eq.dz;
    [~,gaps2b] = calc_gaps([eq2.rbbbs(kbbbs) eq2.zbbbs(kbbbs)],config,[],eq2);
    z_offset = z_offset-1e-3;
    dgapsdz_offset = (gaps2b-gaps2)/1e-3;

    m = [dgapsdscale; dgapsdr_offset; dgapsdz_offset]';
    ngaps = length(gaps);
    w = ones(ngaps,3);
    dgaps = (gaps2-gaps*scale)';
    k = gaps2 < min(gaps*scale);
    w(k,:) = 10;
    if sum(k) > 1
      v = -pinv(w.*m)*(w(:,1).*dgaps);
      if max(abs(v)) > eq2.dr
	v = v/max(abs(v))*eq2.dr;
      end
    else
      v = zeros(3,1);
    end

    scale = scale + v(1)/2;
    r_offset = r_offset + v(2)/2;
    z_offset = z_offset + v(3)/2;

    set(h(1),'XData',r_offset + scale*eq.rbbbs)
    set(h(1),'YData',z_offset + scale*eq.zbbbs)
    set(h(3),'XData',r_offset + scale*iter.rsep)
    set(h(3),'YData',z_offset + scale*iter.zsep)

    drawnow
  end

  config.constraints = 1;
  config.psikn = [0 0.40 0.70 1];
  config.max_iterations = 111;
  config.no_edge_current = true;
  config.no_edge_gradient = true;
  spec = iter.spec;
  spec.locks.ic = nan;
  spec.targets.rsep = r_offset + scale*iter.spec.targets.rsep;
  spec.targets.zsep = z_offset + scale*iter.spec.targets.zsep;
  spec.targets.rx = r_offset + scale*iter.spec.targets.rx;
  spec.targets.zx = z_offset + scale*iter.spec.targets.zx;
  spec.targets.cpasma = scale^2*iter.spec.targets.cpasma;
  dum = 10^floor(log10(spec.targets.cpasma));
  spec.targets.cpasma = round(2*spec.targets.cpasma/dum)/2*dum;
  spec.targets.ic = zeros(config.nc,1);
  spec.weights.ic = ones(config.nc,1)*1e-5;
  spec.cccirc = 1:config.nc;
  spec.targets.rbdef = r_offset + scale*iter.spec.targets.rbdef;
  spec.targets.zbdef = z_offset + scale*iter.spec.targets.zbdef;
  spec.weights.bdef = 0.0;
  %spec.targets.rsep = r_offset + scale*iter.spec.targets.rbdef;
  %spec.targets.zsep = z_offset + scale*iter.spec.targets.zbdef;
  clear gsdesign
  eq2.ic = zeros(config.nc,1);
  eq2.iv = zeros(config.nv,1);
  eq2.rmaxis = r_offset+scale*eq.rmaxis;
  eq2.zmaxis = z_offset+scale*eq.zmaxis;
  eq2.psizr = scale*eq.psizr;
  eq2.pres = scale^2*eq.pres;
  eq2.fpol = scale^2*eq.fpol;
  spec.fig = figure_handle;
  eq = gsdesign(spec,[],config);

  if ishandle(spec.fig)
    set(spec.fig,'Name',['Designed ISS for ' upper(config.tokamak)])
  end

end % End itok
