 %
%  SYNTAX:  test_output_objs
%
%  PURPOSE:  Test consistency of EFIT equilibrium, output_objs data from PCS 
%	calculations, and calculations from objects in toksys vacuum data 
%	structure.  Only tested for D3D and NSTX right now.
%
%  INPUT: 
%	tok_system = data structure created by build_<device>_sys.m
%	tok_data_struct = vacuum data model structure
%	output_objs = data structure created by plasma_outputs.m
%	wait_for_plots = set to 1 to pause after each plot [0]
%	figure_num = figure number for incrementing figures [0]
%
%  OUTPUT:
%	 plots comparing various calculation results
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	9/1/09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)test_output_objs.m	1.1 09/02/09

if(exist('figure_num')~=1)
   figure_num=0;
end
if(exist('wait_for_plots')~=1)
   wait_for_plots=0;
end

if 0	% NSTX testing
   load temp11.mat
%   load out_objs_132185_300_new.mat
end

if 0  % For D3D testing:
   load temp12.mat
%   load out_objs_129870_2302_new.mat
end

equil_data = tok_system.equil_data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing on X point grids:

dr = tok_data_struct.rg(2)-tok_data_struct.rg(1);
dz = tok_data_struct.zg(2)-tok_data_struct.zg(1);
cphi = equil_data.jphi*1e6*dr*dz;

% Use cproj to reduce NSTX coil currents to only those values fitted.
% MAKE SURE that cproj ends up = identity for D3D:

cproj = zeros(size(output_objs.grid_green_fns.GRbgd_c,2),length(tok_system.cc0));
nec = length(tok_data_struct.ecnturn);
for k=1:nec
   cproj(k,k)=1;	% ohmic coil
end
for k=nec+1:size(cproj,1)
   idx = find(tok_data_struct.def_connect.fcid==k-nec);
   if(~isempty(idx))
      cproj(k,idx(1)+nec)=1;
   end
end

vproj = zeros(size(output_objs.grid_green_fns.GRbgd_v,2),length(tok_system.vc0));
vproj(1,1)=1;	% ohmic coil
for k=2:size(vproj,1)
   idx = find(tok_data_struct.def_connect.vvid==k-1);
   if(~isempty(idx))
      vproj(k,1+idx)=tok_data_struct.def_connect.vvfrac(idx);
   end
end

rbgd = unique(output_objs.rbgd);
zbgd = unique(output_objs.zbgd);
rbgd2 = unique(output_objs.rbgd2);
zbgd2 = unique(output_objs.zbgd2);
nr = length(unique(output_objs.rbgd));
nz = length(unique(output_objs.zbgd));
nr2 = length(unique(output_objs.rbgd2));
nz2 = length(unique(output_objs.zbgd2));

next_figure
psi = output_objs.grid_green_fns.Mbgd_c*cproj*tok_system.cc0 + ...
        output_objs.grid_green_fns.Mbgd_v*vproj*tok_system.vc0 + ...
            output_objs.grid_green_fns.Mbgd_p*cphi(:);
psi = reshape(psi,nz,nr);
plot_tok_geo(tok_data_struct,[],equil_data)
hold on
contour(rbgd,zbgd,psi,40)
hold off
axis equal
title('psi computed from output objects vs. efit equilibrium contours')
if(wait_for_plots)
   pause
end

next_figure
psi1 = output_objs.Psi0(sum(output_objs.nptiso)+[1:nr*nz]);
psi1 = reshape(psi1,nz,nr);
plot_tok_geo(tok_data_struct,[],equil_data)
hold on
contour(rbgd,zbgd,psi1,40)
hold off
axis equal
title('psi = Psi0 from output objects vs. efit equilibrium contours, lower grid')
if(wait_for_plots)
   pause
end

next_figure
psi2= output_objs.Psi0(sum(output_objs.nptiso)+nr*nz+[1:nr*nz]);
psi2 = reshape(psi2,nz2,nr2);
plot_tok_geo(tok_data_struct,[],equil_data)
hold on
contour(rbgd2,zbgd2,psi2,40)
hold off
axis equal
title('psi = Psi0 from output objects vs. efit equilibrium contours, upper grid')
if(wait_for_plots)
   pause
end

% TEST the details of individual coil responses onto the X point grid:

next_figure
temp = interp2(equil_data.rg,equil_data.zg,equil_data.psizr,rbgd,zbgd');
psi_interp = reshape(temp,nz,nr);	% <== rtefit equilibrium interpolated onto X pt grid
plot(psi_interp(:))
hold on
plot(psi1(:),'r--')
plot(psi_interp(:)-psi1(:),'c')
hold off
title('efit psi interpolated onto lower grid(b) vs. Psi0 from output_objs(r), and difference(c)')
if(wait_for_plots)
   pause
end

next_figure
temp = interp2(equil_data.rg,equil_data.zg,equil_data.psizr,rbgd2,zbgd2');
psi_interp = reshape(temp,nz,nr);	% <== rtefit equilibrium interpolated onto X pt grid
plot(psi_interp(:))
hold on
plot(psi2(:),'r--')
plot(psi_interp(:)-psi2(:),'c')
hold off
title('efit psi interpolated onto lower grid vs. Psi0 from output_objs')
if(wait_for_plots)
   pause
end

rg=tok_data_struct.rg; zg = tok_data_struct.zg;
if(size(tok_data_struct.def_connect.fcid,1)==1)
   ccid = [1 tok_data_struct.def_connect.fcid+1];
else
   ccid = [1;tok_data_struct.def_connect.fcid+1];
end
for k=1:length(unique(ccid))
   idx = find(ccid == k);
   dd = output_objs.grid_green_fns.Mbgd_c*cproj(:,idx(1))*tok_system.cc0(idx(1));
   d1 = tok_data_struct.mpc(:,idx)*tok_system.cc0(idx);

   next_figure
   plot(dd)
   hold on
   d1 = reshape(d1,33,33);
   temp = interp2(rg,zg,d1,rbgd,zbgd');
   plot(temp(:),'r--')
   hold off
   title(['computed individual coil response ' int2str(k) ' from output_objs (b) and vac (r)'])
   if(wait_for_plots)
      pause
   end

   ratio = dd./temp(:);
   if(max(abs(ratio - 1)) > 0.05)	% only plot of > 5% error
     next_figure
     plot(ratio)
     title(['ccnturn(' int2str(k) ')=' int2str(tok_data_struct.ccnturn(k))])
     title(['ratio of coil response ' int2str(k) ' from output_objs (b) to vac (r)'])
     if(wait_for_plots)
        pause
     end

     next_figure
     plot_tok_geo(tok_data_struct)
     hold on
     dd=reshape(dd,nz,nr);
     contour(rbgd,zbgd,dd,25,'b')
     contour(rg,zg,d1,50,'r--')                           
     contour(rbgd,zbgd,temp,25,'m')
     axis equal
     title(['computed individual coil response ' int2str(k) ' from output_objs (b) and vac (r)'])
     hold off
     if(wait_for_plots)
        pause
     end
   end
end

% END testing for X point grids.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test that isoflux segment flux matches the efit equilibria:

%nr = length(unique(output_objs.rbgd));
%nz = length(unique(output_objs.zbgd));
%rg=tok_data_struct.rg; zg = tok_data_struct.zg;
%if(size(tok_data_struct.def_connect.fcid,1)==1)
%   ccid = [1 tok_data_struct.def_connect.fcid+1];
%else
%   ccid = [1;tok_data_struct.def_connect.fcid+1];
%end

nptiso = output_objs.nptiso;
ntot=0;
for k=1:size(nptiso,1)
  if(nptiso(k)>0)
     idx = ntot+1:ntot+nptiso(k);
     riso = output_objs.riso(idx);
     ziso = output_objs.ziso(idx);
     psi_interp = interp2(equil_data.rg,equil_data.zg,equil_data.psizr,riso,ziso);

     next_figure
     plot(output_objs.Psi0(idx))
     hold on
     plot(psi_interp,'r--')
     hold off
     title(['Segment ' int2str(k) ':Psi0 (b) vs. interpolated from equil_data (r)'])
     if(wait_for_plots)
        pause
     end

if 0
     for j=1:length(unique(ccid))
        idx1 = find(ccid == j);
        dd=output_objs.grid_green_fns.Miso_c(idx1,:)*cproj(:,idx1(1))*tok_system.cc0(idx1(1));
        d1 = tok_data_struct.mpc(:,idx1)*tok_system.cc0(idx1);

        next_figure
        plot(dd)
        hold on
        d1 = reshape(d1,33,33);
        temp = interp2(rg,zg,d1,rbgd,zbgd');
        plot(temp(:),'r--')
        plot(psi_interp(:)*2*pi,'g-.')
        hold off
        title('dd(b) vs. d1(r)')
     end
end

     ntot = ntot+nptiso(k);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0

next_figure
plot_tok_geo(tok_data_struct,[],equil_data);
hold on
contour(rbgd,zbgd,temp,50)

% NOW see if you can make the same plot from the output_objs data

% from Jim's ohmic_dist1.m code:

%  [dfdr,dfdz]= gradient(psizr,rg(2)-rg(1),zg(2)-zg(1));
%  
%  brr= -dfdz./rgg/(2*pi); % br= 1/(2pi*r)*df/dz
%  bzz= +dfdr./rgg/(2*pi);
%  bbb= sqrt(brr.^2 + bzz.^2);
  

Br = output_objs.grid_green_fns.GRbgd_c*cproj*tok_system.cc0 + ...
	output_objs.grid_green_fns.GRbgd_v*vproj*tok_system.vc0 + ...
	    output_objs.grid_green_fns.GRbgd_p*cphi(:);
Bz = output_objs.grid_green_fns.GZbgd_c*cproj*tok_system.cc0 + ...
	output_objs.grid_green_fns.GZbgd_v*vproj*tok_system.vc0 + ...
	    output_objs.grid_green_fns.GZbgd_p*cphi(:);

modB = sqrt(Br.^2 + Bz.^2);
modB = reshape(modB,nz,nr);

plot_tok_geo(tok_data_struct)
hold on
[c,h]= contour(rbgd,zbgd,modB,40);
if ~isempty(c) clabel(c,h); end
axis equal
hold off

% NSTX: in indiv. coil code, we get that the factor dd/temp(:) corresponds to ccnturn(k)
%>> [output_objs.cc0 output_objs.ccnturn]
%
%   -0.0145  240.0000	- corresponds to factor dd/temp(:)
%    0.0131   20.0000	- corresponds to factor dd/temp(:)
%    0.0001   28.0000	- etc ...
%   -0.0033   30.0000
%         0   17.0000
%   -0.0082   24.0000
%   -0.0082   24.0000
%         0   17.0000
%   -0.0049   30.0000
%    0.0037   28.0000
%    0.0078   20.0000
%    0.0000   32.0000
%         0    1.0000
%         0    1.0000
%         0    1.0000
%         0    1.0000
%         0   48.0000
%         0   48.0000
%
tt = cproj*tok_system.cc0

%[tt output_objs.cc0*1e6] = 1.0e+04 *
%   -0.3617   -1.4514
%    1.3081    1.3081
%    0.0138    0.0138
%   -0.3267   -0.3267
%         0         0
%   -0.8199   -0.8199
%   -0.8199   -0.8199
%         0         0
%   -0.4856   -0.4856
%    0.3737    0.3737
%    0.7772    0.7772
%    0.0002    0.0002
%         0         0
%         0         0
%         0         0
%         0         0
%         0         0
%         0         0
%
% Test reproduction of the full flux map - this works OK for NSTX:

fullpsi = tok_data_struct.mpc * tok_system.cc0 + tok_data_struct.mpv * tok_system.vc0 + ...
            mpp_x_vec(tok_data_struct.mpp,cphi(:));
fullpsi = reshape(fullpsi,tok_data_struct.nz,tok_data_struct.nr);
[c,h] = plot_tok_geo(tok_data_struct,[],equil_data);
hold on
contour(tok_data_struct.rg,tok_data_struct.zg,fullpsi,60,'g')
hold off
axis equal


   dd=output_objs.grid_green_fns.Mbgd_c*cproj(:,idx(1))*tok_system.cc0(idx(1));
   d1 = tok_data_struct.mpc(:,idx)*tok_system.cc0(idx);
for k=1:length(tok_system.cc0)
   idx = find(ccid == k);
   mpc = tok_data_struct.mpc(:,idx);
   mpc = reshape(mpc,33,33);
   mpc_temp = interp2(rg,zg,mpc,rbgd,zbgd');
   size(mpc_temp)
   Mbgd_c = output_objs.grid_green_fns.Mbgd_c*cproj(:,idx(1));
   Mbgd_c=reshape(Mbgd_c,18,19);
end

% Test computed response in matlab with stored PCS response (computed by IDL):
name = 'bot_3cm_grid_00041'

[sizedata,cntrl_pts,mutuals,brg,bzg] = read_response(name,6565);

sizedata1=sizedata;
brg1=brg;
bzg1=bzg;
mutuals1 = mutuals;
cntrl_pts1 = cntrl_pts;

load /m/GAtools/tokamaks/d3d/make/after_ADP/d3d_obj_mks_struct_6565.mat
load /m/GAtools/tokamaks/d3d/make/after_ADP/d3d_obj_mks_struct_RTEFIT.mat

name = 'bot_3cm'
[sizedata,control_pts,mutuals,brgreens,bzgreens] = ...
              calc_isoflux_response(name,'/pcshome/walker/pcs_runsa',tok_data_struct);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
