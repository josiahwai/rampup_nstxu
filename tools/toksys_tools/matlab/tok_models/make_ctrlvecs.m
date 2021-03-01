 %
%  USAGE:  >> make_ctrlvecs
%
%  PURPOSE: Script to construct control vectors for M-matrix in PCS.
%		Called from scripts such as
%		make_east_ctrlvecs or make_kstar_ctrlvecs.
%
%  INPUTS:  (Defined in workspace prior to calling make_predictors)
%     Required:
%       objdir = directory containing objects file
%       objfile = objects file name (eg kstar_obj_struct.mat)
%       tokamak = machine name, eg 'KSTAR', 'EAST', (needed for read_gfile_tok)
%	idxcr = PF coil index vector for Rctrlvec
%	idxcz = PF coil index vector for Zctrlvec
%	idxcip = PF coil index vector for Ipctrlvec
%       rrange1 = maj radius range of sub-grid for fields (m) eg [1.5 2]
%       zrange1 = vert range of sub-grid for fields (m) eg [-.5 .5]
%       iusejphi = 1 = use jphi (all current) region for fields, else ???
%       efit_gfile_grid = gfile to use to define grid from jphi if iusejphi=1
%	gfile_data = gfile data to use to define grid from jphi
%	(Only one of efit_gfile_grid or gfile_data should be specified.)
%	rrefk
%	zrefk
%    Optional:
%	separate_figures = set to 1 to get 1 figure per plot, otherwise multiple (default)
%	sign_ip = either +1/-1, assuming standard convention for all Ip, PF coils (default 1)
%	iconstraint = set to 1 to constrain Ip ctrl vector elements all one sign, else 0 (default)
%       ipltbp = 1 to plot Bprobes in plot_geo call (0 to not)
%       ipltfl = 1 to plot flux loops in plot_geo call (0 to not)
% USED???
%       nprtst = # of princ comp of grid to mag set to use in test
%       nsvdz1,2 = # of sing vals to keep in R pred inverse (for P*pred1,2)
%       nsvdr1,2 = # of sing vals to keep in R pred inverse
%       nsvdi1,2 = # of sing vals to keep in Ip pred inverse
%	nsvdk1,2 = # of sing vals to keep in Kappa pred inverse
%
%  OUTPUTS:
%	Zctrlvec = Z-ctrl vector (Vpf_Z = Zctrlvec*PID*Zerror)
%	Rctrlvec = R-ctrl vector (Vpf_R = Rctrlvec*PID*Zerror)
%	Ipctrlvec = Ip-ctrl vector (Vpf_R = Rctrlvec*PID*Zerror)
%	(Sign is such that multiplying ctrlvec by positive value => positive change
%		in controlled parameter, except kappa which has indeterminate sign.)
%    + many  plots evaluating quality of control vectors in producing specified field
%
%  RESTRICTIONS:
%     If iusejphi==1, must have efit_gfile_grid defined (which is source of
%       jphi on which grid is based if iusejphi=1)
%
%  METHOD:  
%	Variety of inversions of grid-to-mag objects, linear fits to selected
%	data from subgrids, subset of magnetics...

%  WRITTEN BY:  Dave Humphreys  ON	7/22/06
%
%  MODIFICATION HISTORY:
%	7/22/06 DAH Deriving make_ctrlvecs from make_predictors.m
%	8/21/06 DAH Removing normalization on ctrl vecs: want unity 
%		    field (T) or flux (Wb) from unity command. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure_num = 3;

if(~exist('sign_ip','var'))
   sign_ip = 1;		% default is positive Ip
end
if(~exist('iconstraint','var'))
   iconstraint = 0;		% default is to allow +/- in Ip control vector
end
if(~exist('separate_figures','var'))
   separate_figures = 1;
end

% Prelims and Constants:
   mu0 = 0.4*pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load, extract objects, select sub-grid, select diagnostics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   eval(['load ',objdir,objfile]);
   struct_to_ws(tok_data_struct);   %unpack objects structure to workspace

%Select sub-grid idxg to base predictor on:
  if iusejphi==1    %must have efit_gfile to get jphi for current region
    if(exist('efit_gfile_grid','var'))
       filename = efit_gfile_grid;
       if(~exist(filename,'file'))
          wait('ERROR make_ctrlvecs: specified file "efit_gfile_grid" is not in matlab path');
          return;
       end
       disp(' reading gfile data from efit_gfile_grid ...');
       gfile_data = read_gfile_tok(filename,tokamak);
    elseif(~exist('gfile_data','var')) % otherwise, an equilibrium data structure must be provided
       wait('ERROR make_ctrlvecs: either gfile_data or efit_gfile_grid must be provided')
       return;
    end
    if(length(rg)~=length(gfile_data.rg))
       fprintf('ERROR make_ctrlvecs: mismatch in efit dimensions between %s and %s\n',objfile,filename);
       wait; return;
    end
    jphi = gfile_data.jphi;
    idxg=find(jphi(:)~=0);
  else
    idxr=find_near(rg,rrange1); idxr=(idxr(1):idxr(2))';
    idxz=find_near(zg,zrange1); idxz=(idxz(1):idxz(2))';
    tmp=zeros(nz,nr);  tmp(idxz,idxr)=ones(length(idxz),length(idxr));
    idxg = find(tmp==1);    %idx vector for selected sub-grid
  end
  ngx = length(idxg);
  next_figure
    iplteq=1; ipltflx=1;
    ipltlim=1;   %turn off, since good limiter in EFIT data only...
    plot_tok_geo(tok_data_struct,plot_options)
    hold on
    plot(rgg(idxg),zgg(idxg),'m.')
    title('Sub-grid Locations','FontSize',15)
    wait('Check that calculation grid (magenta) is what you want for make_ctrlvecs.')
    print -dpsc ctrlvecs_plots.ps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build basic control vectors, make test data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  grmap = gbz2c(idxg,idxcr);  %const Bz over spec grid for Rctrl
  gzmap = gbr2c(idxg,idxcz);  %const Br over spec grid for Zctrl
  gipmap = mpc(idxg,idxcip);  %const flux over spec grid for ohmic (Ipctrl) 
  gkrmap = gbr2c(idxg,idxck); %specified Br over grid for Kctrl
  gkzmap = gbz2c(idxg,idxck); %specified Bz over grid for Kctrl

%Control vectors:
 %Constant fields and flux (R,Z,Ip):
  onegvec = ones(length(idxg),1);
  Rctrlvec1 = pinv(grmap)*onegvec;    %length(idxcr)
  Zctrlvec1 = pinv(gzmap)*onegvec;
  if(iconstraint)	% constrain to make all elements the same sign
     %Ipctrlvec1 = lsqlin(gipmap,onegvec,-eye(length(idxcip)),zeros(length(idxcip),1));
     Ipctrlvec1 = lsqlin(gipmap,onegvec,-eye(length(idxcip)),0*ones(length(idxcip),1));
  else
     Ipctrlvec1 = pinv(gipmap)*onegvec;
  end

 %Quadrupole field (constant Br,Bz gradient = 1 T/m):
  bzofrg = -1*(rg - rrefk);  %dist of Bz across rg (0 @R=rrefk)
  brofzg = 1*(zg - zrefk);  %dist of Br across zg (0 @Z=zrefk)
  bzgg = ones(nz,1)*bzofrg'; bzgg=bzgg(idxg); %quadrupole bz over idxg
  brgg = brofzg*ones(1,nr); brgg=brgg(idxg); %quadrupole br over idxg
  Kctrlvec1 = pinv(gkzmap)*bzgg;  % + 0*pinv(gkrmap)*brgg; 

Ipctrlvec1 = -Ipctrlvec1;
if(sign_ip == 1)
   Zctrlvec1 = -Zctrlvec1;
elseif(sign_ip == -1)
   Rctrlvec1 = -Rctrlvec1;
else
   wait('WARNING make_ctrlvecs: sign_ip is not specified => control vec response sign not defined')
end

% NOT SURE about sign of Kctrlvec1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test control vectors and display fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate fields:
  ccvec=Rctrlvec1;
  psir = reshape(mpc(:,idxcr)*ccvec,nz,nr);    %flux for Rctrl
  brctrl = reshape(gbz2c(:,idxcr)*ccvec,nz,nr);    %field for Rctrl

  ccvec=Zctrlvec1;
  psiz = reshape(mpc(:,idxcz)*ccvec,nz,nr);    %flux for Zctrl
  bzctrl = reshape(gbr2c(:,idxcz)*ccvec,nz,nr);    %field for Zctrl

  ccvec=Ipctrlvec1;
  psiip = reshape(mpc(:,idxcip)*ccvec,nz,nr);    %flux for Ipctrl
  tmp = sqrt((gbr2c(:,idxcip)*ccvec).^2 + (gbz2c(:,idxcip)*ccvec).^2);
  bipctrl = reshape(tmp,nz,nr);    %field for Ipctrl

  ccvec=Kctrlvec1;
  psik = reshape(mpc(:,idxck)*ccvec,nz,nr);    %flux for Kctrl
  bkctrl = reshape(gbr2c(:,idxck)*ccvec,nz,nr);    %field for Zctrl


if(~isnumeric(plot_options))

%Plot fields:

%R,Z: 
iplteq=0; ipltflx=0;   %turn off equilib plotting
ipltlim=1; nlevs=30;   %plot limiter to see plasma region
next_figure
if(~separate_figures)
 nrows=3; ncols=3;
 subplot(nrows,ncols,1)
end
    plot_tok_geo(tok_data_struct,plot_options)
    hold on
    contour(rg,zg,psir,nlevs)
    xlabel(' ')
    title('Rctrlvec Flux','FontSize',14)
if(separate_figures)
    print -dpsc -append ctrlvecs_plots.ps
 next_figure
else
 subplot(nrows,ncols,2)
end
    plot_tok_geo(tok_data_struct,plot_options)
    hold on
    contour(rg,zg,psiz,nlevs)
    xlabel(' ')
    title('Zctrlvec Flux','FontSize',14)
if(separate_figures)
    print -dpsc -append ctrlvecs_plots.ps
 next_figure
else
 subplot(nrows,ncols,3)
end
    plot_tok_geo(tok_data_struct,plot_options)
    hold on
    contour(rg,zg,psiip,nlevs)
    xlabel(' ')
    title('Ipctrlvec Flux','FontSize',14)
if(separate_figures)
    print -dpsc -append ctrlvecs_plots.ps
 next_figure
else
 subplot(nrows,ncols,4)
end
    mesh(rg,zg,psir)
    title('Rctrlvec Flux','FontSize',14)
if(separate_figures)
    print -dpsc -append ctrlvecs_plots.ps
 next_figure
else
 subplot(nrows,ncols,5)
end
    mesh(rg,zg,psiz)
    title('Zctrlvec Flux','FontSize',14)
if(separate_figures)
    print -dpsc -append ctrlvecs_plots.ps
 next_figure
else
 subplot(nrows,ncols,6)
end
    mesh(rg,zg,psiip)
    title('Ipctrlvec Flux','FontSize',14)
if(separate_figures)
    print -dpsc -append ctrlvecs_plots.ps
 next_figure
else
 subplot(nrows,ncols,7)
end
    mesh(rg,zg,brctrl)
    title('Rctrlvec Field (Bz)','FontSize',14)
if(separate_figures)
    print -dpsc -append ctrlvecs_plots.ps
 next_figure
else
 subplot(nrows,ncols,8)
end
    mesh(rg,zg,bzctrl)
    title('Zctrlvec Field (Br)','FontSize',14)
if(separate_figures)
    print -dpsc -append ctrlvecs_plots.ps
 next_figure
else
 subplot(nrows,ncols,9)
end
    mesh(rg,zg,bipctrl)
    title('Ipctrlvec Field','FontSize',14)
    print -dpsc -append ctrlvecs_plots.ps

%Ip, Kappa: 
iplteq=0; ipltflx=0;   %turn off equilib plotting
ipltlim=1; nlevs=30;   %plot limiter to see plasma region
next_figure
if(~separate_figures)
 nrows=3; ncols=2;
 subplot(nrows,ncols,1)
end
    plot_tok_geo(tok_data_struct,plot_options)
    hold on
    contour(rg,zg,psiip,nlevs)
    xlabel(' ')
    title('Ipctrlvec Flux','FontSize',14)
if(separate_figures)
    print -dpsc -append ctrlvecs_plots.ps
 next_figure
else
 subplot(nrows,ncols,2)
end
    plot_tok_geo(tok_data_struct,plot_options)
    hold on
    contour(rg,zg,psik,nlevs)
    xlabel(' ')
    title('Kctrlvec Flux','FontSize',14)
if(separate_figures)
    print -dpsc -append ctrlvecs_plots.ps
 next_figure
else
 subplot(nrows,ncols,3)
end
    mesh(rg,zg,psiip)
    title('Ipctrlvec Flux','FontSize',14)
if(separate_figures)
    print -dpsc -append ctrlvecs_plots.ps
 next_figure
else
 subplot(nrows,ncols,4)
end
    mesh(rg,zg,psik)
    title('Kctrlvec Flux','FontSize',14)
if(separate_figures)
    print -dpsc -append ctrlvecs_plots.ps
 next_figure
else
 subplot(nrows,ncols,5)
end
    mesh(rg,zg,bipctrl)
    title('Ipctrlvec Field','FontSize',14)
if(separate_figures)
    print -dpsc -append ctrlvecs_plots.ps
 next_figure
else
 subplot(nrows,ncols,6)
end
    mesh(rg,zg,bkctrl)
    title('Kctrlvec Field','FontSize',14)
    print -dpsc -append ctrlvecs_plots.ps

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final control vector assembly and reporting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Expand to nc-length vector:
  Rctrlvec=zeros(nc,1); Rctrlvec(idxcr)=Rctrlvec1;
  Zctrlvec=zeros(nc,1); Zctrlvec(idxcz)=Zctrlvec1;
  Ipctrlvec=zeros(nc,1); Ipctrlvec(idxcip)=Ipctrlvec1;
  Kctrlvec=zeros(nc,1); Kctrlvec(idxck)=Kctrlvec1;

%Normalize?:
  Rctrlvec = Rctrlvec;   %/norm(Rctrlvec);
  Zctrlvec = Zctrlvec;   %/norm(Zctrlvec);
  Ipctrlvec = Ipctrlvec; %/norm(Ipctrlvec);
  Kctrlvec = Kctrlvec;   %/norm(Kctrlvec);

%Define present (current based) predictor vecs:
  Rctrlvec_cur = Rctrlvec;
  Zctrlvec_cur = Zctrlvec;
  Ipctrlvec_cur = Ipctrlvec;
  Kctrlvec_cur = Kctrlvec;

%Report control vectors:
disp('Control vectors are: Rctrlvec, Zctrlvec, Ipctrlvec, Kctrlvec')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   disp('All done.')
               
