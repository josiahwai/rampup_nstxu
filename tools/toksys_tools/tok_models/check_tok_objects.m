   function [valsv] = check_tok_objects(tok_data_struct,plot_geo_fun, ...
				efit_gfile,read_gfile_fun)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  USAGE:  >> check_tok_objects(tok_data_struct,plot_geo_fun, ...
%				efit_gfile,read_gfile_fun);
%
%  PURPOSE: Function to do checks of mutuals, Green functions for 
%	tokamak system analogous to D3D Electromagnetic Environment. 
%
%  INPUTS:
%	tok_data_struct = standard structure containing vacuum objects
%	plot_geo_fun = name of geometry plotting script to use 
%			(e.g. 'plot_east_geo')
%	efit_gfile = (opt) name of EFIT gfile for equilibrium plotting
%	read_gfile_fun = (opt) name of gfile reader to use (if efit_gfile)
%
%  OUTPUTS:
%
%    Checks of various objects, plots to confirm...
%	valsv = VV eigenvalues (rad/sec)
%
%  RESTRICTIONS:
%
%
%  METHOD:  
%	cmpares some object values with grid-to-s.t., 
%	does some eigenvalue calculations, displays results to check...
%

%  WRITTEN BY:  Dave Humphreys  ON	3/22/06
%
%  MODIFICATION HISTORY:
%	
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Explicit Inputs:

% Prelims and Constants:
   mu0 = 0.4*pi;
   twopi = 2*pi;
   mcc = 0;    %clear function def of mcc...

% Derived Values:
  struct_to_ws(tok_data_struct);   %unpack structure...
                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Confirmation plot:
  figure(1),clf
  iplteq=0; ipltflx=0; ipltfl=0; ipltbp=0;
  ilabelfc=0; ilabelvv=0; ilabelfl=0; ilabelbp=0;
  if exist('psizr'), iplteq=1; ipltflx=1; end
  if ~isempty(fldata)
    ipltfl=1;   %turn on plotting of FL's
    ilabelfl=0; %turn off labeling of FL indices
  end
  if ~isempty(bpdata)
    ipltbp=1;   %turn on plotting of BP's
    ilabelbp=0; %turn off labeling of BP indices  
  end
  if exist(efit_gfile)
	filename = efit_gfile;
	eval(read_gfile_fun)
        iplteq = 1;
	ipltflx = 1;
  end

  eval(plot_geo_fun)

  %wait('Paused: CR to continue...')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eigenvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VV:
   rvv = diag(resv);
   Avv = -inv(mvv)*rvv;
   [vecs,vals] = eigsort(Avv);   
   valsv = vals;
   disp('Vessel eigenvalues (rad/sec, sec):')   
   [vals(1:10) 1./vals(1:10)]		%print top 10 VV eigenvalues
   figure(2),clf,hold off
   tmp = vvdata';
   subplot(1,2,1)
    mode = vecs(:,1);
    plot_mode(tmp(:,2),tmp(:,1),tmp(:,4),tmp(:,3),mode,tmp(:,5),tmp(:,6));
    title('VV Mode #1')
    xlabel('R [m]')
    ylabel('Z [m]')
    axis image
   subplot(1,2,2)
    mode = vecs(:,2);
    plot_mode(tmp(:,2),tmp(:,1),tmp(:,4),tmp(:,3),mode,tmp(:,5),tmp(:,6));
    title('VV Mode #2')
    xlabel('R [m]')
    ylabel('Z [m]')
    axis image
  % wait('Paused: CR to continue...')

%PF's:
   rcc = diag(resc);
   Acc = -inv(mcc)*rcc;    
   [vecs,vals] = eigsort(Acc);      
   disp('PF eigenvalues (rad/sec,  sec):') 
   [vals 1./vals]			%print all PF eigenvalues
   figure(3),clf,hold off
   tmp = fcdata';
   subplot(1,2,1)
    mode = vecs(:,1);
    plot_mode(tmp(:,2),tmp(:,1),tmp(:,4),tmp(:,3),mode,tmp(:,5),tmp(:,6));
    title('PF Mode #1')
    xlabel('R [m]')
    ylabel('Z [m]')
    axis image
   subplot(1,2,2)
    mode = vecs(:,2);
    plot_mode(tmp(:,2),tmp(:,1),tmp(:,4),tmp(:,3),mode,tmp(:,5),tmp(:,6));
    title('PF Mode #2')
    xlabel('R [m]')
    ylabel('Z [m]')   
    axis image
    %wait('Paused: CR to continue...') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vacuum field patterns 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  idx1 = 1;
  vec=zeros(nc,1); vec(idx1)=1; psit1=reshape(mpc*vec,nz,nr);
  idx2 = fix(nc/2);
  vec=zeros(nc,1); vec(idx2)=1; psit2=reshape(mpc*vec,nz,nr);
  
  figure(4),clf,hold off
  subplot(1,2,1)
    eval(plot_geo_fun)
    hold on
    contour(rg,zg,psit1,30,'m')
    plot(fcdata(2,idx1),fcdata(1,idx1),'rx','MarkerSize',14,'LineWidth',3);
  title('Flux Contours for PFs','FontSize',15)
  subplot(1,2,2)
    eval(plot_geo_fun)
    hold on
    contour(rg,zg,psit2,30,'m')
    plot(fcdata(2,idx2),fcdata(1,idx2),'rx','MarkerSize',14,'LineWidth',3);

  idx1=1;
  vec=zeros(nv,1); vec(idx1)=1; psit1=reshape(mpv*vec,nz,nr);
  idx2 = fix(nv/2);
  vec=zeros(nv,1); vec(idx2)=1; psit2=reshape(mpv*vec,nz,nr);
  
  figure(5),clf,hold off
  subplot(1,2,1)
    eval(plot_geo_fun)
    hold on
    contour(rg,zg,psit1,30,'m')
    plot(vvdata(2,idx1),vvdata(1,idx1),'gx','MarkerSize',14,'LineWidth',3);
  title('Flux Contours for VVs','FontSize',15)
  subplot(1,2,2)
    eval(plot_geo_fun)
    hold on
    contour(rg,zg,psit2,30,'m')
    plot(vvdata(2,idx2),vvdata(1,idx2),'gx','MarkerSize',14,'LineWidth',3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   disp('All done.')
               
              
