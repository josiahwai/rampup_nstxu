function [flux,B,Br,Bz,ivv,t_ivv,imc]= calc_eddy_flux(shotnum,tok_data_struct,rtimes, ...
                                      t,Icoils,psidata,probedata,options,plasma)
 %
%  SYNTAX: [flux,B,Br,Bz,ivv,t_ivv,imc]= 
%                  calc_eddy_flux(shotnum,tok_data_struct,rtimes,t,Icoils,...
%                                   psidata,probedata,options,plasma)
%
%  PURPOSE: Compute flux/field in vacuum region for multiple rtimes using a
%	         fit to measured diagnostics, then compute approximate loop voltage.
%           Fit may optionally use eigenmode approximation of vessel or do full
%           integration of VV eddy currents for fitting vessel effects, include 
%           axisymmetric magnetic materials, and/or include Ip effects.
%
%  INPUT: [default]
%	   shotnum     = number of shot to reconstruct
%	   rtimes      = time steps at which to reconstruct vacuum flux/field
%	   t           = time vector for data signals
%     Icoils      = PF coil current data
%     psidata     = flux loop data
%     probedata   = magnetic probe data    
%     options    
%        cccirc   = PF coils circuit connection mapping
%        Pcc      = mapping from coil states to all PF coils
%                   (only one of cccirc or Pcc should be provided)
%        clevels  = # levels on contour plots [10]
%        idxbp    = indices of B-probes to fit [all=1:nbp]
%        idxfl    = indices of flux loops to fit [all=1:nfl]
%        idxvv    = indices of VV elements to use [all=1:nv]
%        ieddy    = 0: Use eigenmode expansion of VV to fit vv currents
%                   1: Use integrated VV eddy currents for fitting                 
%        navg     = # time points over which measurements are averaged [3]
%        neig     = # eigenmods to reating if ieddy=0 [10]
%        pause_it = pause after plots for inspection (1=yes 0=no) [0] 
%        plotb    = plot magnetic fields (1=yes 0=no) [1]
%        plotfl   = plot magnetic flux (1=yes 0=no) [1]
%        (SEE calc_vv_eddy for further options)
%     plasma   (see calc_vv_eddy)
%
%  NOTE: To reproduce results of calc_vac_flux simply put ieddy=0.
%
%  OUTPUT:
%     flux     = magnetic flux at plasma grid points for each time step	
%     B,Br,Bz  = magnetic field """" 
%     ivv      = calculated vessel currents (if ieddy = 1)
%     t_ivv    = time vector for calculated vessel currents
%     imc      = magnetization currents
%     
%     Plots:
%        1. contour plots of magnetic flux and field over the kstar x-section
%
%  RESTRICTIONS:
%     1. Assumes that all data is uniformly sampled & in same time intervals
%  
%  WRITTEN BY:  Jim Leuer 	ON 	9/22/06
%     Optimization part taken from M. Walkers calc_vac_flux and modified.
%     also different plot features
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @(#)calc_eddy_flux.m	1.7 02/02/10

% Prelims
   struct_to_ws(options);
   struct_to_ws(tok_data_struct);
   
   % Set default optons if not specified
   iblackbg= 0;
   if(exist('cccirc')==1 & exist('Pcc')==1)
      wait('ERROR calc_eddy_flux: only one of cccirc or Pcc may be specified')
      return;
   end
   if(exist('cccirc')~=1 & exist('Pcc')~=1)  cccirc   = 1:nc;           end
   if ~exist('clevels') clevels  = 15;             end
   if ~exist('ieddy')   ieddy    = 0;              end
   if ~exist('navg')    navg     = 3;              end
   if ~exist('neig')    neig     = 10;             end
   if ~exist('pause_it')pause_it = 0;              end
   if ~exist('plasma')  plasma   = [];             end
   if ~exist('plotb')   plotb    = 1;              end 
   if ~exist('plotfl')  plotfl   = 1;              end
   if ~exist('fitcc') 	fitcc	   = 0;			      end
   if ~exist('fitvv')   fitvv    = 1;              end
   if ~exist('doplots') doplots  = 1;              end
   if ~exist('lockmag') lockmag  = 1;              end  
   
   % Setup matrices depending on fitting options
   nv = size(vvdata,2);
   if ~exist('nmc') | nmc == 0
	   nmc = 0;
      mlm = [];
      mmv = [];
      gbm = [];
      imc = [];
   end
   if ~fitvv
      mlv = [];
      mvv = [];      
      gbv = [];
      mpv = [];
      gbz2v = [];
      gbr2v = [];
      nv = 0;
   end
   mlm1 = mlm;
   mmv1 = mmv;
   gbm1 = gbm;
   nmc1 = nmc;  
   if ~ieddy | ~fitvv
      ivv = [];
      t_ivv = [];
   end
      
   fnavg2 = floor(navg/2); 
   if(exist('cccirc')==1)
      coil_constraints = make_constraints(cccirc);  
   else
      coil_constraints = make_constraints(Pcc);  
   end
   ntimes = length(rtimes);
   figure_num= 0;
   fig_start= gcf;
   
   % Below needed for quiver....  Should generalize...
   rp  = 1.8;
   drp = 0.5;
   idpr=  interp1(rg,1:length(rg),[rp, rp-drp/2 , rp+drp/2],'nearest');
   idpz=  interp1(zg,1:length(zg),zeros(size(idpr)),'nearest');
   %
% End Prelims


   % Expand into eigenvectors if not using full vessel model   
   if ~ieddy & fitvv 
      aa=-inv(mvv)*diag(resv);
      [V,d]= eigsort(aa);
      mlv   = mlv*V(:,1:neig);
      mmv   = mmv*V(:,1:neig);
      mpv   = mpv*V(:,1:neig);
      gbv   = gbv*V(:,1:neig);
      gbr2v = gbr2v*V(:,1:neig);
      gbz2v = gbz2v*V(:,1:neig);
      nv = neig;
   end


% Start rtimes Loop
   for i=1:ntimes
      [mm,idxt] = min(abs(rtimes(i)-t));
     
      Icoilss = mean(Icoils(idxt-fnavg2:idxt+fnavg2,:))';
      Imean  = mean(abs(Icoilss));
      if ~isempty(psidata)
         psidataa = mean(psidata(idxt-fnavg2:idxt+fnavg2,:))';
         flmean = mean(abs(psidataa));
         scale_fl = Imean/flmean*2;
      else
         psidataa= [];
         scale_fl = 0;
      end
      if ~isempty(probedata)
         probedataa = mean(probedata(idxt-fnavg2:idxt+fnavg2,:))';
         bpmean = mean(abs(probedataa));
         scale_bp = Imean/bpmean*2;
      else
         probedataa= [];
         scale_bp = 0;
      end         
      if lockmag & i > 1 % magnetization already locke in - do not fit 
         mlm1 = [];
         gbm1 = [];
         nmc1 = 0;
      end 
      
      if ~ieddy % eddy currents not calculated, fit using eigenmodes of VV
         if fitcc
            b = [scale_fl*psidataa; scale_bp*probedataa; Icoilss; ...
                  zeros(size(coil_constraints,1),1)];
            A = [scale_fl*[mlc mlv mlm1]; ...
                 scale_bp*[gbc gbv gbm1]; ...
                 [eye(nc) zeros(nc,neig+nmc1)]; ...
                 [coil_constraints zeros(size(coil_constraints,1),neig+nmc1)]];
         else
            b = [scale_fl*psidataa; scale_bp*probedataa];
            A = [scale_fl*[mlv mlm]; ...
                 scale_bp*[gbv gbm]];
         end
         nvv = 0;
      else % eddy currents precalculated and used in fitting   
         if fitvv
            [ivv,t_ivv] = calc_vv_eddy(shotnum,tok_data_struct,t,Icoils,options,plasma);       
            ivvv =mean(ivv(idxt-navg:idxt+navg,:))';
         end
         if ~fitvv
            ivvv = [];
            ivvmean = 0;
            scale_vv = 0;
         else
            ivvmean= mean(ivvv);
            scale_vv = Imean/ivvmean;
	      end
         if fitcc
            b = [scale_fl*psidataa; scale_bp*probedataa; Icoilss; scale_vv*ivvv; ...
                  zeros(size(coil_constraints,1),1)];
            A = [scale_fl*[mlc mlv mlm1]; ...
                 scale_bp*[gbc gbv gbm1]; ...
                 [eye(nc) zeros(nc,nv+nmc1)]; ...
                 scale_vv*[zeros(nv,nc) eye(nv) zeros(nv,nmc1)]; ...
                 [coil_constraints zeros(size(coil_constraints,1),nv+nmc1)]];
         else
            b = [scale_fl*psidataa; scale_bp*probedataa;scale_vv*ivvv];
            A = [scale_fl*[mlv mlm1]; ...
                 scale_bp*[gbv gbm1]; ...
                 scale_vv*[eye(nv) zeros(nv,nmc1)]];
         end
         % (Last 2 rows are constraints on connected coils.)
         nvv = nv;
      end % if
      if ~fitvv
         nvv = 0;
      end
      if ~fitcc
         b = b - [scale_fl*mlc*Icoilss; scale_bp*gbc*Icoilss;zeros(nvv,1)];
      end
      
         
      % fit currents & calulate flux & fields
      I = pinv(A)*b;
      if ~fitcc
         I = [Icoilss; I];
      end
      if lockmag & i>1 & nmc>0
         I = [I; imc(:,1)]
      end
      psi = [mpc mpv mpm]*I;
      flux(:,:,i) = reshape(psi,nr,nz);
      Br_tmp = [gbr2c gbr2v gbr2m]*I;
      Bz_tmp = [gbz2c gbz2v gbz2m]*I;
      Br(:,:,i) = reshape(Br_tmp,nr,nz);
      Bz(:,:,i) = reshape(Bz_tmp,nr,nz);
      B(:,:,i) = sqrt(Br(:,:,i).^2+Bz(:,:,i).^2);
      imc(:,i) = I(end-nmc+1:end);

      if doplots
         % time slice plots
         figure
         plot(b+[scale_fl*mlc*Icoilss; scale_bp*gbc*Icoilss;zeros(nvv,1)])
         hold on
         grid on
         if ~fitcc & ieddy
            A = [[mlc mlv mlm]*scale_fl; ...  
                 [gbc gbv gbm]*scale_bp; ... 
                 [zeros(nv,nc) eye(nv) zeros(nv,nmc)]*scale_vv; ...
                 [coil_constraints zeros(size(coil_constraints,1),nv+nmc)]];
         elseif ~fitcc & ~ieddy
            A = [[mlc mlv mlm]*scale_fl; ... 
                 [gbc gbv gbm]*scale_bp; ...
                 [coil_constraints zeros(size(coil_constraints,1),nv+nmc)]]; 
         end                 
         plot(A*I,'r')
         xlabel('Data Index (Order: FL BP Ipf Ivv)')
         ylabel('Normalized Value (normalized by mean(Ipf)')
         title(['QUALITY OF FIT: data vs fit (RED=FIT BLUE=DATA), t=' num2str(rtimes(i))])
         if pause_it wait; end

         figure
         plot(I,'r')
         hold on
         if ieddy & fitvv
            plot([Icoilss; ivvv],'b--')
            title(['Ipf & Ivv Quality of Fit (RED=FIT BLUE=DATA) , t=' num2str(rtimes(i))])
            xlabel('Data Index (Order: Ipf Ivv)')
            ylabel('Normalized Value (normalized by mean(Ipf)')
         else
            plot(Icoilss,'b--')
            title(['Ipf Quality of Fit (RED=DATA BLUE=FIT) , t=' num2str(rtimes(i))])
            xlabel('Data Index (Order: Ipf Ivv_modes)')
            ylabel('Normalized Value (normalized by mean(Ipf)')
         end
         grid on
         hold off
         if pause_it wait; end

         figure
         if fitvv & ieddy
            rvv= vvdata(2,:)';
            zvv= vvdata(1,:)';
            ivv_fit= I(nc+1:nc+nv);
            r0= mean(rvv);
            z0= mean(zvv);
            avv= 180/pi*atan2(zvv-z0,rvv-r0);
            plot(avv,ivv_fit,'rx')
            hold on
            grid on
            title(['VV Eddy Currents Vs Angle (Red x = fit, Blue o = model), t=' ...
               num2str(rtimes(i))])
            xlabel(' Poloidal Angle from outer midplane, (Deg)')
            ylabel(' Eddy Current, (A) ')
            i_eddy= sum(I(nc:nc+nv));
            text(.1,.01,['Total fitted VV Eddy Currents = ' num2str(i_eddy)], ...
               'Units','Normalized','color','r')
            if exist('ivvv')==1
               plot(avv,ivvv,'bo')
               text(.1,.05,['Total Model Based VV Eddy Currents = ' num2str(sum(ivvv))], ...
                  'Units','Normalized','color','b')
            end
         elseif fitvv
            Iv = V(:,1:neig)*I(nc+nmc+1:end);
            plot(Iv)
            title('current in vessel elements')
            ylabel('Amps')
            xlabel('element number')
            fprintf('sum of currents in vessel = %f\n',sum(Iv));
         end
         if pause_it wait; end


         if plotfl
            figure
            plot_tok_geo(tok_data_struct)
            title(['vacuum flux, shot=' int2str(shotnum) ',t=' num2str(rtimes(i))])
            hold on
            [c,h]=contour(rg,zg,flux(:,:,i),clevels);
            clabel(c,h)
            axis([1,3,-1.5,1.5])
            if pause_it wait; end
         end % plotfl

         if plotb
            figure
            subplot(1,3,1)
            ipltfl=0; ipltbp=1;
            plot_tok_geo(tok_data_struct)
            title(['vacuum B(G), t=' num2str(rtimes(i))])
            hold on
            [c,h]=contour(rg,zg,B(:,:,i)*1e4,clevels);
            clabel(c,h)
            quiver(rg(idpr), zg(idpz), Br(idpz,idpr,i),Bz(idpz,idpr,i) ,'r')
            dum= diag(B(idpz,idpr,i));
            text(rg(idpr), zg(idpz), int2str(dum*10000),'Color','r')
            %axis([1,3,-1.5,1.5])
            subplot(1,3,2)
            plot_tok_geo(tok_data_struct)
            title(['vacuum Bz(G), t=' num2str(rtimes(i))])
            hold on
            [c,h]=contour(rg,zg,Bz(:,:,i)*1e4,clevels);
            clabel(c,h)
            subplot(1,3,3)
            ipltfl=0; ipltbp=1;
            plot_tok_geo(tok_data_struct)
            title(['vacuum Br(G), t=' num2str(rtimes(i))])
            hold on
            [c,h]=contour(rg,zg,Br(:,:,i)*1e4,clevels);
            clabel(c,h)
            if pause_it wait; end
         end % if plotb
      end % if doplots
   end  %rtimes loop

if 0
   if ntimes > 1
      % plot loop voltage
      figure
      for i=2:ntimes
         psidot = (flux(:,:,i) - flux(:,:,i-1))/(rtimes(i)-rtimes(i-1));
         ipltfl=1; ipltbp=0;
         plot_tok_geo(tok_data_struct)
         hold on
         title(['Vacuum Loop Voltage, t=' num2str(0.5*(rtimes(i)+rtimes(i-1)))])
         [c,h]=contour(rg,zg,-psidot,clevels);
         clabel(c,h);
         axis([1,3,-1.5,1.5]);
         hold off;
      end
   end
end
% End calc_eddy_flux



function constraints = make_constraints(input)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% SYNTAX: constraints = make_constraints(input)
%
% PURPOSE: Converts cccirc or Pcc to independent circuit connection contraints.
%
% INPUT:
%   input = either cccirc or Pcc
%
% OUTPUT:
%   constraints = matrix whose rows define constraints on PF currents (constraints * Ipf = 0)
%
% RESTRICTIONS
%     Only handles pairs of coils if cccirc is input. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(min(size(input))==1)		% cccirc input 
   constraints = [];
   for k = setdiff(unique(abs(cccirc)),0)
      dups = find(abs(cccirc) == k);
      if length(dups) >= 2  % coils tied together
         constraint = zeros(1,length(cccirc));
         constraint(dups(1)) = 1;
         for ii = 2:length(dups)
            constraint(dups(ii)) = -( sign(cccirc(dups(1)))*...
                                       sign(cccirc(dups(ii))));  
         end        
         constraints = [constraints; constraint];
      end
   end
else
   constraints = null(Pcc')';  
end

% End constraints
