function [ivv,t_ivv] ...
          = calc_vv_eddy(shotnum,tok_data_struct,t,icoils,options,plasma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)calc_vv_eddy.m	1.3 04/21/08
%
%  SYNTAX:  [ivv,t_ivv] ...
%          = calc_vv_eddy(shotnum,tok_data_struct,t,icoils,options,plasma)
%
%  PURPOSE: Calculates VV eddy currents from PF coil ramps, with plasma 
%           effects optionally included.
%
%  INPUT: [default]
%	   shotnum         = shot number
%     tok_data_struct = TokSys vacuum object
%     t               = time vector for current waveforms
%     icoils          = coil current matrix (A) (dim = [ntimes,ncoils] 
%     options           
%        cutoff       = Cutoff frequency for filtering I_dot (Hz) [100Hz]
%        cor_offdrf   = dI/dt data correction (0=none 1=offset 2=drift) [0]
%        idxvv        = indices of VV elements to use in calculation [1:nv]
%        inc_pl       = include Ip effects (1=yes 0=no) [0]
%                          (requires plasma structure if 1)
%        ivv0         = initial vv currents (size = length(idxvv)) (A) [0]
%        n_offdrf     = # points to use for offset or drift correction [100]
%        pause_it     = pause after plots for inspection (1=yes 0=no) [0]      
%        ploted       = plot eddy current analysis (1=yes 0=no) [1]
%        resv_fac     = multiplication factor for VV resistance [1]
%        t0           = starting time for VV current integration (s) [0]
%     plasma
%        ip           = plasma current vector (A)
%        r_pl         = plasma major radius (m)
%        a_pl         = plasma minor radius (m)
%        k_pl         = plasma elongations [1.0]
%        d_pl         = plasma triangularity [0]
%        n_pl         = # pts defining plasma edge [100]
%        z_pl         = plasma vertical positions (m) [0]
%
%  OUTPUT:
%       ivv           = VV eddy current vectors (size = [nt,nv]);
%       t_ivv         = VV eddy current time vector (size = nt)
% 
%  RESTRICTIONS: 
%     1. Assumes that all data is uniformly sampled & in same time intervals
% 
%  METHOD:  
%
%  WRITTEN BY:  Jim Leuer 	ON 	9/28/06     taken from calc_eddy_flux
%
%  MODIFICATIONS:
%     2008-03-03  NWE   generalized to use with any tokamak, removed data
%                       retrieval
% =============================================================================

% Prelims
   struct_to_ws(tok_data_struct);
   struct_to_ws(options);

   if ~exist('cutoff')     cutoff     = 100;    end % 100Hz
   if ~exist('cor_offdrf') cor_offdrf = 0;      end
   if ~exist('idxvv')      idxvv      = 1:nv;   end
   if ~exist('inc_pl')     inc_pl     = 0;      end
   if ~exist('ivv0')       ivv0       = zeros(length(idxvv),1); end
   if ~exist('n_offdrf')   n_offdrf   = 100;    end % # pts to calc drift/offset
   if ~exist('pause_it')   pause_it   = 0;      end
   if ~exist('ploted')     ploted     = 1;      end
   if ~exist('resv_fac')   resv_fac   = 1.0;    end
   if ~exist('testit')     testit     = 0;      end
   if ~exist('t0')         t0         = 0;      end

   if inc_pl == 1 
      if ~exist('plasma')
         error('Plasma structure not specified.')
      else
         struct_to_ws(plasma);
         if ~exist('k_pl')    k_pl  = 1.0;      end
         if ~exist('d_pl')    d_pl  = 0.0;      end
         if ~exist('n_pl')    n_pl  = 100;      end
         if ~exist('z_pl')    z_pl  = 0.0;      end  
         [rpl,zpl]=  dee(r_pl,z_pl,a_pl,k_pl,d_pl,(2*pi*(0:n_pl-1)/n_pl)');
         rpl(end+1) = rpl(1); % close the contour
         zpl(end+1) = zpl(1);
      end
   end   
   nv = length(idxvv);   
% End Prelims 

 
% Generate influence matrix from plasma to coils and VV  
   if inc_pl     
      id_pl= find(inpolygon(rgg,zgg,rpl,zpl));
      % here is where you would put in a current distribtion - currently uniform:
      npl= length(id_pl);
      mvpp= sum(mpv(id_pl,idxvv))'/npl; % mutual from plasma to VV elements
      if testit 
         inside_pl= zeros(size(rgg));
         inside_pl(id_pl)= 1;
         spy(inside_pl);
         mcpp= sum(mpc(id_pl,:))'/npl;
         [idr,idz]= find(inside_pl);
         ipf_im= [6850.307143; 7120.175929; 8029.917857; 3260.438409;
	               3260.438431; 1283.91216; 748.403125];
         ipf_im= [ipf_im; ipf_im; 0; 0;];
	      sum(ipf_im.*mcpp); % this should be 3vs 
      end
      ip_dot = diff(ip)./diff(t);
      ip_dot(end+1)=ip_dot(end);
   end
% End Generate Plasma Influence Matrix


% Integrate the VV currents using the PF currents
   % Mvv*Ivv_dot + R*Ivv = -Mvc*Icc_dot => Ivv_dot=Mvv^-1*R*Ivv-Mvv^-1*Mvc*Icc_dot
   M= mvv(idxvv,idxvv);
   R= diag(resv(idxvv))*resv_fac;
   Minv= inv(M);
   if inc_pl
      amat= -Minv*R;
      bmat= -Minv*[mcv(:,idxvv)' mvpp];
      cmat= eye(nv);
      dmat= zeros(nv,nc+1);
   else
      amat= -Minv*R;
      bmat= -Minv*mcv(:,idxvv)';
      cmat= eye(nv);
      dmat= zeros(nv,nc);           
   end
   vvsys= ss(amat, bmat, cmat, dmat);
 
   ic_dot = diff(icoils)./repmat(diff(t),1,nc);
   ic_dot(end+1,:) = ic_dot(end,:); %pad to keep same size
   for ii=1:size(ic_dot,2)
      if cor_offdrf == 1     % use beginning to subtract offset
         ic_dot(:,ii)= ic_dot(:,ii)  - mean(ic_dot(1:n_offdrf,ii));
      elseif cor_offdrf == 2        % remove
         ic_dot(:,ii)= correct_drift(ic_dot(:,ii),t,n_offdrf);
      else % no offset
      end
   end
        
   ic_dotf= [];
   if exist('cutoff')  % apply LP filter to dI/dt
      if cutoff > 0
         for jj=1:size(ic_dot,2)
            [ic_dotff,t]= fft_filter(ic_dot(:,jj), t, cutoff);
            ic_dotf= [ic_dotf ic_dotff];
         end
         if inc_pl [ip_dotf,t]= fft_filter(ip_dot, ip_t, cutoff); end
         else
            ic_dotf= ic_dot;
          if inc_pl ip_dotf= ip_dot; end
      end
   end

   [tss0, idxts]= min(abs(t-t0));
   tss= t(idxts:end);
   ic_dott= ic_dotf(idxts:end,:);
   if inc_pl ip_dott= ip_dotf(idxts:end); end
   if inc_pl
      [ivv,ttt]= lsim(vvsys, [ic_dott ip_dott], tss - tss(1),ivv0);
   else
      [ivv,ttt]= lsim(vvsys, ic_dott, tss - tss(1),ivv0);
   end
   tss= ttt+tss(1);  
   t_ivv= tss; % main time output   
% End Vessel Current Integration


% Plots
   if ploted
      if inc_pl
         next_figure,clf
         subplot_ga(2,1,1)
         title(['Plasma Current & Time derivitive ' int2str(shotnum)])
         plot(ip_t, ip*1e-6, 'r')
         ylabel('I_{pl} (MA)')
         rm_x_label;
         grid on
         subplot_ga(2,1,2)
         plot(ip_t, ip_dotf*1e-6, 'b')
         ylabel('I_{pl} (MA)')
         ylabel('I_{pl\_dot} (MA/s)')
         xlabel('Time (s)')
         grid on
         ax= axis;
         axis([ax(1:2) -4 4]);
      end
      if pause_it wait; end
      
      next_figure, clf
      [E,d] = eigsort(amat);
      subplot(3,1,1)
      plot(-1./d*1e+3)
      title(['Vacuum Vessel Eddy Currents ' int2str(shotnum)])
      ylabel('VV Time Cons. (1/ms)')
      xlabel(' Mode Number')
      grid on
      %
      subplot(3,1,2)
      hold off
      %[idmt,idmtx]= max(abs(ic_dotf)); % find maximum ramp rate
      id= 1;
      plot(t,ic_dot(:,id),'r')
      hold on
      plot(t, ic_dotf(:,id), 'b')
      ylabel('Ics1\_dot (A/s)')
      xlabel('Time (s)')
      grid on
      %id= 1;
      %
      subplot(3,1,3)
      ivv_eddy= sum(ivv')';
      plot(tss,ivv_eddy*1e-3)
      ylabel('Ivv\_edd (kA)')
      xlabel('Time (s)')
      grid on
      hold on
      if pause_it wait; end

     % display coil currents           
      rows = ceil(nc/2);
      columns = 2;
      next_figure, clf     
      for k = 1:nc
         subplot(rows,columns,k)
         plot(t, 1e-3*icoils(:,k))
         ylabel(fcnames(k,:))
         set(gca,'xtick',[])
         grid on
      end
      subplot(rows,columns,1)
         title([' PF Coil Currents (kA); Shot ' int2str(shotnum)]);
      subplot(rows,columns,nc-1)
         set(gca,'xtickMode', 'auto')
         xlabel('Time (s)')
      subplot(rows,columns,nc)
         set(gca,'xtickMode', 'auto')
         xlabel('Time (s)')
      if pause_it wait; end
            
      % display dI/dt
      next_figure, clf 
      for k = 1:nc
         subplot(rows,columns,k)
         plot(t, 1e-3*ic_dotf(:,k))
         ylabel(fcnames(k,:))
         set(gca,'xtick',[])
         grid on
      end
      subplot(rows,columns,1)
         title([' PF Coil Ramp Rate (kA/s); Shot ' int2str(shotnum)]);
      subplot(rows,columns,nc-1)
         set(gca,'xtickMode', 'auto')
         xlabel('Time (s)')
      subplot(rows,columns,nc)
         set(gca,'xtickMode', 'auto')
         xlabel('Time (s)')
      if pause_it wait; end
   end 
      

