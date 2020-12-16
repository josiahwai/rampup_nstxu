function predictors = make_predictors(predictor_inputs)
 %
%  USAGE:  predictors = make_predictors(predictor_inputs)
%
%  PURPOSE: Script to construct plasma geometry predictors from standard
%		environment inputs. Intended to be generic to any standard
%		tokamak system environment.Called from scripts such as
%		make_east_predictors or make_kstar_predictors.
%
%  INPUTS:  (in structure "predictor_inputs")
%	tok_data_struct = vacuum geometry objects structure for the device
%       efit_gfile_test = gfile to use (for testing of predictors only)
%       tokamak = machine name, eg 'KSTAR', 'EAST', (needed for read_gfile_tok)
%       idxbpr = index vector to select desired probes for R-predictors
%       idxbpz = index vector to select desired probes for Z-predictors
%       idxbpip = index vector to select desired probes for Ip-predictors
%       idxbpk = index vector to select desired probes for kappa-predictors
%       idxflr = index vector to select desired flux loops for R-predictors
%       idxflz = index vector to select desired flux loops for Z-predictors
%       idxflip = index vector to select desired flux loops for Ip-predictors
%       idxflk = index vector to select desired flux loops for kappa-predictors
%	idxc = PF coil index vector for fitting/testing with EFIT
%	idxv = VV index vector for testing with VV currents
%       rrange1 = maj radius range of sub-grid for predictor 1 (m) eg [1.5 2]
%       zrange1 = vert range of sub-grid for predictor 1 (m) eg [-.5 .5]
%       iusejphi = flag to use jphi (all current) region for fits
%       efit_gfile_grid = gfile to use to define grid from jphi if iusejphi=1
%       nprtst = # of princ comp of grid to mag set to use in test
%       nsvdz2 = # of sing vals to keep in R pred inverse (for P*pred1,2)
%       nsvdr2 = # of sing vals to keep in R pred inverse
%       nsvdi2 = # of sing vals to keep in Ip pred inverse
%	nsvdk1 = # of sing vals to keep in Kappa pred inverse
%	nrs
%	nzs
%	nkaps
%	nls
%	Rmin,Rmax
%	Zmin,Zmax
%	kapmin,kapmax
%	limin,limax
%	amax
%	Z1,Z2
%	R1,R2
%	cctest
%	icomppfr
%	icomppfz
%	icomppfip
%	icomppfk
%
%  OUTPUTS: (in structure "predictors")
%	Pzpred = Z-predictor (Pzpred*[mag]/Ip0 = predicted Z)
%	Prpred = R-predictor (Prpred*[mag]/Ip0 + Rrel = predicted R)
%	Pipred = Ip-predictor (Pipred*[Bprobes] = predicted Ip)
%	Pkpred = kappa-predictor (Pkpred*[mag]/Ip0??? = predicted kappa)
%	Rpred
%	Zpred
%	Ippred
%	Kpred
%	+ many  plots evaluating quality of predictors
%
%  RESTRICTIONS:
%	Specified objects file must provide units desired in predictor (eg
%	if want predictor based on mks, must use mks object file). Comments
%	in code assume MA-based objects (eg east_objects.mat), but should
%	be consistent with mks/terminal choice if use those object files). 
%     If iusejphi==1, must have efit_gfile_grid defined (which is source of
%       jphi on which grid is based if iusejphi=1)
%
%  METHOD:  
%	Variety of inversions of grid-to-mag objects, linear fits to selected
%	data from subgrids, subset of magnetics...

%  WRITTEN BY:  Dave Humphreys  ON	6/24/05
%
%  MODIFICATION HISTORY:
%	2/9/06 DAH Removing Pred1 completely: useless predictor based on 
%		estimation of plasma current distribution alone without
%		discriminating coil currents.
%	1/12/07 DAH Moved input variables for Kappa predictor into
%		make_*_predictors code; no longer hardwired here in the
%		make_predictors.m code... 
%
%  NOTES:
%	1/20/06 DAH Things to add still: ability to put in different mag idx
%		vectors for each predictor, plotting of selected magnetics
%		for each such new idx vector, (minvert DONE)
%		(probably good idea so can control # of sig vals, although
%		it's slower...), evaluation of effectiveness of magnetics...
%		Add PF coils (DONE) and VV to fit (discriminate conductor 
%		   currents)
%		Remove Pred1 that omits conductor current discrimination
%		Add ability to specify multiple rect regions to union in 
%		   order to make current grid set that contains all expected
%		   plasma region.
%		Also add predictors based on filaments for true shape preds...
%		Add tests of Pred2 with full magdata from sing vals descrip
%		     (maybe just replace the Pred1 tests)
%		Add statistical kappa pred algorithm (gen ensenble of
%		   deterministic cases, not random...)
%
%        6/6/06 DAH Add statistical kappa pred algorithm (gen ensenble of
%                  deterministic cases, not random...)
%	6/24/08 DAH Adding option for relative R pred with input of Rrel value
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @(#)make_predictors.m	1.5 10/13/10

%       ipltlim = 1 to plot limiter in plot_geo call (0 to not)
%       ipltbp = 1 to plot Bprobes in plot_geo call (0 to not)
%       ipltfl = 1 to plot flux loops in plot_geo call (0 to not)

   mcc = 0;	% KLUGE to get matlab to recognize VARIABLE mcc

   struct_to_ws(predictor_inputs);

% Inputs: NOW MUST BE DONE IN SCRIPT THAT CALLS make_predictors
   if exist('Rrel')~=1, Rrel=0; end   %default to Rpred relative to R=0
   
% Prelims and Constants:
   mu0 = 0.4*pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract objects, select sub-grid, select diagnostics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   struct_to_ws(tok_data_struct);   %unpack objects structure to workspace

%Select sub-grid to base predictor on:
  if iusejphi==1    %must have efit_gfile to get jphi for current region
    filename = efit_gfile_grid;
    gfile_data = read_gfile_tok(filename,tokamak);
    jphi = gfile_data.jphi;
    idxg=find(jphi(:)~=0);
  else
    idxr=find_near(rg,rrange1); idxr=(idxr(1):idxr(2))';
    idxz=find_near(zg,zrange1); idxz=(idxz(1):idxz(2))';
    tmp=zeros(nz,nr);  tmp(idxz,idxr)=ones(length(idxz),length(idxr));
    idxg = find(tmp==1);    %idx vector for selected sub-grid
  end
  ngx = length(idxg);
  figure(1),clf,hold off
    iplteq=1; ipltflx=1;
    ipltlim=0;   %turn off, since good limiter in EFIT data only...
    plot_tok_geo(tok_data_struct)
    hold on
    plot(rgg(idxg),zgg(idxg),'m.')
    title('Sub-grid Locations','FontSize',15)

%Select desired diagnostic set, construct grid-to-magnetics mapping:
%All mags for princ comp tests:
   gmp = [mpl(idxg,:)'; gpb(idxg,:)'];   %Green fun for pl to sel mag
%R:
   gmpr = [mpl(idxg,idxflr)'; gpb(idxg,idxbpr)'];   %Green fun for pl to sel mag
   gmcr = [mlc(idxflr,idxc); gbc(idxbpr,idxc)];   %Green fun for PF's to sel mag
   gmvr = [mlv(idxflr,idxv); gbv(idxbpr,idxv)];   %Green fun for VV's to sel mag
%Z:
   gmpz = [mpl(idxg,idxflz)'; gpb(idxg,idxbpz)'];   %Green fun for pl to sel mag
   gmcz = [mlc(idxflz,idxc); gbc(idxbpz,idxc)];   %Green fun for PF's to sel mag
   gmvz = [mlv(idxflz,idxv); gbv(idxbpz,idxv)];   %Green fun for VV's to sel mag
%Ip:
   gmpip = [mpl(idxg,idxflip)'; gpb(idxg,idxbpip)'];   %sel mag
   gmcip = [mlc(idxflip,idxc); gbc(idxbpip,idxc)];   %sel mag
   gmvip = [mlv(idxflip,idxv); gbv(idxbpip,idxv)];   %sel mag
%Kappa:
   gmpk = [mpl(idxg,idxflk)'; gpb(idxg,idxbpk)'];   %Green fun for pl to sel mag
   gmck = [mlc(idxflk,idxc); gbc(idxbpk,idxc)];   %Green fun for PF's to sel mag
   gmvk = [mlv(idxflk,idxv); gbv(idxbpk,idxv)];   %Green fun for VV's to sel mag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build basic predictors, make test data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Ip normalization quantity:
   if imks==1
     ipnorm=1e6;
    else
     ipnorm=1;
   end
   if imks==1, curstr='[A]'; else curstr='[MA]'; end  %label string for Ip

%Basic linear centroid predictors (must divide by Ip for correct prediction):
%Ppred2: Use plasma region + Coils:
 %Z:
   if ifitpfz==1
     if ipinvz
      tmp = pinv([gmpz gmcz]); 
     else
      tmp = minvert([gmpz gmcz],nsvdz2); 
     end
   else
     if ipinvz
      tmp = pinv(gmpz);   
     else
      tmp = minvert(gmpz,nsvdz2);   
     end
   end	
   tmp=tmp(1:ngx,:);  %extract cur pred
   Pzpred2 = zgg(idxg)'*tmp;  %m-MA/mag (if imks=0)
   dum = rank(tmp); disp(['Pzpred kept ',int2str(dum),' singular values'])

%R:
   if ifitpfr==1
     if ipinvr
      tmp = pinv([gmpr gmcr]); 
     else
      tmp = minvert([gmpr gmcr],nsvdr2); 
     end
   else
     if ipinvr
      tmp = pinv([gmpr]); 
     else
      tmp = minvert([gmpr],nsvdr2); 
     end
   end	
   tmp=tmp(1:ngx,:);  %extract cur pred
   Prpred2 = (rgg(idxg)'-Rrel)*tmp;  %m-MA/mag (if imks=0)
   dum = rank(tmp); disp(['Prpred kept ',int2str(dum),' singular values'])

%Ip:
   if ifitpfip==1
     if ipinvip
       tmp = pinv([gpb(idxg,idxbpip)' gbc(idxbpip,idxc)]);
     else
       tmp = minvert([gpb(idxg,idxbpip)' gbc(idxbpip,idxc)],nsvdi2);
     end
   else
     if ipinvip
       tmp = pinv([gpb(idxg,idxbpip)']);
     else
       tmp = minvert([gpb(idxg,idxbpip)'],nsvdi2);
     end
   end	
   tmp=tmp(1:ngx,:);  %extract cur pred
   Pipred2 = sum(tmp); %MA/T (only probes for Ip)
   dum = rank(tmp); 
       disp(['Pipred kept ',int2str(dum),' singular values'])

%Make magnetics test data for entire grid:
 %Test1 = current of 1 (MA or A) in one grid element at a time:
   curdat1 = eye(ngx);   %curr in each sel grid el
   zdat1 = (zgg(idxg)'*curdat1)';   %cur centroids for curdat1
   rdat1 = ((rgg(idxg)'-Rrel)*curdat1)';
   ipdat1 = sum(curdat1)';       %Ip for each el on grid
   ztst1 = (Pzpred2*(gmpz*curdat1))';
   rtst1 = (Prpred2*(gmpr*curdat1))';
   iptst1 = (Pipred2*(gpb(idxg,idxbpip)'*curdat1))';   %estimated Ip for each el on grid
   zerr1 = abs(zdat1-ztst1)/mean(abs(zdat1));  %measure of Z-error
   rerr1 = abs(rdat1-rtst1)/mean(abs(rdat1));  %measure of R-error
   iperr1 = abs(ipdat1-iptst1)/mean(abs(ipdat1));  %measure of Ip-error ongrid
 %Test2 = principal components of mag grid currents
   [U,sig,V] = svd(gmp');    %SVD: cols of U are princ comps of grid currs
   curdat2 = abs(U(:,1:nprtst));
   curdat2 = curdat2./(ones(ngx,1)*(sum(curdat2)));   
   zdat2 = (zgg(idxg)'*curdat2)';   %cur centroids for curdat2
   rdat2 = ((rgg(idxg)'-Rrel)*curdat2)';   
   ipdat2 = sum(curdat2)';       %Ip for each princ component distrib
   ztst2 = (Pzpred2*(gmpz*curdat2))';
   rtst2 = (Prpred2*(gmpr*curdat2))';
   iptst2 = (Pipred2*(gpb(idxg,idxbpip)'*curdat2))';   %estimated Ip for each princ component
   zerr2 = abs(zdat2-ztst2)/mean(abs(zdat2));  %measure of Z-error
   rerr2 = abs(rdat2-rtst2)/mean(abs(rdat2));  %measure of R-error
   iperr2 = abs(ipdat2-iptst2)/mean(abs(ipdat2));  %measure of Ip-error ongrid
   
%Plot basic comparisons: Z,R,Ip
figure(2),clf,hold off
zdat=zdat1; rdat=rdat1; ipdat=ipdat1;
ztst=ztst1; rtst=rtst1; iptst=iptst1;
zerr=zerr1; rerr=rerr1; iperr=iperr1;
nrows=3; ncols=2; 
subplot(nrows,ncols,1)   %Zest
  plot(zdat,ztst,'+')
  hold on
  axis square
  axis0=axis;
  mi=min(axis0(1),axis0(3));
  ma=max(axis0(2),axis0(4));
  plot([mi ma],[mi ma],'c')
  grid on
  ylabel('Ztst [m]')
  xlabel('Zdat [m]')
title('Sub-grid Elements Predictor Test','FontSize',15)
subplot(nrows,ncols,3)   %Rest
  plot(rdat,rtst,'+')
  hold on
  axis square
  axis0=axis;
  mi=min(axis0(1),axis0(3));
  ma=max(axis0(2),axis0(4));
  plot([mi ma],[mi ma],'c')
  grid on
  ylabel('Rtst [m]')
  xlabel('Rdat [m]')
subplot(nrows,ncols,5)   %Ipest
  plot(ipdat,iptst,'+')
  hold on
  axis square
  axis0=axis;
  mi=min(axis0(1),axis0(3));
  ma=max(axis0(2),axis0(4));
  plot([mi ma],[mi ma],'c')
  grid on
  ylabel(['Iptst ',curstr])
  xlabel(['Ipdat ',curstr])
subplot(nrows,ncols,2)   %Zerr
  plot(zdat,zerr*100,'+')
  grid on
  ylabel('Zerr/Zref [%]')
  xlabel('Zdat [m]')
subplot(nrows,ncols,4)   %Rerr
  plot(rdat,rerr*100,'+')
  grid on
  ylabel('Rerr/Rref [%]')
  xlabel('Rdat [m]')
subplot(nrows,ncols,6)   %Iperr
  plot(ipdat,iperr*100,'+')
  grid on
  ylabel('Iperr/Ipref [%]')
  xlabel(['Ipdat ',curstr])


%Pred1's with curdat2 (SVD principal components):
figure(3),clf,hold off
zdat=zdat2; rdat=rdat2; ipdat=ipdat2;
ztst=ztst2; rtst=rtst2; iptst=iptst2;
zerr=zerr2; rerr=rerr2; iperr=iperr2;
nrows=3; ncols=2; 
subplot(nrows,ncols,1)
  plot(zdat,ztst,'+')
  hold on
  axis square
  axis0=axis;
  mi=min(axis0(1),axis0(3));
  ma=max(axis0(2),axis0(4));
  plot([mi ma],[mi ma],'c')
  grid on
  ylabel('Ztst [m]')
  xlabel('Zdat [m]')
title('Principal Component Predictor Test','FontSize',15)
subplot(nrows,ncols,3)
  plot(rdat,rtst,'+')
  hold on
  axis square
  axis0=axis;
  mi=min(axis0(1),axis0(3));
  ma=max(axis0(2),axis0(4));
  plot([mi ma],[mi ma],'c')
  grid on
  ylabel('Rtst [m]')
  xlabel('Rdat [m]')
subplot(nrows,ncols,5)
  plot(ipdat,iptst,'+')
  hold on
  axis square
  axis0=axis;
  mi=min(axis0(1),axis0(3));
  ma=max(axis0(2),axis0(4));
  plot([mi ma],[mi ma],'c')
  grid on
  ylabel(['Iptst ',curstr])
  xlabel(['Ipdat ',curstr])
subplot(nrows,ncols,2)
  plot(zdat,zerr*100,'+')
  grid on
  ylabel('Zerr/Zref [%]')
  xlabel('Zdat [m]')
subplot(nrows,ncols,4)
  plot(rdat,rerr*100,'+')
  grid on
  ylabel('Rerr/Rref [%]')
  xlabel('Rdat [m]')
subplot(nrows,ncols,6)
  plot(ipdat,iperr*100,'+')
  grid on
  ylabel('Iperr/Ipref [%]')
  xlabel(['Ipdat ',curstr])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kappa predictor calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Beginning kappa predictor calculation..')

% Kappa pred 1: Statistical ensemble (but deterministic, not random):

%-------------------
if 0   %DISABLED: these are defined in make_*_predictors now
  %Inputs:
	R1 = 1.36;   %inboard limiter for plasmas in kappa data set
	R2 = 2.35;   %outboard limiter for plasmas in kappa data set
	Z1 = -0.9;   %lower limiter for plasmas in kappa data set
	Z2 = 0.9;    %upper limiter for plasmas in kappa data set
	amax = 0.5;  %max allowed minor radius in kappa data set
	nrs = 5;
	Rmin = 1.75;  %min maj rad of plasmas in kappa data set (1.6)
	Rmax = 1.85;  %max maj rad of plasmas in kappa data set (2.0)
	nzs = 5;
	Zmin = -0.05;
	Zmax = 0.05;
	nkaps = 20;  
	kapmin = 1.0; 
	kapmax = 2.0;
	nlis = 1;
        limin = 1.0;
        limax = 1.0;
	nsvdk1 = 30;
end  %end if 0 disabling inputs
%-------------------

  %Magnetics mapping:
    gmpkap = [mpl(:,idxflk)'; gpb(:,idxbpk)'];   %Green fun for pl to sel mag

  %Generate data with makelpcur:
    nkapdata = nrs*nzs*nkaps*nlis;
    nmag = length(idxflk) + length(idxbpk);
    Rs = linspace(Rmin,Rmax,nrs)';
    Zs = linspace(Zmin,Zmax,nzs)';
    kaps = linspace(kapmin,kapmax,nkaps)';
    lis = linspace(limin,limax,nlis)';
    kappadat = zeros(nkapdata,1);    %vector of all kappa data in set
    gindat = zeros(nkapdata,1);    %vector of all kappa data in set
    magdatkap = zeros(nmag,nkapdata); %matrix of all magdata
    ii = 0;   %initialize counter
    for iir=1:nrs
      for iiz=1:nzs
	for iik=1:nkaps
    	 for iil=1:nlis
 	  ii = ii+1;
	  R=Rs(iir); Z=Zs(iiz); K=kaps(iik); li=lis(iil);
  	  a = min([abs(R-R1) abs(R-R2) amax]);  %limit in/out
	  if (Z+K*a)>Z2, a=(Z2-Z)/K; end       %limit up
	  if (Z-K*a)<Z1, a=(Z-Z1)/K; end       %limit down
  	  gin = R-a-R1;   %inner gap
	  cphi = makelpcur(zgg,rgg,1,Z,R,a,K,li); %make current
	  gindat(ii) = gin;
	  kappadat(ii) = kaps(iik);
	  magdatkap(:,ii) = gmpkap*cphi(:);
	 end
	end
      end
    end
    if ipinvk
       Pkpred1 = kappadat'*pinv(magdatkap);
    else
       Pkpred1 = kappadat'*minvert(magdatkap,nsvdk1);
    end

%Test kappa predictor1:
  kaptst = Pkpred1*magdatkap;
  figure(4),clf
    plot(kappadat,kaptst,'x')
    hold on
    axis0=axis;
    mi=min(axis0(1),axis0(3));
    ma=max(axis0(2),axis0(4));
    plot([mi ma],[mi ma],'r')
    axis([mi ma mi ma]), axis('square')
    grid on
    title('Kappa Predictor 1 Comparison','FontSize',15)
    ylabel('Kappa Estimate','FontSize',14)
    xlabel('Kappa Data','FontSize',14)


% Current-distrib based kappa:
    %NOT DONE

disp('Ending kappa predictor calculation..')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test predictors against EFIT jphi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get EFIT data:
   filename = efit_gfile_test;
   gfile_data = read_gfile_tok(filename,tokamak);
   struct_to_ws(gfile_data);   %need for plot_geo later...
   dz = zg(2)-zg(1);
   dr = rg(2)-rg(1);
   cphi = jphi*dz*dr*ipnorm; %actual current in grid elements [MA] or A
   ip0 = sum(cphi(:));   %cphi turned into A if mks...
   if ~exist('cc') 
	disp('No cc exists: using cctest') 
        cc = cctest;  %musst have defined in make_*_predictors
   end
   cc = cc*ipnorm;    %turn cc into A if mks...
	disp('cctest needed but dim not matching mcc...');
   if length(cc)~=size(mcc,1); 
      disp('cc dim not matching mcc...')
      if length(cctest)~=size(mcc,1)   %cctest must be size(mcc,1)
	disp('cctest needed but dim not matching mcc...');
	return
      end
      ccx = cc;
      cc = cctest;     
      cc(1:length(fcid)) = ccx(fcid).*fcturn;
   end
   if imks==1
      cc = cc./fcnturn;    %turn cc into terminal current...
   end

%Generate magnetics and centroid position from EFIT current:
   %magdat3 = gmp*cphi(idxg);       %generated mag data from EFIT cphi 
   %bpdat3ip = gpb(idxg,idxbpip)'*cphi(idxg);  %generated Bp data from EFIT cphi
   magdat4r = gmpr*cphi(idxg)+0*gmcr*cc(idxc); %generated data f/EFIT cphi+cc
   magdat4z = gmpz*cphi(idxg)+0*gmcz*cc(idxc); %generated data f/EFIT cphi+cc
   bpdat4ip = gpb(idxg,idxbpip)'*cphi(idxg)+0*gbc(idxbpip,idxc)*cc(idxc);  %cphi+cc
   magdat4k = gmpk*cphi(idxg)+0*gmck*cc(idxc); %generated data f/EFIT cphi+cc
   zcefit = zgg(:)'*cphi(:)/ip0;   %actual Zcur from EFIT
   rcefit = (rgg(:)'-Rrel)*cphi(:)/ip0;   %actual Rcur from EFIT RELATIVE TO Rrel!!!
   idx = find(cphi(:)~=0);         %plasma current region
   hplasma = max(max(zgg(idx))) - min(min(zgg(idx)));
   wplasma = max(max(rgg(idx))) - min(min(rgg(idx)));
   kappaefit = hplasma/wplasma;

%Apply predictor to EFIT  magnetics data:
   ztst2_4 = Pzpred2*magdat4z/ip0; %Test of EFIT with Pred2 vs total EFIT mags
   rtst2_4 = Prpred2*magdat4r/ip0;
   ktst1_4 = Pkpred1*magdat4k/ip0;

%Apply Ip predictor to EFIT Bprobe data:
   iptst2_4 = Pipred2*bpdat4ip;

%Report EFIT test results:
   disp('Predictions using cphi+cc pred + mags from EFIT cphi+cc (magdat4):')
   disp(['Zcur(EFIT,pred2_4) = ',num2str(zcefit),' ',num2str(ztst2_4)])
   disp(['Rcur(EFIT,pred2_4) = ',num2str(rcefit),' ',num2str(rtst2_4)])
   disp(['Ip(EFIT,pred2_4) = ',num2str(ip0),' ',num2str(iptst2_4)])
   disp(['Kappa(EFIT,pred2_4) = ',num2str(kappaefit),' ',num2str(ktst1_4)])

%Plot plasma boundary over sub-grid to check:
 figure(1),hold on
   iplteq=1; ipltflx=1;
   ipltlim=1; 
   plot_tok_geo(tok_data_struct)
   title('Sub-grid Locations','FontSize',15)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tests with equilibrium coil currents cc and VV currents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Entering tests with cc and vv currents...')
 %Test5 = PF and VV currents only, no plasma (should give ~zero result for all)
  %%cols = Bp data for PF's, VV's:
   bpdat5ip = [gbc(idxbpip,:)*diag(cc) gbv(idxbpip,:)]; 
   ztst5 = (Pzpred2*[gmcz*diag(cc) gmvz])'/ip0;  %must normalize by source of magdat5...
   rtst5 = (Prpred2*[gmcr*diag(cc) gmvr])'/ip0;
   iptst5 = (Pipred2*bpdat5ip)';   %estimated Ip for each princ component

figure(5),clf,hold off
ztst=ztst5; rtst=rtst5; iptst=iptst5;
nrows=3; ncols=1;
subplot(nrows,ncols,1)
  plot(ztst,'+')
  grid on
  ylabel('Ztst [m]')
  xlabel('Index (PF+VV)')
title('Test with PF+VV Currents Alone','FontSize',15)
subplot(nrows,ncols,2)
  plot(rtst,'+')
  grid on
  ylabel('Rtst [m]')
  xlabel('Index (PF+VV)')
subplot(nrows,ncols,3)
  plot(iptst,'+')
  grid on
  ylabel(['Iptst ',curstr])
  xlabel('Index (PF+VV)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final predictor assembly and reporting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Make PF compensation matrix(reorder PF,IC list for PCS order): 
  reorder_pfic =[1:18];	% configuration dependent - don't do this here.
  compcmat = [mlc(:,reorder_pfic); gbc(:,reorder_pfic)];

%Make predictors fill vector for full set of magnetic diagnostics:
  Prpred=zeros(nfl+nbp,1);
     Prpred(idxflr)=Prpred2(1:length(idxflr));
     Prpred(nfl+idxbpr)=Prpred2(length(idxflr)+(1:length(idxbpr)));
     Rpred = [-icomppfr*(Prpred'*compcmat)'; Prpred];   %pre-pend coil currents

  Pzpred=zeros(nfl+nbp,1);
     Pzpred(idxflz)=Pzpred2(1:length(idxflz));
     Pzpred(nfl+idxbpz)=Pzpred2(length(idxflz)+(1:length(idxbpz)));
     Zpred = [-icomppfz*(Pzpred'*compcmat)'; Pzpred];   %pre-pend coil currents

  Pipred=zeros(nfl+nbp,1);
     Pipred(nfl+idxbpip)=Pipred2(1:length(idxbpip));
     Ippred = [-icomppfip*(Pipred'*compcmat)'; Pipred];   %pre-pend coil currents

  Pkpred=zeros(nfl+nbp,1);
     Pkpred(idxflk)=Pkpred1(1:length(idxflk));
     Pkpred(nfl+idxbpk)=Pkpred1(length(idxflk)+(1:length(idxbpk)));
     Kpred = [-icomppfk*(Pkpred'*compcmat)';  Pkpred];   %pre-pend coil currents

descriptions = struct(...
'Prpred','Z-predictor (Pzpred*[mag]/Ip0 = predicted Z)', ...
'Pzpred','R-predictor (Prpred*[mag]/Ip0 = predicted R)', ...
'Pipred','Ip-predictor (Pipred*[Bprobes] = predicted Ip)', ...
'Pkpred','kappa-predictor (Pkpred*[mag]/Ip0 = predicted kappa)', ...
'Rpred','same as Prpred but with coil currents at beginning?', ...
'Zpred','same as Pzpred but with coil currents at beginning?', ...
'Ippred','same as Pipred but with coil currents at beginning?', ...
'Kpred','same as Pkpred but with coil currents at beginning?', ...
'inputs','user input parameters to construction of predictors');

predictors = struct(...
'Prpred',Prpred, ...
'Pzpred',Pzpred, ...
'Pipred',Pipred, ...
'Pkpred',Pkpred, ...
'Rpred',Rpred, ...
'Zpred',Zpred, ...
'Ippred',Ippred, ...
'Kpred',Kpred, ...
'inputs',predictor_inputs, ...
'descriptions',descriptions);

%Report predictors:
disp('Predictors are: Rpred, Zpred, Ippred, Kpred')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   disp('All done.')
               
