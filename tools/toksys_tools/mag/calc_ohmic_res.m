  function rcc= calc_ohmic_res(mcc,icc)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PURPOSE: Calculates Best Coil resistances to produce constant flux decay
%
%  USAGE:  rcc= calc_ohmic_res(mcc,icc);
%
%  INPUTS: [default]
%           
%     mcc=      coil mutual inductance matrix
%     icc=	Optimum initial magnetization current vector
%
%  OUTPUTS:
%
%     rcc =   optimum resistance matrix
%
%  METHOD:
%     [I]^-1 is diaganal matrix made up of 1/I_i currents
%     R/gamma =  [I]^-1 [M] I
%
% LIMITATIONS:
%    Assumes the initial current distrubution is Ohmic (i.e. produces constant
%            flux plateau in plasma region.
%            icc is typically generated using ohmic_dist
%
%    see also ohmic_dist.m

%  WRITTEN BY:  Jim Leuer  ON   1/20/2006
% could be easily modified to produce arbitrary flux distribution
%
  rcc= diag(1./icc)*mcc*icc;

  return

% =================================
% testing

  clf

  if exist('tok_data_struct')~=1
    thordisk=getenv('THORDISK');
    str= ['load ' thordisk '/walker/tokamaks/EAST/make_objects/east_obj_struct'];
    eval(str)
%    load /home/walker/tokamaks/EAST/make_objects/east_obj_struct
  end
  
  mcc= tok_data_struct.mcc;
  nc= 1:length(mcc)-2; % remove control coils
  nc= [1 2 3 5 6 7 8 9 10 12 13 14]
  mcc= mcc(nc,nc);
  resc= tok_data_struct.resc;
  resc= resc(nc);
  ro= mean(tok_data_struct.rg);
  zo= 0;
  a= ro/3;
  k= 2;
  d= 0;
  

  theta= 0:pi/50:2*pi;
  theta= theta(1:end-1);
  
  [r,z]= dee(ro,zo,a,k,d,theta);
  r=r';
  z=z';
  
  plot(r,z)
  axis image

  fcdata= tok_data_struct.fcdata;
  fcdata= fcdata(:,nc);
  
  rcc= fcdata(2,:)';
  zcc= fcdata(1,:)';
  drcc= fcdata(4,:)';
  dzcc= fcdata(3,:)';
  hold on
  plot_box(rcc,zcc,drcc,dzcc) 
  axis auto
  axis equal

% [i_ohmic,fl,br,bz,rbs,zbs,err] = ...
%         ohmic_dist(r,z,rcc,zcc,drcc,dzcc,mk_bs_uniform,minimize_i)
  

 [i_ohmic,fl,br,bz,rbs,zbs,err] = ohmic_dist(r,z,rcc,zcc,drcc,dzcc);

  bb= sqrt(br.^2+bz.^2);
  
  icc= i_ohmic;
  rcc= calc_ohmic_res(mcc,icc); 
  
  amat= -mcc^-1*diag(rcc);
  bmat= zeros(length(amat),1);
  cmat= eye(size(amat));
  dmat= zeros(size(cmat,1),1);
  
 
  tvec = (0:.01:1)';
  u= zeros(length(tvec),1);

  sys = ss(amat,bmat,cmat,dmat);

  [y,t] = lsim(sys,u,tvec,icc);
  
  dy= y ./ (ones(length(tvec),1)*y(1,:));

 % validation complete: below shows all lines overlay which is ideal solution
  plot(t,dy); % this plots all current vectors normalized to initial current
  


