
 function [vvdata,idirnu,iacnu,fracnu]=lim2vv(limdata,dperp,dw,idir,frac,iac);

 %
%  SYNTAX: [vvdata,dirnu,iacnu,fracnu]=lim2vv(limdata,dperp,dw,idir,frac,iac);
%
%  PURPOSE:  
%	Calculate vessel elements from limiter specification given by
%	limdata = set of vertices in limiter, between which will be located
%	vessel elements. 
%
%  INPUTS:
%	limdata = set of vertices in limiter, between which will be located
%		vessel elements. limdata has format [[r1 z1];[r2 z2]...]
%	dperp = amount vessel elements are displaced from limiter
%	dw = width of vessel elements
%	idir = (optional) direction indices corr to each vessel element 
%		(nvv = nlim - 1). Default corresponds to +1; reverses 
%		direction of wall displacement if set element of idir to -1.
%	frac = (optional) fraction by which to scale each VVelement
%	iac = (optional) indices telling if AC or AC2 type element.
%		Default = all AC (index 1) elements. Note that automatically
%		sets vertical or horizontal elements to type 0...
%
%  OUTPUTS:
%	vvdata = vvdata array with standard format
%	idirnu = new idir vector (actually used)
%	iacnu = new iac vector (which was used in specifying actual vvdata)
%	fracnu = new frac vector (actually used)
%	plots of corresponding geometry
%
%  RESTRICTIONS:
%
%  METHOD:  

%  WRITTEN BY:  Dave Humphreys 	ON	1/11/00
%
%  MODIFICATION HISTORY:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)lim2vv.m	1.2 02/23/09

% Prelims:
   mu0 = 0.4*pi;

% Derived Values:
   xlim = limdata(2,:);
   ylim = limdata(1,:);
   nlim = length(xlim);
   nvv = nlim - 1;
   if nargin<4, idir=ones(nvv,1); end
   if nargin<5, frac=ones(nvv,1); end   
   if nargin<6, iac=ones(nvv,1); end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iterate over VV elements (pairs of limiter vertices)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 vvdata = zeros(6,nvv);
 for ii=1:nlim-1
   m=1; minv=1;  %initializing to non-zero value for correct logic below
   x1=xlim(ii); x2=xlim(ii+1); x0=0.5*(x1+x2);
   y1=ylim(ii); y2=ylim(ii+1); y0=0.5*(y1+y2);
   len = sqrt((y2-y1)^2 + (x2-x1)^2); 
   if x2-x1==0
     minv=0; ang=pi/2;    %ang = cw from horizontal?
    else 
     m=(y2-y1)/(x2-x1); ang=atan2(y2-y1,x2-x1);
   end
   if minv==0    %then vertical element (iac=2)
     x0nu=x0+dperp*idir(ii); y0nu=y0;   %if vertical element
    else
     x0nu=x0+dperp*sin(ang)*idir(ii); y0nu=y0-dperp*cos(ang)*idir(ii);
   end
   vvdata(1,ii)=y0nu; vvdata(2,ii)=x0nu;
   if nargin<6|iac(ii)==0  %if no iac input OR input=0, determine iac f/ geom
     if (ang>=pi/4&ang<=3*pi/4)|(ang<=-pi/4&ang>=-3*pi/4)
       iac(ii)=2;   %then it's an AC2 (flat bottom)
      else
      iac(ii)=1;   %then it's an AC (flat side)
     end
   end   %end if nargin:  otherwise, leave as input
   if m==0, iac(ii)=1; vvdata(6,ii)=0; vvdata(5,ii)=0; end
   if minv==0, iac(ii)=2; vvdata(6,ii)=0; vvdata(5,ii)=0; end
   if iac(ii)==1
     vvdata(3,ii)=abs(dw/cos(ang));
     vvdata(4,ii)=abs(frac(ii)*len*cos(ang));
     vvdata(5,ii)=ang*180/pi; 
   end
   if iac(ii)==2
     vvdata(3,ii)=abs(frac(ii)*len*sin(ang));
     vvdata(4,ii)=abs(dw/sin(ang));
     vvdata(6,ii)=ang*180/pi; 
   end
 end   %end for ii
 idirnu = idir;
 iacnu = iac;
 fracnu = frac;

%VVdata object for inserting into vv.data file (R,z order):
 vv_data=[vvdata(2,:);vvdata(1,:);vvdata(4,:);vvdata(3,:);vvdata(5:6,:)];
 disp('vv_data = ')
 disp(vv_data')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot resulting geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure(1),clf,hold off
  plot_toksys_geo
  hold on
  plot(limdata(2,:),limdata(1,:),'yx')
  plot(vvdata(2,:),vvdata(1,:),'gv')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   disp('All done.')

