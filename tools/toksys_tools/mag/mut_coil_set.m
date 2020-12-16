  function mut= mut_coil_set(rc,zc,dr,dz)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAXX:
%         mut= mut_coil_set(rc,zc,dr,dz)
%
%  PURPOSE:  Generate the mutual inductance between all coils in a set
%            Coils are coaxial and axisymmetric, have rectangular 
%            cross-section and carry uniform current density
%
%  INPUT:
%	rc  = Radial center of coil [m]
%	zc  = Axial center of coil [m]
%	dr  = Radial width of coil [m]
%	dz  = Axial height of coil [m]
%
%  OUTPUT:
%	mut = [n,n] matrix of mutual inductance of coils [Henries]
%
%  RESTRICTIONS:
%       See rectl for Self Inductance Restriction: err<1% for delz/(2*delr) < 1
%       Coils with aspect ratio approaching 1 are problematic

%  METHOD:  
%
%  WRITTEN BY: Jim Leuer 8/30/99  
%  uses fine and integration used in mutrect.for
%  see \field\coil_mut_ind.m for use:  read coil set and plot
%  See JAL D3D Book 19
%  jal16May2007 fix of transpose problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  3rd ord.  LEGENDRE QUADRATURE
   wgt= [ 0.1739275, 0.3260725, 0.3260725, 0.1739275]; 
   del= [-0.4305680,-0.1699905, 0.1699905, 0.4305680];

   n= length(rc);
   amu= 4*pi*1e-7;

% do diagonal self inductance

    mut= rectl(rc,dz,dr)*1e-6;
    mut= diag(mut);

% do off diagonal elements
   for ii= 1:n-1 %sending coil
     ri= (rc(ii) - 0.5*dr(ii))*ones(16,1);
     ro= (rc(ii) + 0.5*dr(ii))*ones(16,1);
     zl= (zc(ii) - 0.5*dz(ii))*ones(16,1);
     zu= (zc(ii) + 0.5*dz(ii))*ones(16,1);

     for jj=ii+1:n % receiving coil
       rcc= rc(jj);
       zcc= zc(jj);
       drr= dr(jj);
       dzz= dz(jj);
       id= 0;
       for iii= 1:4
         zpp= zcc + del(iii)*dzz;
         for jjj= 1:4
           id= id+1;
           zp(id,:)= zpp;
           rp(id,:)= rcc + del(jjj)*drr;
           cu(id,:)= wgt(iii)*wgt(jjj); 
         end % jjj
       end % iii

       [hr,hz,fl] =fine(ri,ro,zl,zu,cu,rp,zp);
       mut(ii,jj)= sum(fl)*amu;
       mut(jj,ii)= mut(ii,jj);
     end % for j

   end % for i

  return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------------------------------------------------
% Test 1  two coils on top of each other

  rc= [1;1]; zc=[-.25;.25]; dr=[.5;.5]; dz=[.5;.5]

  m= 1e+6*mut_coil_set(rc,zc,dr,dz)

%  m=
%
%   2.03039935700305   1.12770797982673
%   1.12770797982673   2.03038251408063

% Check below using fortran:

%$ run cforce2
%INPUT Number of integration grid: Nr,Nz4,4
% INPUT From Coil: Rcfrom,Zcfrom,Drfrom,Dzfrom,ampfrom: 1,.25,.5,.5,1
% INPUT To Coil: Rcto(-=>from),Zcto,Drto,Dzto,Ampto: 1,-.25,.5,.5,1
% Mutave,Mutcent  =   1.1269728E-06  1.1189655E-06 -0.7105145    
%
%$ run drectl
%
%INPUT: a(radius), c(radial Build), b(height) 1,.5,.5
%
%
%dRECTL REPORTS (Inductances):
% a(rad), b(hgt), c(rad):  1.00000     0.500000     0.500000    
% subroutine:                 L           dL/da     d^2L/da^2   
%    RECTL                0.203038E-05
%    RECTL1               0.203038E-05 0.320444E-05 0.131312E-05
%    RECTL2(itype=1)      0.203038E-05 0.320435E-05 0.131309E-05
%    RECTL2(itype=2)      0.203038E-05 0.261182E-05 0.342198E-06


