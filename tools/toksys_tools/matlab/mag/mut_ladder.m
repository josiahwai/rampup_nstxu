  function mut= mut_ladder(a,b,r,n)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAXX:
%         mut= mut_ladder(a,b,r,n)
%
%  PURPOSE:  Generate the mutual inductance between cells in a ladder shaped
%            component made of n cells in the a direction.
%            Each cell has dimensions in the x by y direction of a by b
%            Cells are replicated in the x direction resulting in a ladder
%            structure of length n*a in the x direction and width b in y.
%            Each cell is adjacent to the next cell in the x direction and
%            the rung is common to the left and right cells
%
%  INPUT:
%	a   = side dimension of one cell in the length ladder direction, [m]
%	b   = sied dimension of one cell in the width direction, [m]
%       r   = radius or wire cross (vector) [m]
%       n   = number of rectangular cells
%
%  OUTPUT:
%	mut = [n,n] matrix of mutual inductance of cells in ladder [Henries]
%             (Order is 1st cell on left to last cell on right)
%
%  RESTRICTIONS:
%        At present all cells are same size; but could be modified to do
%        variable cell size in x direction by changing vector x internally

%  METHOD:  
%
%  WRITTEN BY: Jim Leuer 8/25/99  
%  Utilizes summation principle to determine mutual from self-inductance
%  of common components. Reference: Grover
%  this is real tricky - see Book 19 8/29/99 See Final Ladder
%  Note for all positive currents (counter-clockwise) in the ladder cells
%  the mutual inductance is negative. This is because the flux in cell B
%  from positive current in cell A is negative relative to the flux in cell B
%  from its self current Ib. That is: Flux_ba = Mba*Ia = (negative). So for
%  positive Ia, Mba must be negative. 
%  in addition the number of turns does not come into play as is does in
%  solenoidal system since the coils are adjacent rather than axial.
%  Configuration for 2 cells (A&B) seperated by cell C with arrows represnting
%  currents (counter clockwise currents are positive)
%
%   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%  v           ^ v           ^ v           ^
%  v           ^ v           ^ v           ^
%  v           ^ v           ^ v           ^
%  v    A      ^ v    C      ^ v    B      ^
%  v           ^ v           ^ v           ^
%  v           ^ v           ^ v           ^
%  v>>>>>>>>>>>> >>>>>>>>>>>>> >>>>>>>>>>>>>
%
% This configureation can be used to show that
% 
% 2*Mab= Labc + Lc - Lac - Lcb
%  See JAL D3D Book 19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   dx= a;
   dy= b;
   xp= dx*(1:n); % modify xp and xm to make variable rung ladder
   xm= xp - dx;
   mut= zeros(n,n);

% do diagonal self inductance

  for i=1:n
    mut(i,i)= induc_rect(dy,xp(i)-xm(i),r);
  end

% do off diagonal elements
   for i=1:n-1
     dxa= xp(i)-xm(i);
     for j=i+1:n
       dxb=    xp(j)-xm(j);
       dxc=    xm(j)-xp(i);
       dxabc=  xp(j)-xm(i);
       dxac=   dxa+dxc;
       dxbc=   dxb+dxc;
       labc=   induc_rect(dy,dxabc,r);
       lac=    induc_rect(dy,dxac,r);
       lbc=    induc_rect(dy,dxbc,r);
       lc=     0;
       if j-i-1 > eps
          lc = induc_rect(dy,dxc,r);
       end
       mut(i,j)= 0.5*( labc + lc - lac  - lbc );
       mut(j,i)= mut(i,j);
     end % for j
   end % for i

  return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------------------------------------------------
% Test 1 3 coil (only tests implementation not theory)

  a= 1; b=1; r=0.1; n=3;

  mut= mut_ladder(a,b,r,n);

  la=   induc_rect(a,b,r);
  lab=  induc_rect(2*a,b,r);
  labc= induc_rect(3*a,b,r);

  mut

  m11= la;
  m21= 0.5*(lab-2*la);
  m31= 0.5*(labc+la-2*lab);

  test= [m11,m21,m31;m21,m11,m21;m31,m21,m11]
          
