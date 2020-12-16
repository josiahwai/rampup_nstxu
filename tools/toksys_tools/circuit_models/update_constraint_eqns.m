function [LHS,RHS,ier] = update_constraint_eqns(LHS,RHS,nvar,addLHS,addRHS)
 %
%  SYNTAX: update_constraint_eqns 
%
%  PURPOSE:  Update the constraint equations in dynamic_equil.m.
%
%  INPUT:
%    LHS, RHS = matrices defining previous constraint equation (LHS*x=RHS)
%    nvar = number of variables in optimization (size of x)
%    addLHS, addRHS = additional rows to add to constraint equation
%		(# columns in addLHS must be < = # columns in LHS)
%
%  OUTPUT:
%    LHS, RHS = matrices for updated constraint equations
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	7/1/09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)update_constraint_eqns.m	1.2 05/05/10

ier = 0;
[l1,l2]=size(LHS);
[r1,r2]=size(RHS);
[al1,al2]=size(addLHS);
[ar1,ar2]=size(addRHS);
tempLHS = addLHS; tempRHS = addRHS;

if(~isempty(LHS) & l2~=nvar)
   wait('ERROR update_constraint_eqns - size(LHS,2) must be = nvar')
   ier = 1;
   return
end

if(al1~=ar1)
   wait('ERROR update_constraint_eqns - size(addLHS,1) must be = size(addRHS,1)')
   ier = 2;
   return
end

if(l1~=r1)
   wait('ERROR update_constraint_eqns - size(LHS,1) must be = size(RHS,1)')
   ier = 3;
   return
end

if(al2>nvar)
   wait('ERROR update_constraint_eqns - size(addLHS,2) must be <= nvar')
   ier = 4;
   return
end

if(ar2>1)
   wait('ERROR update_constraint_eqns - RHS must be column vector')
   ier = 5;
   return
end

if(al2<nvar)
   idiff = nvar-al2;
   tempLHS = [tempLHS zeros(al1,idiff)];
end

LHS = [LHS; tempLHS];
RHS = [RHS; tempRHS];
