function [Amod,nodelist,branchlist,n] = echelon(A,nodelist,branchlist,n)
 %
%  SYNTAX: [Amod,nodelist,branchlist] = echelon(A,nodelist,branchlist)
%
%  PURPOSE:  Echelon algorithm.  Calculation to define branches in tree and
%	co-tree and the matrix mapping currents in co-tree branches to currents
%	in branches of tree.
%
%  INPUT:
%   A          = node incidence matrix
%   nodelist   = node numbers corresponding to each row of A
%   branchlist = branch indices into the netlist corresponding to columns of A
%
%  OUTPUT:
%   Amod       = modified incidence matrix, of the form [I Dc]
%   nodelist   = node numbers corresponding to each row of Amod
%   branchlist = netlist branch indices corresponding to columns of Amod
%
%  RESTRICTIONS:
%
%  METHOD:  Reduction to row-echelong form based on Gaussian elimination.
 
%  WRITTEN BY:  Mike Walker 	ON 	11/20/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
debug = 0;

%   n = counter of recursive calls to routine (debug variable only)
if(nargin<4)
   n = 0;
end

Amod = A;
n=n+1;

if(isempty(A))
   return;
end


% Step 1: find nonzero column
if(min(size(Amod))==1)
   ii = min(find(abs(Amod)));
else
   ii = min(find(max(abs(Amod))));
end
if(isempty(ii))
   return;
end;
%if(ii>1)
%   temp = Amod(:,1);
%   Amod(:,1) = Amod(:,ii);
%   Amod(:,ii) = temp;
%   temp = branchlist(1);
%   branchlist(1) = branchlist(ii);
%   branchlist(ii) = temp;
%end
%
% Step 2:
jj = min(find(Amod(:,ii)~=0));		% find nonzero row
if(debug)
   fprintf('n=%d, ii = %d, jj = %d\n',n,ii,jj);
end
if(jj>1)     % interchange rows
   temp = Amod(1,:);
   Amod(1,:) = Amod(jj,:);
   Amod(jj,:) = temp;
   temp = nodelist(1);
   nodelist(1) = nodelist(jj);
   nodelist(jj) = temp;
end

% Step 3:
Amod(1,:) = Amod(1,:)/Amod(1,ii);

% Step 4:
for i=2:size(Amod,1)
   Amod(i,:) = Amod(i,:) - Amod(i,ii)*Amod(1,:);
end

% Step 5:
if(debug)
   figure(1),clf,spy(Amod)
end
idx2 = setdiff(1:size(Amod,2),ii);
[Amod(2:end,idx2),nodelist(2:end),branchlist(idx2),n] = ...
	echelon(Amod(2:end,idx2),nodelist(2:end),branchlist(idx2),n);

n=n-1;
if(n==1)

% re-order to get identity in first half of matrix

   for(k=1:size(Amod,1))
      ii = min(find(abs(Amod(k,:))));
      if(ii>k)
if 0
         temp = Amod(:,k);
         Amod(:,k) = Amod(:,ii);
         Amod(:,ii) = temp;
         temp = branchlist(k);
         branchlist(k) = branchlist(ii);
         branchlist(ii) = temp;
end
         temp = Amod(:,[k:ii-1 ii+1:end]);
         Amod(:,k:end) = [Amod(:,ii) temp];
         temp = branchlist([k:ii-1 ii+1:end]);
         branchlist(k:end) = [branchlist(ii) temp];
      end
   end

% remove nonzero values above the diagonal

   for k=2:size(Amod,1)
      for j=1:k-1
         Amod(j,:) = Amod(j,:) - Amod(j,k)*Amod(k,:);
      end
   end

end
