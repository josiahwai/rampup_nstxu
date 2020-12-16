function [A,nodelist,sorted_netlist,Midx] = make_incidence(netlist)
 %
%  SYNTAX:  [A,nodelist,sorted_netlist,Midx] = make_incidence(netlist)
%
%  PURPOSE:  Construct incidence matrix.  Sort so that branches representing
%   voltage sources come first, then capacitors, then resistors.
%
%  INPUT:
%
%  OUTPUT:
%    A              = incidence matrix
%    nodelist       = set of all nodes specified in the netlist
%    sorted_netlist = netlist with elements sorted according to the order in [1]
%    Midx           = indices of inductive ('M') conductors in sorted_netlist
%			(in same order as original netlist)
%
%  RESTRICTIONS:
%
%  METHOD:  Described in [1] M.L. Walker, "A General Purpose Circuit Modeling 
%	Code for Tokamak Plasma Magnetic Control, General Atomics Engineering 
%	Physics Memo EPMmlw070120a, Jan. 20, 2007
 
%  WRITTEN BY:  Mike Walker 	ON 	11/20/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)make_incidence.m	1.3 04/19/11

debug=0;

struct_to_ws(netlist);

nodelist = union(N1,N2);
nodelist = sort(nodelist);
A = zeros(length(nodelist),length(N1));

% Rows of incidence matrix correspond to nodes in graph, columns to the elements 
% between nodes. A value of +1 means that current flows out of node and into 
% element, while -1 means current flows into the node from the element.  We use 
% update by +1 or -1 to allow for the case where both ends of element are 
% connected to the same node, so that net current flow from elt is 0.

Vidx=[]; Cidx=[]; Ridx=[]; Iidx=[]; Lidx=[]; Midx=[];
for inode = 1:length(nodelist)
   idx1 = find(N1==nodelist(inode));
   idx2 = find(N2==nodelist(inode));
if 1
   A(inode,idx1) = A(inode,idx1)-1;
   A(inode,idx2) = A(inode,idx2)+1;
else
   A(inode,idx1) = -1;
   A(inode,idx2) = 1;
end
end

A0=A;

% Process one element at a time.
% (When tested, this gives exactly the same as above, but is 
%  easier to understand.)

A = zeros(length(nodelist),length(N1));
for k=1:length(N1)
   A(N1(k)+1,k) = A(N1(k)+1,k) - 1;
   A(N2(k)+1,k) = A(N2(k)+1,k) + 1;
end

% Putting L before M here guarantees that if there is a linear dependence 
% between inductor currents, the L elements will go into the tree first and
% therefore the M elements will become states.

for ibranch=1:length(N1)
   if(strcmp(names(ibranch,1:1),'V'))
      Vidx = [Vidx ibranch];
   end
   if(strcmp(names(ibranch,1:1),'C'))
      Cidx = [Cidx ibranch];
   end
   if(strcmp(names(ibranch,1:1),'R'))
      Ridx = [Ridx ibranch];
   end
   if(strcmp(names(ibranch,1:1),'L'))
      Lidx = [Lidx ibranch];
   end
   if(strcmp(names(ibranch,1:1),'M'))
      Midx = [Midx ibranch];
   end
   if(strcmp(names(ibranch,1:1),'I'))
      Iidx = [Iidx ibranch];
   end
end

% Similarly, using reverse order for M objects means that earlier numbered 
% coils and vessel elements would preferentially be used as states.
Midx = fliplr(Midx);

idx = union(union(union(union(union(Vidx,Cidx),Ridx),Iidx),Lidx),Midx);
idx = setdiff(1:length(N1),idx);
if(~isempty(idx))
   wait(['ERROR make_incidence: invalid element: ' names(idx(1),:)])
end
idx_sort = [Vidx Cidx Ridx Lidx Midx Iidx];
Midx = Midx+length(Vidx)+length(Cidx)+length(Ridx)+length(Lidx);

A = A(:,idx_sort);
sorted_netlist.names = netlist.names(idx_sort,:);
sorted_netlist.N1 = netlist.N1(idx_sort);
sorted_netlist.N2 = netlist.N2(idx_sort);
sorted_netlist.values = netlist.values(idx_sort);

if(debug)
   fid = fopen('debug.dat','w');
   for k=1:length(sorted_netlist.N1)
      fprintf(fid,'%s %d %d %f\n',sorted_netlist.names(k,:), ...
	sorted_netlist.N1(k),sorted_netlist.N2(k),sorted_netlist.values(k));
   end
  fclose(fid);
end
