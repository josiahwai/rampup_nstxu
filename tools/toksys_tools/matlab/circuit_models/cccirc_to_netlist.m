function [netlist,inode1,inode2] = ...
 		 cccirc_to_netlist(cccirc,tok_data_struct,netlist_file);
 %
%  SYNTAX:  [netlist,inode1,inode2] = ...
%		 cccirc_to_netlist(cccirc,tok_data_struct,netlist_file)
%
%  PURPOSE:  Convert circuit connection defined by cccirc, vvgroup, and vvcirc
%		 into netlist.
%
%  INPUT:
%	cccirc = either cccirc vector or structure containing any combination
%			of cccirc, vvgroup, and vvcirc
%	tok_data_struct = vacuum data objects structure
%	netlist_file = name of netlist file to create (optional, not created
%			if no name specified)
%
%  OUTPUT:
%	netlist = similar to netlist input to Spice code, except that the 
%	   "value" of M entries is the row corresponding to M and R matrices for
%          that coil.  Fictitious connections to ground are added to simplify 
%	   the nodal analysis that occurs in building the system model.
% 
%  RESTRICTIONS: Code has been added to handle non-default vvgroup,vvcirc and the more
%	complex code has been verified to work for default vvgroup, vvcirc but not 
%	extensively tested for non-default settings yet.
 
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	11/19/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% What are outputs inode1, inode2 used for?

if(nargin==3)
   if(~ischar(netlist_file))
      wait('ERROR cccirc_to_netlist: netlist_file name must be string')
      return;
   end
   write_file=1;
else
   write_file=0;
end
if(isstruct(cccirc))
   struct_to_ws(cccirc);
end

if(exist('cccirc')~=1)
   cccirc = 1:tok_data_struct.nc;
end
if(exist('vvgroup')~=1)
   vvgroup = 1:tok_data_struct.nv;
end
if(exist('vvcirc')==1)
   temp = vvgroup(:,1);
   temp = setdiff(unique(temp),0);
   if(length(temp) ~= max(temp))
      wait('ERROR build_tokamak_system: vvgroup must define contiguous numbered groups')
      return;
   end
   temp1 = setdiff(union(temp,abs(vvcirc)),0);
   if(length(temp1)~=length(temp))
      wait('ERROR build_tokamak_system: vvcirc and vvgroup inconsistent')
      return;
   end
else
   vvcirc = setdiff(unique(vvgroup),0);
end
if(size(vvgroup,1)<=size(vvgroup,2))
   vvgroup = vvgroup';
end
if(any(vvgroup(:,1)<0))
   wait('ERROR cccirc_to_netlist: vvgroup has negative indices');
   return;
end

% Check consistency of inputs vvcirc and vvgroup:
idx_vvcirc = unique(setdiff(abs(vvcirc),0));
idx_vvgroup = unique(setdiff(vvgroup(:,1),0));
if(length(setdiff(idx_vvcirc,idx_vvgroup))>0)
   wait('ERROR cccirc_to_netlist: vvcirc references indices not defined by vvgroup')
   return;
elseif(length(setdiff(idx_vvgroup,idx_vvcirc))>0)
   wait('ERROR cccirc_to_netlist: vvgroup defines indices not referenced by vvcirc')
   return;
end

% Key idea here is to use a "coil" element in the same way that resistors,
% inductors, capacitors, etc. are used, but with "value" defined by a row
% in the circuit equation model.

% Put the netlist processing as a layer below the "cccirc" specification, so
% that easy to model tokamaks can still use the cccirc, while more difficult
% ones need to get into the details of the netlist.  (Can they be optimally
% combined?)

% Step 1: use standard initial coil modeling to build M and R

M = [	tok_data_struct.mcc tok_data_struct.mcv; ...
	 tok_data_struct.mcv' tok_data_struct.mvv];
R = diag([tok_data_struct.resc; tok_data_struct.resv]);
nc = tok_data_struct.nc;
nv = tok_data_struct.nv;

% Step 2: in netlist, need to define nodes at either end of coil
%	(here is where connected coils get defined by assigning same node
%	number to ends of different coils)
% Use M# to indicate an "inductor component"

if(isfield(tok_data_struct,'ecdata') & ~isempty(tok_data_struct.ecdata))
   coil_names = strvcat(tok_data_struct.ecnames, tok_data_struct.fcnames);
else
   coil_names = tok_data_struct.fcnames;
end

ncirc = length(unique(abs(setdiff(cccirc,0))));	% number of coil circuits
if(length(cccirc) > size(coil_names,1))
   wait('ERROR cccirc_to_netlist: length(cccirc) > # coils defined in ecnames+fcnames')
   return;
end
inode1(1) = 0; inode2(1) = 1;
row = 0;	% element counter

% Make sure to put current that is desired to represent the state last in
% the netlist, since that will ensure it ends up in the co-tree.

netlist = ''; Vnetlist = '';
maxnode = 0;
for k=1:ncirc
   idx = find(abs(cccirc)==k);
   if(~isempty(idx) & all(cccirc(idx)<0))
      wait(['ERROR cccirc_to_netlist: circuit ' int2str(k) ...
					' has no positive indices']);
      return;
   end
   if(length(idx)>0)
      row = row+1;	% increment each time a line is added to netlist
      startrow = row;
      staterow = 0;	% reset index

      ii = idx(1);
      if(cccirc(ii)>0)
         inode1(row) = 0; 
         inode2(row) = maxnode+1;	%new node
         staterow = row;
      else
         inode1(row) = maxnode+1;	%new node
         inode2(row) = 0; 
      end
      maxnode = maxnode+1;
      node1str = int2str(inode1(row));
      node2str = int2str(inode2(row));
      kstr = int2str(ii);
      Mstr = ['M' coil_names(ii,:)];
      str = [Mstr ' ' node1str ' ' node2str ' ' kstr '     % PF coil ' kstr];
      netlist = strvcat(netlist,str);
      for j=2:length(idx)
         row = row+1;	% increment each time a line is added to netlist
         i1 = idx(j);
         if(cccirc(i1)>0)
            inode1(row) = maxnode;
            inode2(row) = maxnode+1;
            if(~staterow), staterow = row; end;
         else
            inode1(row) = maxnode+1;
            inode2(row) = maxnode;
         end
         maxnode = maxnode+1;
         node1str = int2str(inode1(row));
         node2str = int2str(inode2(row));
         kstr = int2str(i1);
         Mstr = ['M' coil_names(i1,:)];
         str = [Mstr ' ' node1str ' ' node2str ' ' kstr '     % PF coil ' kstr];
         netlist = strvcat(netlist,str);
      end
      if(length(idx)>1)		% re-sort to put current state last
         temp = netlist(startrow:row,:);
         rowidx = setdiff(startrow:row,staterow) - startrow + 1;
         netlist(startrow:row-1,:) = temp(rowidx,:);
         netlist(row,:) = temp(staterow-startrow+1,:);
      end

      vnode1(k) = 0;
      vnode2(k) = maxnode;
      node1str = int2str(vnode1(k));
      node2str = int2str(vnode2(k));
      kstr = int2str(k);
      str=['V' kstr ' ' node1str ' ' node2str ' ' ...
				 kstr '     % PF voltage ' kstr];
      Vnetlist = strvcat(Vnetlist,str);
   else
      wait('ERROR cccirc_to_netlist: cccirc is missing an index')
      return;
   end

end

netlist = strvcat(netlist,Vnetlist); 
row = row + size(Vnetlist,1);		% total number of lines in netlist

% Close all vessel coils at the ground to make later network analysis easy.

if 1 	% If vvgroup and/or vvcirc exist:

nvvgroup = length(unique(setdiff(vvgroup(:,1),0)));

nvvcirc = length(unique(abs(setdiff(vvcirc,0))));	% number of vessel circuits

for k=1:nvvcirc
   idx = find(abs(vvcirc)==k);
   if(~isempty(idx) & all(vvcirc(idx)<0))
      wait(['ERROR cccirc_to_netlist: circuit ' int2str(k) ...
					' has no positive indices']);
      return;
   end
   if(length(idx)>0)
      ii = idx(1);
      if(length(idx)>1)
         if(vvcirc(ii)>0)
            start_node = 0; 
            end_node = maxnode+1;	%new node
         else
            start_node = maxnode+1; 
            end_node = 0;
         end
         maxnode = maxnode+1;
      else
         start_node = 0; 
         end_node = 0;	
      end

% loop over indices in vvgroup that match:
      idxvvgroup = find(vvgroup(:,1)==ii);
      for j=1:length(idxvvgroup)
         row = row+1;
         k = idxvvgroup(j);
         inode1(row) = start_node; 
         inode2(row) = end_node;
         node1str = int2str(inode1(row));
         node2str = int2str(inode2(row));
         kstr1 = int2str(k+nc); kstr2 = int2str(k);
         Mstr = ['Mvv' int2str(k)];
         str = [Mstr ' ' node1str ' ' node2str ' ' kstr1 '    % vessel element ' kstr2];
         netlist = strvcat(netlist,str);
      end

      for j=2:length(idx)
         i1 = idx(j);
         if(j<length(idx))
            if(vvcirc(i1)>0)
               start_node = maxnode;
               end_node = maxnode+1;	% new node
            else
               start_node = maxnode+1;
               end_node = maxnode;
            end
            maxnode = maxnode+1;
         else			% special case for last element in series
            if(vvcirc(i1)>0)
               start_node = maxnode; 
               end_node = 0;	
            else
               start_node = 0; 
               end_node = maxnode;	
            end
         end

% loop over indices in vvgroup that match:
         idxvvgroup = find(vvgroup(:,1)==i1);
         for m=1:length(idxvvgroup)
            row = row+1;
            k = idxvvgroup(m);
            inode1(row) = start_node; 
            inode2(row) = end_node;
            node1str = int2str(inode1(row));
            node2str = int2str(inode2(row));
            kstr1 = int2str(k+nc); kstr2 = int2str(k);
            Mstr = ['Mvv' int2str(k)];
            str = [Mstr ' ' node1str ' ' node2str ' ' kstr1 '     % vessel element ' kstr2];
            netlist = strvcat(netlist,str);
         end
      end

   else
      wait('ERROR cccirc_to_netlist: vvcirc is missing an index')
      return;
   end
end

else		% otherwise, if all vessel connections are default

for k=1:tok_data_struct.nv
   row = row+1;
   inode1(row) = 0;
   inode2(row) = 0;
   node1str = int2str(inode1(row));
   node2str = int2str(inode2(row));
   kstr = int2str(k+nc);
   Mstr = ['Mvv' int2str(k)];
   str=[Mstr ' ' node1str ' ' node2str ' ' kstr '    % vessel ' int2str(k)];
   netlist = strvcat(netlist,str);
end

end

if(write_file)
   fid = fopen(netlist_file,'w');
   for k=1:size(netlist,1)
      fprintf(fid,'%s\n',netlist(k,:));
   end
   fclose(fid);
end

