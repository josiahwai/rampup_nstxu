function model = model_from_netlist(netlist,Mmat,Rvec,Iout,Vout,verbose)
 %
%  SYNTAX: model = model_from_netlist(netlist,Mmat,Rvec,Iout,Vout,verbose)
%
%  PURPOSE: 
%   Produces model that is generalization of that produced by the simpler 
%   tokamak/plasma modeling schemes, of the form:
%       Mhat*dx/dt + Rhat*x = Vhat*u + What*du/dt
%                         y = Chat*x + Dhat*u + DWhat*du/dt
%   where the generalized state and input vectors are:
%	x = [i_inductors; v_capacitors],
%	u = [Vsource; Isource]	
%	y = [current_outputs; voltage_outputs]
%
%  INPUT:
%    netlist = netlist specification of circuit connections
%    Mmat = mutual inductance matrix for all "M" type elements in netlist
%    Rvec = resistance vector for all "M" type elements in netlist
%     (Values of netlist for "M" objects are indices into these objects.)
%    Iout = array of strings defining branches for which to output currents
%		(optional, default = no current outputs)
%    Vout = n x 2 matrix, each row containing node numbers [N1 N2], with
%		voltage output defined as V(N1)-V(N2)
%		(optional, default = no voltage outputs)
%    verbose = set to >0 to print model-building diagnostics to terminal
%		(optional, default = 0)
%
%  OUTPUT: (in structure "model")
%    model = structure containing model objects; see descriptions field for
%		description of contents of each field
%
%  RESTRICTIONS: 
%       (1) Does not yet handle circuits containing capacitors or current sources.
%	(2) Right now it can't handle outputs Vout if circuit is not specified as
%		a connected graph (i.e. no floating potentials).
%
%  METHOD: Described in: M.L.Walker, A General Purpose Circuit Modeling Code 
%     for Tokamak Plasma Magnetic Control, Engineering Physics Memo 
%     EPM070120a, January 20, 2007

%  WRITTEN BY:  Mike Walker 	ON 	12/28/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nargin < 6)
   verbose = 0;
end

%TO DO:
% figure out where to put current source inputs, if they exist.
% fix method to produce 0 columns in Vhat
% is resistance in inductive elements handled correctly?
% fix input order to [V,I]
% what is Iout_mat?
% what is Icond_idx?
% add names of input signals to output of function
% need to confirm that Qbe (Wilson) is same as Rohrer calculation
% What are the "exception" conditions?  Do we handle them correctly? (or at
%  least flag them?)  Will they ever occur?
% If there are only two branches on a node, add logic to treat as one branch
%  (for inductors and resistors).
% Sort resistors for improved numerics.

%%%%%%%%%%%%%%%%%%
% Error checking:
%%%%%%%%%%%%%%%%%%

for k=1:size(netlist.names,1)
   name = netlist.names(k,:);
   idx = strmatch(name,netlist.names);
   if(length(idx)>1)
      fprintf(['ERROR model_from_netlist: Netlist elts in rows %d,' ...
      '%d have the same name "%s". \n'],idx(1),idx(2),deblank(name));
      wait('Names must be unique.');
      return;
   end
end

if(nargin<5)
   Vout = [];
elseif(~isempty(Vout))
   if(size(Vout,2)~=2)
      wait('ERROR model_from_netlist: Vout must be size n x 2')
      return;
   end
end
if(nargin<4)
   Iout = [];
end

[s1,s2]=size(Mmat);
s3=length(Rvec);
if(s1~=s2)
   wait('ERROR model_from_netlist: Mmat must be square')
   return;
end
if(s1~=s3)
   wait('ERROR model_from_netlist: length(Rvec) must match size(Mmat)')
   return;
end

Midx = find(netlist.names(:,1:1)=='M')';
[mi,ii]=max(netlist.values(Midx));
if(mi>s3)
   fprintf('ERROR model_from_netlist: invalid value %d for netlist element %s.', ...
		mi, netlist.names(Midx(ii),:));
   wait('   Value must be the index of a row/column of Mmat.')
   return;
end
if(length(unique(netlist.values(Midx))) < length(Midx))
   srtd =  sort(netlist.values(Midx));
   dsrtd = diff(srtd);
   ii = srtd(find(dsrtd==0));
   fprintf('\nERROR model_from_netlist: value %d is repeated for multiple netlist "M" elements.\n', ii);
   wait
   return;
end

idx = find(netlist.names(:,1)=='V');
if(length(idx)~=length(unique(netlist.values(idx))))
   wait('ERROR model_from_netlist: duplicate voltage source indices')
   return;
end

%%%%%%%%%%%%%%%%%%%%%%
% Begin the real work.
%%%%%%%%%%%%%%%%%%%%%%

Rvec0 = Rvec;	% save original content of Rvec

% index vectors identifying classes of elements
%   alpha   = capacitive chords
%   beta    = resistive chords
%   gamma   = inductive chords
%   delta   = capacitive branches
%   epsilon = resistive branches
%   zeta    = inductive branches
%   phi     = current source chords
%   theta   = voltage source branches
% chords is terminology for branches contained in the co-tree, while branches
% in the tree are still called branches.

% Before processing the netlist, see if any resistances can be combined into
% "inductive" elements.

netlist_change=1;
while(netlist_change)
   netlist_change=0;

% find resistor element
   idxR = find(netlist.names(:,1:1)=='R');
   idx_remove = [];	% indices of elements to remove from netlist

% See if 1st or 2nd node of resistor appears in another element = inductor or resistor
% If so, it can be combined with that element.

   for k=1:length(idxR)
      idxr = idxR(k);
      N1r = netlist.N1(idxr);
      N2r = netlist.N2(idxr);
      idx1 = union(find(netlist.N1==N1r),find(netlist.N2==N1r));
      idx2 = union(find(netlist.N1==N2r),find(netlist.N2==N2r));
      idx1 = setdiff(idx1,idxr);
      idx2 = setdiff(idx2,idxr);

      type1 = netlist.names(idx1,1:1);
      type2 = netlist.names(idx2,1:1);
      if(length(idx1)==1 & ismember(type1,['R','L','M']))
         if(type1=='R')
            netlist.values(idx1)=netlist.values(idx1)+netlist.values(idxr);
         elseif(type1=='L')
            netlist.values(idx1)=netlist.values(idx1)+netlist.values(idxr)*i;
         elseif(type1=='M')
            Rvec(netlist.values(idx1)) = ...
		 Rvec(netlist.values(idx1)) + netlist.values(idxr);
         end
         if(netlist.N1(idx1)==N1r)
            netlist.N1(idx1)=N2r;
         elseif(netlist.N1(idx1)==N2r)
            netlist.N1(idx1)=N1r;
         elseif(netlist.N2(idx1)==N1r)
            netlist.N2(idx1)=N2r;
         elseif(netlist.N2(idx1)==N2r)
            netlist.N2(idx1)=N1r;
         end
         idx_remove = [idx_remove idxr];
         netlist_change=1;
         break;
      elseif(length(idx2)==1 & ismember(type2,['R','L','M']))
         if(type2=='R')
            netlist.values(idx2)=netlist.values(idx2)+netlist.values(idxr);
         elseif(type2=='L')
            netlist.values(idx2)=netlist.values(idx2)+netlist.values(idxr)*i;
         elseif(type2=='M')
            Rvec(netlist.values(idx2)) = ...
		 Rvec(netlist.values(idx2)) + netlist.values(idxr);
         end
         if(netlist.N1(idx2)==N1r)
            netlist.N1(idx2)=N2r;
         elseif(netlist.N1(idx2)==N2r)
            netlist.N1(idx2)=N1r;
         elseif(netlist.N2(idx2)==N1r)
            netlist.N2(idx2)=N2r;
         elseif(netlist.N2(idx2)==N2r)
            netlist.N2(idx2)=N1r;
         end
         idx_remove = [idx_remove idxr];
         netlist_change=1;
         break;
      end
   end

   idx = setdiff(1:length(netlist.N1),idx_remove);
   netlist.names = netlist.names(idx,:);
   netlist.N1 = netlist.N1(idx);
   netlist.N2 = netlist.N2(idx);
   netlist.values = netlist.values(idx);
end

allnodes = union(netlist.N1,netlist.N2);
numNode = length(allnodes);
numV = length(find(netlist.names(:,1)=='V'));

if 1 	% TESTING 1

% All logic is easier if the set of nodes is contiguous and starts at 0.
netlist0 = netlist;
[netlist1,orig_nodes] = make_nodes_contiguous(netlist);

% Need to modify Vout, which requests voltages in terms of original nodes.
% At this point, orig_nodes contains an ordered list of all original nodes.
% These are converted by make_nodes_contiguous into a contiguous list of 
% nodes starting at 0, so that the order of the nodes is unchanged.

Vout_mod = Vout;
for k=1:length(Vout(:))
   idx = find(orig_nodes==Vout(k));
   if(isempty(idx))
      fprintf(['ERROR model_from_netlist: output voltage requested relative ' ...
	'to node %d, which does not exist\n'],Vout(k));
      wait;
      return;
   end
   Vout_mod(k) = idx-1;
end

allnodes = union(netlist1.N1,netlist1.N2);
numNode = length(allnodes);
numV = length(find(netlist.names(:,1)=='V'));

if 1	% TESTING 2

% If necessary, modify netlist so that it represents a connected graph.
% The way we do this is to require every subgraph to touch the node 0 (often 
% interpreted as "ground").  If there is a subgraph that is unconnected to 
% node 0, just modify one of nodes of that subgraph to be node 0.

nodelist = allnodes;

cost_matrix = Inf*ones(numNode+1);
for k=1:length(netlist1.N1)
   cost_matrix(netlist1.N1(k)+1,netlist1.N2(k)+1)=1;
   cost_matrix(netlist1.N2(k)+1,netlist1.N1(k)+1)=1;
end

netlist_change=1;
while(netlist_change)
   netlist_change=0;

   numNode=length(nodelist);
   for k=1:numNode
      nodenum = nodelist(k);
      [path,cost]=dijkstra(numNode,cost_matrix,1,nodenum+1);
      if(isinf(cost))
         idx1 = find(netlist1.N1==nodenum);
         idx2 = find(netlist1.N2==nodenum);
         netlist1.N1(idx1) = 0;
         netlist1.N2(idx2) = 0;
         nodelist = setdiff(nodelist,nodenum);
         numNode = numNode-1;
         netlist_change=1;
         break;
      end
   end
   netlist1 = make_nodes_contiguous(netlist1);
   nodelist = union(netlist1.N1,netlist1.N2);
   cost_matrix = Inf*ones(numNode+1);
   for k=1:length(netlist1.N1)
      cost_matrix(netlist1.N1(k)+1,netlist1.N2(k)+1)=1;
      cost_matrix(netlist1.N2(k)+1,netlist1.N1(k)+1)=1;
   end
end

allnodes = nodelist;

end 	% TESTING 2

netlist = netlist1;

end 	% TESTING 1

% Pre-sort mutual inductance elements so that they appear in ascending order
% of value (corresponding to indices in M, R).
Icond_idx = strmatch('M',netlist.names);
[temp,idx] = sort(netlist.values(Icond_idx));
temp = netlist; 
temp.names(Icond_idx,:)  = netlist.names(Icond_idx(idx),:);
temp.N1(Icond_idx)     = netlist.N1(Icond_idx(idx));
temp.N2(Icond_idx)     = netlist.N2(Icond_idx(idx));
temp.values(Icond_idx) = netlist.values(Icond_idx(idx));
netlist = temp;

% Pre-sort voltage source indices so that they appear in ascending order
% of value (corresponding to voltage source index defined by user).
Vidx = find(netlist.names(:,1:1)=='V');
[temp,idx] = sort(netlist.values(Vidx));
temp = netlist; 
temp.names(Vidx,:)  = netlist.names(Vidx(idx),:);
temp.N1(Vidx)     = netlist.N1(Vidx(idx));
temp.N2(Vidx)     = netlist.N2(Vidx(idx));
temp.values(Vidx) = netlist.values(Vidx(idx));
netlist = temp;

% Icond_names = netlist.names(Icond_idx,:);

[A,nodelist,sorted_netlist,Icond_idx] = make_incidence(netlist);
A = A(2:end,:);		% Why do we do this???

%   At this point, nodelist defines the nodes corresponding to each row, and 
% sorted_netlist defines the branches corresponding to each column of A, after
% sorting to put in order: voltage sources, capacitors, resistors, inductors,
% current sources, with inductors sorted in reverse order using the value (index
% in M,R) so that the selection of inductor branches will come from lower ordered
% indices.
%   The sorted_netlist object is a valid "netlist" that represents the user specified
% circuit connections, so can be used as the starting point for what follows.
%   By construction, all conductors that close on themselves will produce
% columns with only zeros in A.

% Now apply the echelon algorithm to define the tree and co-tree ...

[Amod,nodelist,branchlist] = echelon(A,nodelist,[1:size(A,2)],0);
nbranches = length(branchlist);

% At this point, Amod is of the form [I -F'], branchlist contains the 
% original (branch) indices of the re-ordered columns of A, and nodelist 
% defines the node numbers corresponding to each row of Amod.

tree_len = size(Amod,1);	% (not counting node 0)
tree_idx = branchlist(1:tree_len);

if(verbose>1)
  disp('Branches of tree are:')
  sorted_netlist.names(tree_idx,:)
end

% The first #nodes-1 branches are in the tree. The remainder are in the co-tree.
Atree = Amod(:,1:tree_len);
if(norm(Atree-eye(size(Atree)))~=0)
   wait('ERROR : Atree should be identity')
   return;
end
Acotree = Amod(:,tree_len+1:end);

% Sort netlist to put circuit elements (branches) in correspondence with the
% columns of the incidence matrix

cotree_idx = setdiff(1:size(A,2),tree_idx);
netlist2.names = sorted_netlist.names([tree_idx cotree_idx],:);
netlist2.N1 = sorted_netlist.N1([tree_idx cotree_idx]);
netlist2.N2 = sorted_netlist.N2([tree_idx cotree_idx]);
netlist2.values = sorted_netlist.values([tree_idx cotree_idx]);

% At this point, the tree should have all of the voltage sources, a few 
% resistors, and maybe some inductors, in that order.  The co-tree
% should have (maybe) some resistors, and the rest of the inductors,
% in that order (assumes no capacitors or current sources).

% Check that all voltage sources show up as first n branches (must be in the tree):

if(any(netlist2.names(1:numV,1:1)~='V'))
   wait('ERROR model_from_netlist: voltage sources not in tree')
end

% Partition for calculation.  Make sure to preserve the correspondence between 
% elements in netlist2 and the columns of Amod.

rvector = zeros(nbranches,1);		% length = # circuit elements
Ridx = find(netlist2.names(:,1)=='R');
rvector(Ridx) = netlist2.values(Ridx);

epsilon = Ridx(find(Ridx<=tree_len)); 	% define resistor (R elts) in tree
neps = length(epsilon);	

beta = Ridx(find(Ridx>tree_len)); 	% define resistor (R elts) in co-tree
nbeta = length(beta);

L = zeros(nbranches,nbranches);		% size = # circuit elements

% Coil/vessel conductors (M) have both resistance and inductance, defined
% by the vector Rvec and matrix Mmat.

Midx = find(netlist2.names(:,1:1)=='M')';
if(length(Midx)>length(Rvec))
   fprintf('ERROR model_from_netlist: # netlist M entries=%d > size Rvec,Mmat=%d\n', ...
			length(Midx), length(Rvec))
   wait
   return;
end
idx = netlist2.values(Midx);
rvector(Midx) = rvector(Midx) + Rvec(idx);
L(Midx,Midx) = Mmat(idx,idx);
num_M_conds = size(Mmat,2);
if(num_M_conds < length(Midx))
   wait('ERROR model_from_netlist: # of mutual (M) elements in netlist > size of input M matrix\n');
   return;
elseif(num_M_conds > length(Midx))
   % This is to remind users to remember to handle the Ip entry correctly:
   disp('WARNING model_from_netlist: size of input M matrix > # of mutual (M) elements in netlist\n');
   wait('  Make sure that Ip entry in netlist points to last row in M matrix. ')
end

% Isolated inductors (L) also have both inductance and resistance, but no
% coupling to other inductors:

Lidx = find(netlist2.names(:,1:1)=='L');
for k=1:length(Lidx)
   idx = Lidx(k);
   L(idx,idx) = imag(netlist2.values(idx));
   R(idx,idx) = real(netlist2.values(idx));
end

% Collect all inductors:

LMidx = union(Lidx,Midx);

zeta = LMidx(find(LMidx<=tree_len)); 	% define inductors in tree
nzeta = length(zeta);

gamma = LMidx(find(LMidx>tree_len)); 	% define inductors in co-tree
ngamma = length(gamma);

% Capacitors:

C = zeros(nbranches,1);
Cidx = find(netlist2.names(:,1:1)=='C');
C(Cidx) = netlist2.values(Cidx);

delta = Cidx(find(Cidx<=tree_len));	% define capacitors in tree
ndelta = length(delta);

alpha = Cidx(find(Cidx>tree_len)); 	% define capacitors in co-tree
nalpha= length(alpha);

% Here, the objects rvector, L, and C all have dimension equal to the set of
% all branches.  Index into the appropriate object by type of branch.

% Inductors:

Iidx = find(netlist2.names(:,1:1)=='I');

phi = Iidx(find(Iidx>tree_len)); 	% current sources in co-tree
nphi = length(phi);
if(nphi ~= length(Iidx))
   wait('ERROR model_from_netlist: not all current sources in co-tree')
   return;
end

Vidx = find(netlist2.names(:,1:1)=='V');
V_original_values = netlist2.values(Vidx);

% code here is PROBABLY OBSOLETE:
% Make indices for voltage sources contiguous:
% ********** CHECK **********
%   Shouldn't this rearrange other fields in netlist2 as well?
% ********** CHECK **********
for k=1:numV
   change = 1;
   while(change)	% if an index is missing shift down those above it
      change = 0;
      idx = find(netlist2.values(Vidx)==k);
      if(isempty(idx))
         netlist2.values(Vidx(k:end)) = netlist2.values(Vidx(k:end))-1;
         change = 1;
      end
   end
end

theta = Vidx(find(Vidx<=tree_len)); 	% voltage sources in tree
ntheta   = length(theta);
if(ntheta ~= length(Vidx))
   wait('ERROR model_from_netlist: not all voltage sources in tree')
   return;
end
if(max(netlist2.values(Vidx))~=ntheta)
   wait('ERROR model_from_netlist: error in voltage source values')
   return;
end

% At this point, the tree and co-tree have been defined so we can fix the order
% of elements pointed to by gamma and zeta back to ascending order.  Since the
% element list, which is maintained by netlist2, is still in reverse order, we can
% just reverse the order of elements in the gamma and zeta lists to point to them
% in reverse order.

gamma = flipud(gamma); zeta = flipud(zeta);

F = -Acotree';
Fat = F(alpha-tree_len,theta);
Fad = F(alpha-tree_len,delta);
Fae = F(alpha-tree_len,epsilon);
Faz = F(alpha-tree_len,zeta);

Fbt = F(beta-tree_len,theta);
Fbd = F(beta-tree_len,delta);
Fbe = F(beta-tree_len,epsilon);
Fbz = F(beta-tree_len,zeta);

Fgt = F(gamma-tree_len,theta);
Fgd = F(gamma-tree_len,delta);
Fge = F(gamma-tree_len,epsilon);
Fgz = F(gamma-tree_len,zeta);

Fpt = F(phi-tree_len,theta);
Fpd = F(phi-tree_len,delta);
Fpe = F(phi-tree_len,epsilon);
Fpz = F(phi-tree_len,zeta);

Rbb = diag(rvector(beta));	% separate resistors in co-tree
Rgg = diag(rvector(gamma));	% resistances of state inductors (co-tree)
Rzz = diag(rvector(zeta));	% resistances of non-state (in tree)
Ree = diag(rvector(epsilon));	% separate resistors in tree

Caa = C(alpha,alpha);
Cdd = C(delta,delta);

Lgg = L(gamma,gamma);
Lgz = L(gamma,zeta);
Lzg = L(zeta,gamma);
Lzz = L(zeta,zeta);

Qgg=[eye(ngamma) Fgz]*[Lgg Lgz; Lzg Lzz]*[eye(ngamma); Fgz'];
Qbb = Rbb + Fbe*Ree*Fbe';
Qbbinv = inv(Qbb);
Qbe = Qbbinv*Fbe*Ree;

Pdd = Cdd + Fad'*Caa*Fad;
Rbbinv = inv(Rbb);
Pee = inv(Ree) + Fbe'*Rbbinv*Fbe;
Peeinv = inv(Pee);

% Generalization of the usual tokamak model structure:

if(verbose>0)
fprintf('system description:\n');
fprintf('%d current states = ',ngamma);
for k=1:ngamma
   fprintf('%s ',netlist2.names(gamma(k),:)); 
end
fprintf('\n');
fprintf('%d voltage states = ',ndelta);
for k=1:ndelta
   fprintf('%s ',netlist2.names(delta(k),:)); 
end
fprintf('\n');
fprintf('%d current sources = ',nphi);
for k=1:nphi
   fprintf('%s ',netlist2.names(phi(k),:)); 
end
fprintf('\n');
fprintf('%d voltage sources = ',ntheta);
for k=1:ntheta
   fprintf('%s ',netlist2.names(theta(k),:)); 
end
fprintf('\n');
end

Mhat = [Qgg zeros(ngamma,ndelta);
	zeros(ndelta,ngamma) Pdd];

R11 = Rgg + Fgz*Rzz*Fgz' + Fge*Peeinv*Fge';
R12 = Fgd - Fge*Qbe'*Fbd;	% different from Rohrer
R21 = -Fgd' + Fbd'*Qbe*Fge';	% different from Rohrer
R22 = Fbd'*Qbbinv*Fbd;
Rhat = [R11 R12; R21 R22]; 

% matrix that multiplies inputs i_phi and v_theta

V11 = -Fgz*Rzz*Fpz' - Fge*Peeinv*Fpe';
V12 = -Fgt + Fge*Qbe'*Fbt;
V21 = Fpd' - Fbd'*Qbe*Fpe';
V22 = -Fbd'*Qbbinv*Fbt;
Vhat = [V11 V12; V21 V22];

% matrix that multiplies derivatives of inputs i_phi and v_theta

What = [-(Lgz+Fgz*Lzz)*Fpz' zeros(ngamma,ntheta);
	 zeros(ndelta,nphi) -Fad'*Caa*Fat];

if 0  % Seems to be NOT USED
% Construct mapping to get all 'M' currents from states.
% Indices gamma represent all currents (states) in co-tree, and zeta
% represents all currents (non-states) in tree.

nIout = size(Icond_names,1);
Iout_mat = zeros(nIout,length(branchlist)-tree_len);
for k=1:nIout
   Icond_names(k,:);
   idx = strmatch(deblank(Icond_names(k,:)),netlist2.names,'exact');
%    wait
   if(idx>tree_len)
      Iout_mat(k,idx-tree_len)=1;
   else
      Iout_mat(k,:)= -Acotree(idx,:);
   end
%   Iout_mat(:,1:35)
%   wait
end

if(verbose>0)
fprintf('%d current outputs = ',nIout);
for k=1:nIout
   fprintf('%s ',Icond_names(k,:)); 
end
fprintf('\n');
end

end 	% if 0

nIstates = length(gamma);
P = [eye(nIstates); Fgz'];

% At this point, the indices in gamma define the current states, and indices in zeta define
% the currents that are linear combinations of those states.  The "hat" objects have these
% currents in the order [gamma, zeta].  For those currents that correspond to "M" type inductors,
% need to map back to the order used in the original M and R matrices.

% Istates_idx = netlist2.values(gamma);	% states only
% Icond_idx = netlist2.values([gamma; zeta]);	% all toroidal conductor currents

% Construct Pxx = mapping from states to full set of inductor currents/cap. voltages', etc:
Pxx = zeros(size(Mmat,1),ngamma);
for k=1:length(gamma)
   Pxx(netlist2.values(gamma(k)),k) = 1;
end
FgzT = Fgz'; % Dependence of non-states on the states is contained in Fgz.
idx=0;
for k=1:length(zeta)
   idx = idx+1;
   Pxx(netlist2.values(zeta(k)),:) = FgzT(idx,:);
end

% Sort the inputs according to the indices specified for each current and voltage
% source in the input netlist. Missing values will produce zero columns in Vhat.

temp = Vhat; Vhat = zeros(size(Vhat));
Vhat(:,Iidx) = temp(:,1:nphi);
Vhat(:,netlist2.values(Vidx)) = temp(:,nphi+[1:ntheta]);

if(length(Midx)==length(LMidx))		% if all inductors are "M" type ...
   test = norm(Mhat - Pxx'*Mmat*Pxx);
   if(test~=0)
      wait('ERROR model_from_netlist: failed test ...')
      return
   end
end

Istates = netlist2.names(gamma,:);
Istate_idx = netlist2.values(gamma);
len = size(netlist2.names,2);
inputs = char(zeros(nphi+ntheta,len));
inputs(Iidx,:)= netlist2.names(phi,:);
inputs(netlist2.values(Vidx),:) = netlist2.names(theta,:);

% Generate output matrices, if specified.
Chat = []; Dhat = []; dDhat = []; outputs='';
if(~isempty(Iout) | ~isempty(Vout))
   Mhatinv = inv(Mhat);
end
if(~isempty(Iout))
   for k=1:size(Iout,1);
      outputs = strvcat(outputs,['current through ' Iout(k,:)]);
   end
   temp = [zeros(nalpha,ngamma) -Fad'*Caa*Fad];
   CIa = -temp*Mhatinv*Rhat; 
   DIa = temp*Mhatinv*Vhat;
   dDIa= temp*Mhatinv*What + [zeros(nalpha,nphi) -Fad'*Caa*Fat];

   CIb = -[Qbe*Fge' Qbbinv*Fbd]; DIb = -[Qbe*Fpe' Qbbinv*Fbt];
   dDIb = zeros(nbeta,nphi+ntheta);

   CIg = eye(ngamma,ngamma+ndelta); DIg = zeros(ngamma,nphi+ntheta);
   dDIg = zeros(ngamma,nphi+ntheta);

   CIp = zeros(nphi,ngamma+ndelta); DIp = eye(nphi,nphi+ntheta); 
   dDIp = zeros(nphi,nphi+ntheta);

   CIt = Fat'*CIa + Fbt'*CIb + Fgt'*CIg + Fpt'*CIp;
   DIt = Fat'*DIa + Fbt'*DIb + Fgt'*DIg + Fpt'*DIp;
   dDIt = Fat'*dDIa + Fbt'*dDIb + Fgt'*dDIg + Fpt'*dDIp;

   CId = Fad'*CIa + Fbd'*CIb + Fgd'*CIg + Fpd'*CIp; 
   DId = Fad'*DIa + Fbd'*DIb + Fgd'*DIg + Fpd'*DIp; 
   dDId = Fad'*dDIa + Fbd'*dDIb + Fgd'*dDIg + Fpd'*dDIp; 

   CIe = Fbe'*CIb + Fge'*CIg + Fpe'*CIp;
   DIe = Fbe'*DIb + Fge'*DIg + Fpe'*DIp;
   dDIe = Fbe'*dDIb + Fge'*dDIg + Fpe'*dDIp;

   CIz = Fgz'*CIg + Fpz'*CIp;
   DIz = Fgz'*DIg + Fpz'*DIp;
   dDIz = Fgz'*dDIg + Fpz'*dDIp;

   j=0;	% KLUGE to get past matlab
   for k=1:size(Iout,1)
      Iout_index = strmatch(Iout(k,:),netlist2.names);
      if(ismember(Iout_index, alpha))
         j = find(Iout_index==alpha);
         Chat=[Chat;CIa(j,:)]; Dhat=[Dhat;DIa(j,:)]; dDhat=[dDhat;dDIa(j,:)];
      elseif(ismember(Iout_index, beta))
         j = find(Iout_index==beta);
         Chat=[Chat;CIb(j,:)]; Dhat=[Dhat;DIb(j,:)]; dDhat=[dDhat;dDIb(j,:)];
      elseif(ismember(Iout_index, gamma))
         j = find(Iout_index==gamma);
         Chat=[Chat;CIg(j,:)]; Dhat=[Dhat;DIg(j,:)]; dDhat=[dDhat;dDIg(j,:)];
      elseif(ismember(Iout_index, phi))
         j = find(Iout_index==phi);
         Chat=[Chat;CIp(j,:)]; Dhat=[Dhat;DIp(j,:)]; dDhat=[dDhat;dDIp(j,:)];
      elseif(ismember(Iout_index, theta))
         j = find(Iout_index==theta);
         Chat=[Chat;CIt(j,:)]; Dhat=[Dhat;DIt(j,:)]; dDhat=[dDhat;dDIt(j,:)];
      elseif(ismember(Iout_index, delta))
         j = find(Iout_index==delta);
         Chat=[Chat;CId(j,:)]; Dhat=[Dhat;DId(j,:)]; dDhat=[dDhat;dDId(j,:)];
      elseif(ismember(Iout_index, epsilon))
         j = find(Iout_index==epsilon);
         Chat=[Chat;CIe(j,:)]; Dhat=[Dhat;DIe(j,:)]; dDhat=[dDhat;dDIe(j,:)];
      elseif(ismember(Iout_index, zeta))
         j = find(Iout_index==zeta);
         Chat=[Chat;CIz(j,:)]; Dhat=[Dhat;DIz(j,:)]; dDhat=[dDhat;dDIz(j,:)];
      else
         wait(['ERROR model_from_netlist: current output ' deblank(Iout(k,:))...
		' not in netlist model'])
      end
   end

end

if(~isempty(Vout))
   for k=1:size(Vout,1)
      str = ['Vnode' int2str(Vout(k,1)) '-Vnode' int2str(Vout(k,2))];
      outputs = strvcat(outputs,str);
   end

   CVt = zeros(ntheta,ngamma+ndelta); DVt = [zeros(ntheta,nphi) eye(ntheta)];
   dDVt = zeros(ntheta,nphi+ntheta);

   CVd = [zeros(ndelta,ngamma) eye(ndelta)]; DVd = zeros(ndelta,nphi+ntheta);
   dDVd = zeros(ndelta,nphi+ntheta);

   CVe = Peeinv*[Fge' -Fbe'*Rbbinv*Fbd]; DVe = Peeinv*[Fpe' -Fbe'*Rbbinv*Fbt];
   dDVe = zeros(neps,nphi+ntheta);

   temp = [Lzg+Lzz*Fgz' zeros(nzeta,ndelta)];
   CVz = [Rzz*Fgz' zeros(nzeta,ndelta)] - temp*Mhatinv*Rhat;
   DVz = [Rzz*Fpz' zeros(nzeta,ntheta)] + temp*Mhatinv*Vhat;
   dDVz = [Lzz*Fpz' zeros(nzeta,ntheta)] + temp*Mhatinv*What;

   CVa = -(Fat*CVt + Fad*CVd);
   DVa = -(Fat*DVt + Fad*DVd);
   dDVa = -(Fat*dDVt + Fad*dDVd);

   CVb = -(Fbt*CVt + Fbd*CVd + Fbe*CVe);
   DVb = -(Fbt*DVt + Fbd*DVd + Fbe*DVe);
   dDVb = -(Fbt*dDVt + Fbd*dDVd + Fbe*dDVe);

   CVg = -(Fgt*CVt + Fgd*CVd + Fge*CVe + Fgz*CVz);
   DVg = -(Fgt*DVt + Fgd*DVd + Fge*DVe + Fgz*DVz);
   dDVg = -(Fgt*dDVt + Fgd*dDVd + Fge*dDVe + Fgz*dDVz);

   CVp = -(Fpt*CVt + Fpd*CVd + Fpe*CVe + Fpz*CVz);
   DVp = -(Fpt*DVt + Fpd*DVd + Fpe*DVe + Fpz*DVz);
   dDVp = -(Fpt*dDVt + Fpd*dDVd + Fpe*dDVe + Fpz*dDVz);

% Build cost matrix for shortest path (dijkstra) algorithm. Add one to all
% node indices in case a node value of 0 is used.
   cost_matrix = Inf*ones(numNode+1);
   for k=1:length(netlist2.N1)
      cost_matrix(netlist2.N1(k)+1,netlist2.N2(k)+1)=1;
      cost_matrix(netlist2.N2(k)+1,netlist2.N1(k)+1)=1;
   end

   for k=1:size(Vout,1)
      idx1 = find(Vout_mod(k,1)==union(netlist2.N1,netlist2.N2));
      idx2 = find(Vout_mod(k,2)==union(netlist2.N1,netlist2.N2));
      if(isempty(idx1) | isempty(idx2))
         wait(['ERROR model_from_netlist: one of voltage output nodes ' ...
		int2str(Vout(orig_nodes(k),:)) ' is not in netlist model'])
      end

% Find the shortest path between the specified nodes.

      [path,cost]=dijkstra(numNode,cost_matrix,Vout_mod(k,1)+1,Vout_mod(k,2)+1);
      path = path-1;	% shift from indices back to node numbers 

% Next, identify branches that connect the nodes in path.

      clear branch vsign
      for m=1:length(path)-1
         for j=1:length(netlist2.N1)
            if(netlist2.N1(j)==path(m) & netlist2.N2(j)==path(m+1))
               branch(m) = j;
               vsign(m) = 1;
               break;
            elseif(netlist2.N2(j)==path(m) & netlist2.N1(j)==path(m+1))
               branch(m) = j;
               vsign(m) = -1;
               break;
            end
         end
         if(length(branch)<m)
            wait('ERROR model_from_netlist: programming error - find & fix');
         end
      end

% Finally, "sum the voltages" along the path, each with the appropriate sign.

      Csum = 0; Dsum = 0; dDsum = 0;
      for m=1:length(branch)
         Vout_index = branch(m);
         if(ismember(Vout_index, alpha))
            j = find(Vout_index==alpha);
            CVtemp = CVa(j,:); DVtemp = DVa(j,:); dDVtemp = dDsum+dDVa(j,:);
         elseif(ismember(Vout_index, beta))
            j = find(Vout_index==beta);
            CVtemp = CVb(j,:); DVtemp = DVb(j,:); dDVtemp = dDsum+dDVb(j,:);
         elseif(ismember(Vout_index, gamma))
            j = find(Vout_index==gamma);
            CVtemp = CVg(j,:); DVtemp = DVg(j,:); dDVtemp = dDsum+dDVg(j,:);
         elseif(ismember(Vout_index, phi))
            j = find(Vout_index==phi);
            CVtemp = CVp(j,:); DVtemp = DVp(j,:); dDVtemp = dDsum+dDVp(j,:);
         elseif(ismember(Vout_index, theta))
            j = find(Vout_index==theta);
            CVtemp = CVt(j,:); DVtemp = DVt(j,:); dDVtemp = dDsum+dDVt(j,:);
         elseif(ismember(Vout_index, delta))
            j = find(Vout_index==delta);
            CVtemp = CVd(j,:); DVtemp = DVd(j,:); dDVtemp = dDsum+dDVd(j,:);
         elseif(ismember(Vout_index, epsilon))
            j = find(Vout_index==epsilon);
            CVtemp = CVe(j,:); DVtemp = DVe(j,:); dDVtemp = dDVe(j,:);
         elseif(ismember(Vout_index, zeta))
            j = find(Vout_index==zeta);
            CVtemp = CVz(j,:); DVtemp = DVz(j,:); dDVtemp = dDVz(j,:);
         else
            CVtemp = 0; DVtemp = 0; dDVtemp = 0;
         end
         Csum = Csum+vsign(m)*CVtemp; 
         Dsum = Dsum+vsign(m)*DVtemp; 
         dDsum= dDsum+vsign(m)*dDVtemp;
      end
      Chat = [Chat; Csum]; Dhat = [Dhat; Dsum]; dDhat = [dDhat; dDsum];
   end
end

% Fix Vhat so that input channel numbers are as specified in the 
% original netlist that was input by user.

temp = Vhat; clear Vhat;
Vhat(:,V_original_values) = temp;

if(~isempty(Dhat))
   temp = Dhat; clear Dhat;
   Dhat(:,V_original_values) = temp;
end

temp = What; clear What;
What(:,V_original_values) = temp;

if(~isempty(dDhat))
   temp = dDhat; clear dDhat;
   dDhat(:,V_original_values) = temp;
end

temp = inputs; clear inputs;
inputs(V_original_values,:) = temp;
inputs = deblank(inputs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the output data structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

descriptions = struct( ...
'Mhat','circuit modified mutual inductance matrix', ...
'Rhat','circuit modified resistance matrix', ...
'Vhat','modified mapping from current/voltage sources to circuit voltages', ...
'What','mapping from current/voltage source derivatives to circuit voltages',...
'Chat','matrix multiplying states to produce current, voltage outputs',...
'Dhat','matrix multiplying inputs to produce current, voltage outputs', ...
'DWhat','matrix multiplying d(inputs)/dt to get current, voltage outputs', ...
'Pxx','mapping from states to full set of inductor currents/cap. voltages', ...
'rxx','matrix of all resistance in series with netlist M elements', ...
'Istates','string array defining current states', ...
'Istate_idx','indices in inductors (both M and L) of current states',...
'inputs','string array defining model inputs', ...
'outputs','string array defining model outputs', ...
'netlist','original netlist input, after sorting and combining Rs');

model = struct( ...
'Mhat',Mhat, ...
'Rhat',Rhat, ...
'Vhat',Vhat, ...
'What',What, ...
'Chat',Chat, ...
'Dhat',Dhat, ...
'DWhat',dDhat, ...
'Pxx',Pxx, ...
'rxx',diag(Rvec), ...
'Istates',Istates, ...
'Istate_idx',Istate_idx,...
'inputs',inputs, ...
'outputs',outputs, ...
'netlist',netlist2, ...
'descriptions',descriptions);

function [netlist,old_node_numbers] = make_nodes_contiguous(netlist_in)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  netlist_out = make_nodes_contiguous(netlist_in)
%
%  PURPOSE: Change node numbers so that the set of node numbers starts at 0 and
%		is contiguous. 
%
%  INPUT:
%	netlist_in = input_netlist
%
%  OUTPUT:
%	netlist_out = output netlist, equivalent to netlist_in except node numbering
%	old_node_numbers = sorted list of original node numbers before modification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

netlist = netlist_in;
all_nodes = union(netlist.N1,netlist.N2);
mm =  min(all_nodes);
old_node_numbers = all_nodes;
if(mm ~= 0)
   netlist.N1 = netlist.N1-mm;
   netlist.N2 = netlist.N2-mm;
end

all_nodes = union(netlist.N1,netlist.N2);
changed=1;
while(changed)
   changed=0;
   for k=1:length(all_nodes)
      if(all_nodes(k)~=k-1)
         idx1 = find(netlist.N1>k-1);
         idx2 = find(netlist.N2>k-1);
         netlist.N1(idx1) = netlist.N1(idx1)-1;
         netlist.N2(idx2) = netlist.N2(idx2)-1;
         changed=1;
         break;
      end
   end
   all_nodes =union(netlist.N1,netlist.N2);
end

%  save debug_netlist_model.mat

