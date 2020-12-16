function [netlist,numNode,numV] = load_netlist(fname)
 %
%  SYNTAX:  [netlist,numNode,numV] = load_netlist(fname)
%
%  PURPOSE:  Takes a netlist (similar to SPICE), parses it to derive the
%		circuit equations, then solves them symbolically.  
%
%  INPUT:
%	fname = name of netlist file to load
%
%  OUTPUT:
%	netlist  = netlist information in a data structure
%	numNode  = total number of nodes in model
%	numV     = number of voltage sources
 
%  RESTRICTIONS:
%
%  METHOD: 
%
%  WRITTEN BY:  Mike Walker 	ON 	11/15/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @(#)load_netlist.m	1.3 07/28/10

names = '';
values = [];

[Name N1 N2 arg3]=textread(fname,'%s %s %s %s ','commentstyle','matlab');    

clear netlist
netlist.names = '';
%numNode=0;  %Number of nodes, not including ground (node 0).
numV = 0;
allnodes = [];
for k=1:length(Name)
   netlist.names = strvcat(netlist.names,Name{k});
   tt = N1(k); netlist.N1(k) = str2num(tt{:});
   tt = N2(k); netlist.N2(k) = str2num(tt{:});
   tt = arg3(k); tt=tt{:}; 
   if(tt(1:1)=='(')
      i1=findstr(',',tt); L = str2num(tt(2:i1-1));
      i2=findstr(')',tt); R = str2num(tt(i1+1:i2-1));
      netlist.values(k) = L + R*i;
   else
      netlist.values(k) = str2num(tt);
   end
   allnodes = union(union(allnodes,str2num(N1{k})),str2num(N2{k}));
%   numNode=max(str2num(N1{k}),max(str2num(N2{k}),numNode));
   if(strcmp(netlist.names(k,1:1),'V'))
      numV=numV+1;
   end
end
numNode = length(allnodes);
