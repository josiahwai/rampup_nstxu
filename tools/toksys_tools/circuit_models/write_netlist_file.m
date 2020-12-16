function write_netlist_file(netlist,netlist_file)
 %
%  SYNTAX:  write_netlist_file(netlist,netlist_file)
%
%  PURPOSE:  Write out a formatted ascii netlist file from a netlist data
%	structure in matlab.
%
%  INPUT:
%    netlist = netlist data structure variable
%    netlist_file = name of netlist file to write
%
%  OUTPUT:
%    ascii file containing netlist
 
%  RESTRICTIONS:
%
%  METHOD:  
% 
%  WRITTEN BY:  Mike Walker 	ON 	9/13/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(netlist_file,'w');

for k = 1:length(netlist.N1)
   if(round(netlist.values(k))==netlist.values(k))
      fprintf(fid,'%s\t%d\t%d\t%d\n',netlist.names(k,:), ...
		netlist.N1(k),netlist.N2(k),round(netlist.values(k)));
   else
      fprintf(fid,'%s\t%d\t%d\t%f\n',netlist.names(k,:), ...
		netlist.N1(k),netlist.N2(k),netlist.values(k));
   end
end

fclose(fid);

