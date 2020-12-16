function save_whos_data(whos_data)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX: save_whos_data(whos_data)
%
%  PURPOSE:  Save to file the information you get when typing "whos". 
%
%  INPUT:
%	whos_data = data structure you get by typing
%			>> whos_data = whos
%
%  OUTPUT: Local file called "whos_data.dat".
 
%  RESTRICTIONS:  
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	9/23/99
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = whos_data;
fid = fopen('whos_data.dat','w');
for k=1:size(whos_data,1)
  fprintf(fid,'%s\t%d\t%d\t%d\t%s\n',s(k).name,s(k).size,s(k).bytes,s(k).class);
end
fclose(fid);
