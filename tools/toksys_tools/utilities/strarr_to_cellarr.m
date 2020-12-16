function cellarr = strarr_to_cellarr(strarr)
 %
%  SYNTAX:  cellarr = strarr_to_cellarr(strarr)
%
%  PURPOSE:  Convert string array with multiple strings on each line, delimited
%	by spaces or tabs, to cell array.
%
%  INPUT:
%     strarr = either a string array or (an integer) file id corresponding to
%		a data file (open for reading) containing the string array.
%
%  OUTPUT:
%     cellarr = cell array containing strings defined by strarr
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	1/21/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)strarr_to_cellarr.m	1.2 01/21/10

if(ischar(strarr))

   for k=1:size(strarr,1)
      line = strarr(k,:);
      icnt=0; tok='';
      while(length(line)>0)
         icnt = icnt+1;
         [temp,line] = strtok(line);
         tok = strvcat(tok,temp);
      end
      for j=1:size(tok,1)
         cellarr{k,j} = remove_space(tok(j,:),0);
      end
   end

else

   fid = strarr;
   k=0;
   while(1)
      line = fgetl(fid);
      if(line==-1), break; end;
      k=k+1;
      icnt=0; tok='';
      while(length(line)>0)
         icnt = icnt+1;
         [temp,line] = strtok(line);
         tok = strvcat(tok,temp);
      end
      for j=1:size(tok,1)
         cellarr{k,j} = remove_space(tok(j,:),0);
      end

   end
end
