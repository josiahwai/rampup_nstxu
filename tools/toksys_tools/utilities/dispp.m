  function count= dispp(array,fid,ndig);
%=============================================================================
% Prints numeric or character array like disp(array) but writes to file (fid).
%        Column heading and line spliting is not performed
%        Useful for making "diary" screen data go to file instead of to screen
%        also does multiple prints to screen and file or multiple files 
%
% SYNTAX:
%         dispp(array,fid);           % prints to fid (fid must be open file id)
%         count=dispp(array,fid);     % returns number of characters printed
%         dispp(array,fid,ndig);      % uses num2str(array,ndig) for numeric arr
%
% INPUT: [default]
%    array=   printed array (just like in disp(array))
%    fid=     open file id (use fid= fopen(file_name,'w'); [0] = screen
%             if fid is an array dispp writes to all files in fid list
%    ndig=    number of digits for numeric [default uses: num2str(array)]
%
% OUTPUT: 
%          prints data to file just like disp(array) does to screen
%          count=    number of array elements printed
%                 or number characters printed including \n line feeds
%
% Limitations: At present array must be all one type: number or char array
%              Column heading and line spliting is NOT performed like disp()
%=============================================================================

% =======================================
% Jim Leuer 05/04/2006 Leuer@fusion.gat.com
% =======================================
% default input

  if nargin <= 0
     help dispp
     count=[];
     return
  elseif nargin <= 1
     count= length(array(:));
     disp(array);
     return
  end
  
  if isempty(fid)
     fid=1; % print to screen
  end

% check for open files
  nfid= 1:length(fid);
  id=   ones(size(nfid));
  for ii=nfid
   [f1,p1,m1]= fopen(fid(ii)); % check for open file  
   if isempty(f1)
      id(ii)= 0;
      disp(['CAUTION dispp: file not opened, fid= ' int2str(fid(ii))])
   end
  end

  fid1= fid(find(id));  
  if isempty(fid1)
      count=[];
      disp(['ERROR dispp: Exiting; No files open, fid= ' int2str(fid)])
      return
   end
  
  if isnumeric(array)
    if nargin>=3
       array= num2str(array,ndig);
    else
       array= num2str(array);
    end     
  end
% ========================
  if ischar(array)
     len= size(array,2);
     f1= char('%'*ones(1,len));
     f2= char('c'*ones(1,len));
     f3= [f1;f2];
     f3= f3(:)';
     f3= [f3,'\n'];
     for ii=1:length(fid1)
         count= fprintf(fid1(ii),f3,array');
     end
  else
      disp(['ERROR dont understand array type ']);
      whos array;
      count=[];
  end
  
% ========================
 return

% ========================
% testing:

  array= ['abc';'def'];
  count= dispp(array,1);
  array= [1:5;6:10];
  count= dispp(array,1);
  count= dispp(array,1,8);
  array= [1:pi:50;51:pi:100];
  count= dispp(array,1);
  count= dispp(array,1,10);

  fidd= fopen('dum.dum','w');
  fid= [1, fidd, 999]; %
  count= dispp(array,fid,10);
  fclose(fidd);

