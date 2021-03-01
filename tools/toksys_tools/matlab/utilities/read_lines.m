   function [tit,narray,array]= read_lines(file_name)
%  function [tit,narray,array]= read_lines(file_name)
%  reads data from file file_name which is expected to be lines.dat
%  It reads sets of data until the eof is found or a n=-99 is incoundered
%
%  Format of file is:
%
%  Title for data set 1
%    n, data, data, .... where: n= number of records to follow in data set
%      x, y, z, ...        data to be read into array stored by row
%       : read n data records
%      x, y, z, ...        data to be read into array stored by row
%
%   repete above for next data set
%
%  output:
%    tit             = array containing titles read for each data set
%    narray(:,1)     = Number of rows in array for each data set
%    narray(:,2:end) = other data read in the 2nd data set record
%    array           = Data storage of all data sets stored sequential by row
%    length(narray(:,1)= Number of data sets

% Jim Leuer, General Atomics, 11-4-95
% jal fixed so uses 1st data as size of rows for remaining reads
% 1-15-98 narray=[] and array=[] initialize to prevent matlab5 caution
% 9-98 check for null or NULL in filename corrected

  error= -1;
%                                       open file
%  e= exist(file_name);
%  if (e~=2) file_name= 'lines.dat', end
  narray= [];
  array= []; 
  if strcmp(file_name,'null') | strcmp(file_name,'NULL') 
    return
  end
  e= exist(file_name);
  if (e~=2) 
    disp([' %ERROR read_lines: cant find file: ',file_name])
    tit= [' %ERROR read_lines: cant find file: ',file_name];
    narray= [];
    array=[];
    return
  end
%  disp([' Opening file: ',file_name]);
  fid= fopen(file_name);
  if (fid <= 0)
    error= -10;
    disp([' %ERROR read_lines: couldnt open file: ',file_name])
    tit= [' %ERROR read_lines: couldnt open file: ',file_name];
    narray= [];
    array=[];
    return
  end

  done= 0;
  icount= 0;
% ------------------------------------- begin loop on data sets in file
  while done== 0
    line1= fgetl(fid);
    if ~isstr(line1), done= 1; break, end
%    disp(['line1: ',line1]);
    line2= fgetl(fid);
    if ~isstr(line2), done= 2; break, end
%    disp(['line2: ',line2]);
    n= sscanf(line2,'%f');
    if (n <= 0) done= 2.5; break, end
    icount= icount+1;
    if(icount <= 1) 
       tit= line1;
    else
       tit= str2mat(tit,line1);
    end
    narray= [narray;n'];

%   read 1st line of data and use this length for all other records
    line3= fgetl(fid);
    i= 1;
    if ~isstr(line3), done= 2.8; break, end
    id= findstr(line3,','); % allow numbers to be seperated by commas
    if ~isempty(id) line3(id)= ' '; end
    num= sscanf(line3,'%f');
    num_len= length(num);
    array= [array;num'];
%      disp(['i: ',num2str(i),' line3: ',line3]);
%      disp(num');
%
    for i= 2:n(1)
      line3= fgetl(fid);
      if ~isstr(line3), done= 3; break, end
%      disp(['i: ',num2str(i),' line3: ',line3]);
      id= findstr(line3,','); % allow numbers to be seperated by commas
      if ~isempty(id) line3(id)= ' '; end
      num= sscanf(line3,'%f');
      num_l= length(num);
      if num_l < num_len, num=[num;zeros(num_len-num_l,1)]; end
      if num_l > num_len, num= num(1:num_len,1); end
%      disp(num');
      array= [array;num'];
    end
  end  

% ------------------------------------- end loop on data sets
  error= fclose(fid);
  return
