  function point_list = read_point(point_file,rmemptyln);
%
% point_list=read_point(point_file); read point names into an array of characters
%
% Input:
%        point_file=	File name with point names in it
%        rmemptyln=     remove empty lines (optional) [1]
%
% Note: Each point name should start on a seperate line
%        line length no longer restricted
%       Comments can be included in file by starting the line with % or ?
%       For points you dont want to include in read simply put % at start of line
%
% Output:
%       point_list= list of point names found in file, one point per row

% read_point.m leuer 5-94

%  echo read_point on

  if nargin <2
     rmemptyln=1;
  end
  
  fid = fopen(point_file);
  if fid <= 0
    disp([' %Error: read_point: file not found: ',point_file])
  else
    count= 0;
    for i=1:10000
      line= fgetl(fid);
      if ~isstr(line), break, end
      line_leg= length(line);
%      line_leg= min(80,line_leg);
      fchar= 0;
      if line_leg~=0
%         line= line(1,1:line_leg);
         fchar= line(1,1);
      end
      if fchar =='%' | fchar=='?'
%        disp([' skiping over line ',int2str(i),' line: ',line]);
      else
         if line_leg~=0
            count = count + 1;
            id= findstr(line,' ');
	    if ~isempty(id)
	       line_leg= max(1,id(1)-1);
	    end   	       
	    point_list(count,1:line_leg)= line(1:line_leg);
	 elseif ~rmemptyln
            count = count + 1;
	 end
      end
    end
%    disp([' %OK read_point read: ',int2str(count),' pointnames from: ',...
%             point_file]);
    fclose(fid);
  end

  if ~isempty(point_list)
    point_list= deblank(point_list); %remove ending blanks
  end
% improved to read blank lines and read more than 10 characters jal 10-02
% truncate at line at first blank jal 3/03
