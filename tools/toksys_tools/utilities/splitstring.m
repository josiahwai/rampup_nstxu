function rv = splitstring( str, varargin )
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  len =splitstring( str, varargin )
%
%  PURPOSE: ARRAY = SPLITSTRING( STR, DELIM, ALLOWEMPTYENTRIES ) splits the
%		character string STR, using the delimiter DELIM (which must be a
%		character array). ARRAY is a cell array containing the resulting
%		strings. If DELIM is not specified, space delimiter is assumed (see
%		ISSPACE documentation). ALLOWEMPTYENTRIES should be a logical single
%		element, specifying weather empty elements should be included in the
%		results. If not specified, the value of ALLOWEMPTYENTRIES is false.
%
%  INPUT:
%	string = text string to find length of
%	varargin is either DELIM or ALLOWEMPTYENTRIES
%
%  OUTPUT:
%	len = length of character string
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  M.J. Lanctot on 02/28/2013 (taken from online matlab library)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delim = '';
AllowEmptyEntries = false;

if numel(varargin) == 2
        delim = varargin{1};
        AllowEmptyEntries = varargin{2};
elseif numel(varargin) == 1
        if islogical(varargin{1})
                AllowEmptyEntries = varargin{1};
        else
                delim = varargin{1};
        end
end

if isempty(delim)
        delim = ' ';
        ind = find( isspace( str ) );
else
        ind = strfind( str, delim );
end

startpos = [1, ind+length(delim)];
endpos = [ind-1, length(str)];

rv = cell( 1, length(startpos) );
for i=1:length(startpos)
        rv{i} = str(startpos(i):endpos(i));
end

if ~AllowEmptyEntries
        rv = rv( ~strcmp(rv,'') );
end
