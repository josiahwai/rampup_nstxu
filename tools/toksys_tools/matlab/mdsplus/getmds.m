function [data,varargout] = getmds(shot,name,range_min,range_max,tree,server,ical);
 %
%  SYNTAX OPTIONS: 
%      [data,tvec,ier] = getmds(shot,name,range_min,range_max,tree,server,ical)
%      [data,xvec,tvec,ier] = getmds(shot,name,range_min,range_max,tree,server,ical)
%
%  PURPOSE:  Get mdsplus data from the MDS database. 
%
%  INPUT:
%    shot 	= shot number
%    name	= name of data - must be in single quotes
%    range_min	= minimum value of independent variable(s) - use NaN to specify no limit:
%                  if scalar, this is min time (s)
%                  if vector, last entry is min time (s), other entries are minima 
%                      for other independent variables
%    range_max	= maximum value of independent variable(s) - use NaN to specify no limit:
%                  if scalar, this is max time (s)
%                  if vector, last entry is max time (s), other entries are maxima
%                      for other independent variables
%    tree	= which tree to use (optional, default = 'EFIT01')
%    server     = which server to use (optional). Options: 'DIII-D'(default),
%		    'JET', 'NSTX', 'THOR'  
%    ical       = return data in units of digitizer counts (ical=0), volts into
%		digitizer (ical=2), or physics units (ical=1).  Only used for
%		data from 'THOR', 'data-server3' (EAST), and 'pcs_kstar' (KSTAR)
%
%  OUTPUT:
%    data	= data for name defined by input (dimensions correspond to dimensions
%                 specified in range_min)
%   (xvec..tvec) = variable number of arguments containing vectors of independent 
%                   variables specified in range_min and range_max 
%    ier	= error code = 0 if OK, else > 0
%
% 
%  EXAMPLES
%    getmds(106649,'rmaxis',1.5,5,'EFIT01','DIII-D')    
%    getmds(57035,'MAGN/IPLA',1.0,3.0,'PPF','JET')
%    getmds(109070,'ip',1.0,3.0,'wf','NSTX')
%    getmds(999991,'PCIP',0.0,8.0,'PCS_KSTAR','THOR');
%
%		NEED 3D examples
%
%  RESTRICTIONS:  
%  (1) Only works for range_min/max up to length 2 at the moment (make a request for more).
%  (2) Use the lower level routines starting with "mds" to manipulate the MDS
%      data directly. 

%  METHOD: 
%
%  WRITTEN BY:  Mike Walker 	ON 	6/4/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)getmds.m	1.12 08/08/13

% Error checking:

if (nargin < 4)
   disp('ERROR getmds: Must have at least 4 input arguments')
   return;
end;

ndim = length(range_min);
if(ndim~=length(range_max))
   disp('ERROR getmds: range_min and range_max must be equal length')
   return;
end

% test for existence of optional arguments

if (nargin < 5)
    tree = 'EFIT01';
end;
if (nargin < 6)
    server = 'DIII-D';
end;

if(ndim>1)
   xmin = range_min(1); tmin = range_min(2);
   xmax = range_max(1); tmax = range_max(2);
else
   tmin = range_min; tmax = range_max;
end

if(isempty(tmin))
   tmin=-1e+6;
end
if(isempty(tmax))
   tmax=1e+6;
end
if(nargin < 7)
   ical=1;
end
data=[];tvec=[];ier=0;
if(ical~=0 & ical~=1 & ical~=2)
   fprintf('ERROR getmds: invalid ical = %d\n',ical);
   ier=1; return;
end

upcase_server = upper(server);

% process default DIII-D server case
if ((nargin < 6) | ((nargin >= 6) & (strcmp(upcase_server,'DIII-D')))) 

    status = mdsconnect('atlas.gat.com');	% returns status=1 if success, else -1
    if status ~= 1
       disp(['ERROR getmds: unable to open server ' server])
       ier=1;
       data = [];
       varargout = getmds_varargout(nargout,ndim,ier,[],[]);
       return;
    end

    ier=0;
    [shoto,status] = mdsopen(tree,shot);  % returns status=1 on success, else 0
    if ~mod(status,2)
        ier=1;
        fprintf(['ERROR getmds: unable to open ' tree ' for shot '  ...
	int2str(shot) '\n'])
        varargout = getmds_varargout(nargout,ndim,ier,[],[]);
        return;
    end

    % get data vector:
    dname = ['\' name];
    [data,status] = mdsvalue(dname);
    if ~mod(status,2)
        ier=2;
        fprintf('ERROR getmds: unable to get data %s from %s\n',dname,tree)
        data = [];
        varargout = getmds_varargout(nargout,ndim,ier,[],[]);
        return;
    end

    % get time vector:
    if(ndim==1)
       tname = ['dim_of(\' name ')'];
    else
       tname = ['dim_of(\' name ',1)'];
    end
    [tvec,status] = mdsvalue(tname);
    if ~mod(status,2)
        ier=3;
        fprintf(['ERROR getmds: unable to get time data from ' tree '\n'])
        tvec = [];
        varargout = getmds_varargout(nargout,ndim,ier,[],[]);
        return;
    end
    if(min(tvec)>10 | max(tvec)>100)	% time units must be milliseconds
       tvec = tvec*1e-3;	% convert to seconds
    end

    % get x vector:
    if(ndim>1)
       xname = ['dim_of(\' name ',0)'];
       [xvec,status] = mdsvalue(xname);
       if ~mod(status,2)
           ier=3;
           fprintf(['ERROR getmds: unable to get x data from ' tree '\n'])
           xvec = [];
           varargout = getmds_varargout(nargout,ndim,ier,tvec,[]);
           return;
       end
    end

end;

% process non DIII-D server case
if ((nargin >= 6) & (~strcmp(upcase_server,'DIII-D')))

   if (strcmp(upcase_server,'JET'))
        status = mdsconnect('mdsplus.jet.efda.org');
        if status ~= 1
           disp(['ERROR getmds: unable to open server ' server])
           ier=1;
           data = [];
           varargout = getmds_varargout(nargout,ndim,ier,[],[]);
           return;
        end
        cmdString = '_s = JET("';
        cmdString = strcat('_s = JET("',tree,'/',name,'",',int2str(shot),')');
        data = mdsvalue(cmdString);
        tvec = mdsvalue('DIM_OF(_s)');
   elseif (strcmp(upcase_server,'NSTX'))
        status = mdsconnect('birch.pppl.gov:8501');
        if status ~= 1
           disp(['ERROR getmds: unable to open server ' server])
           ier=1;
           data = [];
           varargout = getmds_varargout(nargout,ndim,ier,[],[]);
           return;
        end
        
        ier=0;
        [shoto,status] = mdsopen(tree,shot);
        if ~mod(status,2)
           [shoto,status] = mdsopen(tree,shot);	% seems to work better to do it twice
        end
        if ~mod(status,2)
            ier=1;
            fprintf(['ERROR getmds: unable to open ' tree ' for shot '  ...
            int2str(shot) '\n'])
            data = [];
            varargout = getmds_varargout(nargout,ndim,ier,[],[]);
            return;
        end

        % get data vector:
        dname0=name;
        if(name(1:1)=='.')
           dname = name;	% relative reference
        elseif(strcmp(name(1:1),'\'))
           dname = name;
        else
           dname = ['\' name];
        end
        [data,status] = mdsvalue(dname);
        if ~mod(status,2)
            ier=2;
            fprintf('ERROR getmds: unable to get data %s from %s\n',...
                    deblank(dname),tree);
            data = [];
            varargout = getmds_varargout(nargout,ndim,ier,[],[]);
            return;
        end

        % get time vector:
        if(ndim==1)
           tname = ['dim_of(' dname ')'];
        else
           tname = ['dim_of(\' name ',1)'];
        end
        [tvec,status] = mdsvalue(tname);
        if ~mod(status,2)
            ier=3;
            fprintf(['ERROR getmds: unable to get time data from ' tree '\n'])
            tvec = [];
            varargout = getmds_varargout(nargout,ndim,ier,[],[]);
            return;
        end

        if(ndim>1) 		% get x vector
           xname = ['dim_of(\' name ',0)'];
           [xvec,status] = mdsvalue(xname);
           if ~mod(status,2)
               ier=3;
               fprintf(['ERROR getmds: unable to get x data from ' tree '\n'])
               xvec = [];
               varargout = getmds_varargout(nargout,ndim,ier,tvec,[]);
               return;
           end
        end

   elseif (strcmp(upcase_server,'THOR') | ...
	strcmp(upcase_server,'202.127.205.8') | ...
	strcmp(upcase_server,'203.230.125.212:8005') | ...
     (strcmp(upcase_server,'202.127.205.59') & strcmp(upper(tree),'PCS_EAST')))

      status = mdsconnect(upcase_server);
        if status ~= 1
           disp(['ERROR getmds: unable to open server ' server])
           ier=1;
           data = [];
           varargout = getmds_varargout(nargout,ndim,ier,[],[]);
           return;
        end

      ier=0;
      [shoto,status] = mdsopen(tree,shot);
      if ~mod(status,2)
         ier=1;
         fprintf(['ERROR getmds: unable to open ' tree ' for shot '  ...
      	 int2str(shot) '\n'])
           data = [];
           varargout = getmds_varargout(nargout,ndim,ier,[],[]);
         return;
      end

      if(strcmp(upper(tree),'PCS_EAST') | ...		% EAST PCS tree
	 strcmp(upper(tree),'PCS_KSTAR') | ...		% KSTAR PCS tree
		strcmp(upper(tree),'RDATA'))		% KSTAR non-PCS tree
         dname = ['\' name];
      else
         dname = name;
      end

    % get data vector:
      dname0 = upper(dname);
      if(ical==1)
         dname = upper(dname);
      elseif(ical==2)
         dname = [upper(dname) ':RAW'];
      elseif(ical==0)
         dname = [upper(dname) ':RAW'];
      end
      [data,status] = mdsvalue(dname);
      if ~mod(status,2)
         ier=2;
         fprintf(['ERROR getmds: unable to get data from ' tree '\n'])
         data = [];
         varargout = getmds_varargout(nargout,ndim,ier,[],[]);
         return;
      end

    % get time vector:
      if(ndim==1)
         tname = ['dim_of(' dname0 ')'];
      else
         tname = ['dim_of(' dname0 ',1)'];
      end
      [tvec,status] = mdsvalue(tname);
      if ~mod(status,2)
          ier=3;
          fprintf(['ERROR getmds: unable to get time data from ' tree '\n'])
          tvec = [];
          varargout = getmds_varargout(nargout,ndim,ier,[],[]);
          return;
      end

      if(ndim>1)			  % get x vector
         xname = ['dim_of(\' name ',0)'];
         [xvec,status] = mdsvalue(xname);
         if ~mod(status,2)
            ier=3;
            fprintf(['ERROR getmds: unable to get x data from ' tree '\n'])
            xvec = [];
            varargout = getmds_varargout(nargout,ndim,ier,tvec,[]);
            return;
         end
      end

      if(ical==2)
         attr_name = [upper(name) ':OFFSET'];
         [offset,status] = mdsvalue(attr_name);

         attr_name = [upper(name) ':FBITS'];
         [counts_to_volts,status] = mdsvalue(attr_name);

         data = (data+offset) * counts_to_volts;
         
         attr_name = [upper(name) ':SCALE'];
         [scale,status] = mdsvalue(attr_name);

      end


   else		% unknown server
        status = mdsconnect(server);
        if status ~= 1
           disp(['ERROR getmds: unable to open server ' server])
           ier=1;
           data = [];
           varargout = getmds_varargout(nargout,ndim,ier,[],[]);
           return;
        end
        
        ier=0;
        [shoto,status] = mdsopen(tree,shot);
        if ~mod(status,2)
            ier=1;
            fprintf(['ERROR getmds: unable to open ' tree ' for shot '  ...
            int2str(shot) '\n'])
            data = [];
            varargout = getmds_varargout(nargout,ndim,ier,[],[]);
            return;
        end

        % get data vector:
        dname = ['\' name];
        [data,status] = mdsvalue(dname);
        if ~mod(status,2)
            ier=2;
            fprintf(['ERROR getmds: unable to get data from ' tree '\n'])
            disp(data)
            data = [];
            varargout = getmds_varargout(nargout,ndim,ier,[],[]);
            return;
        end

        % get time vector:
        if(ndim==1)
           tname = ['dim_of(\' name ')'];
        else
           tname = ['dim_of(\' name ',1)'];
        end
        [tvec,status] = mdsvalue(tname);
        if ~mod(status,2)
            ier=3;
            fprintf(['ERROR getmds: unable to get time data from ' tree '\n'])
            tvec = [];
            varargout = getmds_varargout(nargout,ndim,ier,[],[]);
            return;
        end 

        if(ndim>1)			% get x vector
           xname = ['dim_of(\' name ',0)'];
           [xvec,status] = mdsvalue(xname);
           if ~mod(status,2)
               ier=3;
               fprintf(['ERROR getmds: unable to get x data from ' tree '\n'])
               xvec = [];
               varargout = getmds_varargout(nargout,ndim,ier,tvec,[]);
               return;
           end
        end


   end;
end

% return only the data requested:
if(isempty(tvec) | tvec(end) < tmin | tvec(1) > tmax)
     ier = 4;
     data=[]; xvec=[]; tvec=[];
     fprintf(['ERROR getmds: no data between tmin and tmax \n'])
elseif((ndim>1) & (isempty(xvec) | max(xvec) < xmin | min(xvec) > xmax))
    ier = 4;
    data=[]; xvec=[]; tvec=[];
    fprintf(['ERROR getmds: no data between xmin and xmax \n'])
elseif((ndim>1) & (isnan(xmin) || isnan(xmax))) % Don't truncate x-series
    [mm,m1]=min(abs(tvec-tmin));
    [mm,m2]=min(abs(tvec-tmax));
    tvec = tvec(m1:m2);
    data = data(:,m1:m2);
elseif (isnan(tmin) || isnan(tmax)) % Don't truncate t-series
    if(ndim>1)
       [nn,n1]=min(abs(xvec-xmin));
       [nn,n2]=min(abs(xvec-xmax));
       xvec = xvec(n1:n2);
       data = data(n1:n2,:);
    end
else
    [mm,m1]=min(abs(tvec-tmin));
    [mm,m2]=min(abs(tvec-tmax));
    tvec = tvec(m1:m2);
    if(ndim>1)
       [nn,n1]=min(abs(xvec-xmin));
       [nn,n2]=min(abs(xvec-xmax));
       xvec = xvec(n1:n2);
       data = data(n1:n2,m1:m2);
    else
       data = data(m1:m2);
    end
end

if(ndim>1)
   varargout = getmds_varargout(nargout,ndim,ier,tvec,xvec);
else
   varargout = getmds_varargout(nargout,ndim,ier,tvec);
end

mdsdisconnect;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = getmds_varargout(nargs,ndim,ier,tvec,xvec)

out = cell(1);
if(ndim>1)
   out{1} = xvec;
   out{2} = tvec;
else
   out{1} = tvec;
end
if(nargs > ndim+1)
   out{nargs-1} = ier;
end
