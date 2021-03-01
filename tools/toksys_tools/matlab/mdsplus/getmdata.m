function [ data_vals, time_vals, x_vals, ier, units, connObj, tree_found ] = getmdata( shotnum, ptname, treename, connOrTokamak,tokamak)
%  SYNTAX OPTIONS: 
%      [data_vals,time_vals,ier,units] = getmdata(shotnum, ptname, treename, connOrTokamak, getxdata)
%
%  PURPOSE:  Get mdsplus or ptdata data with reusable connections through mdsPlus library.
%
%  INPUT:
%    shotnum =      shot number
%    ptname	=       name of data - must be in single quotes
%    treename =     which tree to use for mdsplus data, or 'ptdata' (optional)
%                       Default Options depend on the server:
%                           Server:      Default treename:
%                           'D3D'     -> 'ptdata'
%                           'JET'     -> 'PPF'
%                           'NSTX'    -> 'wf'
%                           'KSTAR'   -> 'pcs_kstar'
%                           'EAST'    -> 'pcs_east'
%                           'EAST@GA' -> 'pcs_east'
%
%
%    connOrTokamak = tokamak name, server, or reusable mdsplus library connection object.
%                   (optional, default = 'd3d')
%                    -mdsplus library objects can be made with getmdsconn
%                    -Server address can be input as a string
%                    -Tokamak name will lookup default server addresses:
%                       Options:   (default = 'd3d')
%                           'D3D'     -> 'atlas.gat.com'
%                           'JET'     -> 'mdsplus.jet.efda.org'
%                           'NSTX'    -> 'birch.pppl.gov:8501'
%                           'KSTAR'   -> '203.230.126.212:8005'
%                           'EAST'    -> '202.127.204.12:8000'
%                           'EAST@GA' -> 'eastdata.gat.com'
%
%  OUTPUT:
%    data_vals	=   Data for ptname defined by input (complete data record)
%    time_vals =    Time values
%    ier	=       error code = 0 if OK, else > 0
%    units =        Units string for ptname if defined in the MDSplus DB.
%
% 
%  EXAMPLES
%    getmdata(176991, 'ip', 'ptdata', 'd3d')
%    getmdata(57035, 'MAGN/IPLA', 'PPF', 'jet')
%    getmdata(57035,'MAGN/IPLA','PPF','JET')
%    getmdata(109070,'ip','wf','NSTX')
%    getmdata(999991,'PCIP','PCS_KSTAR','kstar')
%    getmdata(73502,'pcrl01','pcs_east','east')
%    getmdata(73502,'pcrl01','pcs_east','east')
%
%  EXAMPLE re-using MDSplus connection:
%   [connObj, ier] = getmdsconn('d3d');
%   [ data_vals, time_vals, ier, units ] = getmdata(176991, 'ip', 'ptdata', connObj)
%   [ data_vals, time_vals, ier, units ] = getmdata(176991, 'ipmeas', 'efit01', connObj)
%   [ data_vals, time_vals, ier, units ] = getmdata(176991, 'li', 'efitrt1', connObj)
%
%   or:
%   [ data_vals, time_vals, ier, units, connObj ] = getmdata(176991, 'ip', 'ptdata', 'd3d')
%   [ data_vals, time_vals, ier, units ] = getmdata(176991, 'ipmeas', 'efit01', connObj)
%   [ data_vals, time_vals, ier, units ] = getmdata(176991, 'li', 'efitrt1', connObj)
%   
%
%
%  RESTRICTIONS:  
%  (1) Only works for 1-D and 2-D data right now
%
%  WRITTEN BY:  Jayson Barr ON 	8/16/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ier = 1;
    data_vals=nan(0);
    time_vals=nan(0);
    x_vals=nan(0);
    units = '';
    connObj=nan(0);
    tree_found=treename;
    
    
    % Parse input args:
    
    if nargin < 4 || isempty(connOrTokamak)==1
        tokamak='d3d';
        connOrTokamak = tokamak;
    end
    
    
    if nargin < 2 || ~ischar(ptname) || isempty(ptname)
        ier=1;
        return
    end
    
    % Make a new MDSplus connection if none was provided:
    if isjava(connOrTokamak) % It's a pre-existing connection object:
        if nargin<5
            tokamak='';
        end
        connObj = connOrTokamak;
    elseif ischar(connOrTokamak) % It's a tokamak name or server, connect:
        tokamak = connOrTokamak;
        [connObj, ier] = getmdsconn(connOrTokamak);
    end
    
   
    [~, ~, default_trees, ~, ~, ~] = get_defaults_for_tok(tokamak);
    
    
    % Check for and remove escape chars in the pointname:  
    %   Ex: if pointname includes a : in it, user must escape like this \:
    ptname(ptname=='\')=[];
    
    if isempty(treename)
        if ischar(tokamak) && ~isempty(tokamak)
            [~, ~, default_trees, ~, ~, ~] = get_defaults_for_tok(tokamak);
        else
            default_trees=[];
        end
        for j=1:length(default_trees)
            try_this_tree=default_trees{j};
            [ data_vals, time_vals, x_vals, ier, units, connObj, tree_found ] = getmdata( shotnum, ptname, try_this_tree, connObj, tokamak);
            if ier==0
                tree_found=try_this_tree;
                return
            end
        end
        
        [try_this_tree,ier] = ptname_find_mdstree(ptname, connObj);
        if ier==0 && ischar(try_this_tree) && ~isempty(try_this_tree)
            [ data_vals, time_vals, x_vals, ier, units, connObj, tree_found ] = getmdata( shotnum, ptname, try_this_tree, connObj, tokamak);
            if ier==0
                tree_found=try_this_tree;
                return
            end
        end
        
        return
    end
    
    % Special case for accessing ptdata:
    if strcmpi('ptd',treename)
        treename = 'ptdata';
    end
    
    if shotnum>900000 && strcmpi(treename,'ptdata')
        % Look at using mdsplus.Data instead of mdsplus.Connection object
        [data_vals,time_vals,ier] = getptd(shotnum,ptname,-inf,inf);
        x_vals=nan(0);
        units='';
        return
    end
    
    if strcmpi('ptdata',treename)
        if shotnum < 900000
            ptname = ['PTDATA("' ptname '",' num2str(shotnum) ')'];
        else
            ptname = ['PTDATA2("' ptname '",' num2str(shotnum) ')'];
        end
    elseif ~strcmpi(ptname(1),'\')
        ptname = ['\' ptname];
    end
    
    % Open the tree:
    if ~strcmpi('ptdata',treename)
        try
            connObj.openTree(treename, shotnum)
        catch err
            ier = 2;
            return;
        end
    end
    
    % Pull the pointname's data:
    try
        data_vals = connObj.get([ptname]);
        time_vals = connObj.get(['dim_of(' ptname ')']);
        
        data_vals = NATIVEvalue(data_vals);
        time_vals = NATIVEvalue(time_vals);
        
        try
            x_vals = connObj.get(['dim_of(' ptname ',1)']);
            x_vals = NATIVEvalue(x_vals);
        catch
            x_vals = NaN(0);
        end
    catch err
        data_vals=nan;
        time_vals=nan;
        x_vals=nan;
        ier=3;
        return
    end
    
    % Get the units:
    d_units = '';
    t_units = '';
    x_units=nan;
    if nargout > 4
        try
            d_units = connObj.get(['units(' ptname ')']);
            d_units = NATIVEvalue(d_units);
        catch err
            d_units = '';
        end
        try
            t_units = connObj.get(['units(dim_of(' ptname ',0))']);
            t_units = NATIVEvalue(t_units);
        catch err
            t_units = '';
        end
        try
            x_units = connObj.get(['units(' ptname ',1))']);
            x_units = NATIVEvalue(x_units);
        catch err
            x_units = nan;
        end
        if isnan(x_units)
            units=d_units;
        else
            units={d_units,x_units};
        end
    end
    
    % CHeck for errors:
    if isempty(data_vals) || isempty(time_vals)
        ier=3;
        return
    elseif length(data_vals) < 3
        ier=3;
        return
    elseif length(data_vals) ~= length(time_vals)
        ier=3;
        return
    elseif (time_vals(3)-time_vals(2))<=0
        ier=3;
        return
    end
    
    % Convert to doubles:
    data_vals = double(data_vals);
    time_vals = double(time_vals);
    x_vals = double(x_vals);
    
    % Check if time in s or ms: 
    if max(time_vals)>1000 || strcmpi(t_units,'ms') % Must be in ms:
        time_vals = time_vals .* 1e-3;
        t_units='s';
    end
    
    % Sanitize any weird points at the end of the array:
    test = time_vals(2:end) - time_vals(1:end-1);
    idx_zero = find(test<0,1);
    idx_nan = find(isnan(time_vals),1);
    
    if ~isempty(idx_zero)
        idx_zero = idx_zero + 1;
    else
        idx_zero = nan;
    end
    
    if isempty(idx_nan)
        idx_nan = nan;
    end
    
    if ~isnan(idx_zero) && ~isnan(idx_nan)
        idx_cut = min(idx_zero,idx_nan);
    elseif ~isnan(idx_zero)
        idx_cut = idx_zero;
    elseif ~isnan(idx_nan)
        idx_cut = idx_nan;
    else
        idx_cut = nan;
    end
    
    if ~isnan(idx_cut)
        time_vals = time_vals(1:(idx_cut-1));
        data_vals = data_vals(1:(idx_cut-1),1:end);
    end
    
    % Jump to relevant time period:
    j_start = find(-0.1 < time_vals, 1);
    time_vals = time_vals(j_start:end,1:end);
    data_vals = data_vals(j_start:end,1:end);
    
    ier=0;
    
    return
end