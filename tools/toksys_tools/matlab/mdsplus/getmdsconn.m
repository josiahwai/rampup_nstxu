function [connObj, ier] = getmdsconn(serverTok)
    %
    %  SYNTAX OPTIONS: 
    %      [connObj, ier] = getmdsconn(serverTok)
    %
    %  PURPOSE:  Establish a connection to an MDSplus server, return a
    %            connection object.
    %
    %  INPUT:
    %    serverTok  = Either a hostname or a tokamak shortname:
    %                   'DIII-D', 'D3D' -> 'atlas.gat.com'
    %       		    'JET'  -> 'mdsplus.jet.efda.org'
    %                   'NSTX' -> 'birch.pppl.gov:8501'
    %                   'KSTAR'-> '203.230.126.212:8005'
    %                   'EAST' -> '202.127.204.12:8000'
    %                   'EAST@GA' -> 'eastdata.gat.com'
    %
    %  OUTPUT:
    %    connObj	= connection object used for getting data 
    %    ier	= error code = 0 if OK, else > 0
    %
    % 
    %  EXAMPLES:
    %    [connObj, ier] = getmdsconn('d3d')
    %    [connObj, ier] = getmdsconn('NSTX')
    %    [connObj, ier] = getmdsconn('east@ga')
    %
    %
    %  RESTRICTIONS:  
    %  (1) 
    %
    %  METHOD: 
    %
    %  WRITTEN BY:  Jayson Barr 	ON 	2/14/18
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    import MDSplus.Connection

    hosts = {...% Tokamak       % Server                % Aliases
                'd3d',          'atlas.gat.com',        {'dIIId','d3d','dIII-d','d3-d'};...
                'd3dsim',       'iris.gat.com',         {'dIIIdsim','d3dsim','sim'};...
                'jet',          'mdsplus.jet.efda.org', {'jet'};...
                'nstxu',        'birch.pppl.gov:8501',  {'nstxu','nstx-u','nstx'};...
                'kstar',        '203.230.126.229:8005', {'kstar','k-star'};...
                'eastga',       'eastdata.gat.com',     {'eastga','east@ga','gaeast','eastdata'};...
                'east',         '202.127.204.12:8000',  {'east','eastasipp','east@asipp','asippeast'};...
            };

    [n_tokamaks,~] = size(hosts);
    %tokamak='';
    hostname='';
    for j=1:n_tokamaks
        %[test_tok, test_host, test_aliases] = hosts{j,:};
        [~, test_host, test_aliases] = hosts{j,:};
        [~, n_aliases] = size(test_aliases);
        
        for k=1:n_aliases
            if strcmpi(serverTok, test_aliases(k))
                %tokamak = test_tok;
                hostname = test_host;
                break
            end
            if ~isempty(hostname)
                break
            end
        end
        if ~isempty(hostname)
            break
        end
    end
    
    if isempty(hostname)
        %tokamak = serverTok;
        hostname = serverTok;
    end
    
    ier = 0;
    try
        connObj = Connection(hostname);
    catch err
        connObj = nan(0);
        ier=1;
    end
    
end
