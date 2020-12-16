  function  [shot,tms,dr,server] = shot_from_gfile(gfile)
 %
%  SYNTAX:
%  [shot,tms,dr]         = shot_from_gfile(gfile) % pointer to file
%  [shot,tms,tree,server] = shot_from_gfile(gfile) % pointer to MDS+ tree
%
%  PURPOSE: look at string in gfile name and extract shot, tms, etc
%
%  INPUT: <default>
%    gfile   = string contining gfile name or "dot" pointer to MDS tree ex:
%              '/u/leuer/efit/d3d/shot131498/g131498.02000' % file pointer
%              'D3D.EFITRT1.g131498.02000'                  % MDS+ pointer
%
%  OUTPUT: <default>
%    shot  = shot number
%    tms   = shot time [ms]
%    dr    = directory if pointer to gfile
%    tree  = EFIT Tree if pointer to MDS tree <'EFIT01'>
%    server= MDS server for efit data <'D3D'>

% NOTE:
%           FULL Pointer to efit trees is: '\D3D::TOP.MHD.EFIT'
 
%  WRITTEN BY:  Jim Leuer    ON      17Mar2010 
% ==========================================================================

  if nargin==0
     disp('% shot_from_gfile needs at least a "gfile" argument')
     help shot_from_gfile
     return
  end

  id= strfind(gfile,'.'); % "dot" Format: SERVER.TREE.SHOTNUM.TIME	

  if isempty(id)
     disp(['%ERROR: shot_from_gfile: minimum form "shot.tms" ',gfile])
     shot=[]; tms=[]; dr=[]; server=[];

  else % isempty(id)

     [dr, nam, ext]= parse_filename(gfile);

% ---------------------
% 1) Time in ms
     id1= strfind(ext,'.'); 
     if length(ext)~= 6 | isempty(id1)
       disp([' %ERROR shot_from_gfile file extention needs format ".02600" ',ext])
     end
     tms= str2num(ext(2:end));
     if tms<0 | tms>99999
        disp([' %ERROR shot_from_gfile:TIME extention error ',ext]);
        tms=[];
     end

% ---------------------
% 2) Shot Number
     sshot= nam(end-5:end); % string of shotnum last 6 digits of nam
     shot= sscanf(sshot,'%d'); 
     if shot<0 | shot>999999
        disp([' %ERROR shot_from_gfile:SHOT number error ',sshot]);
        shot=[];
     end

     if length(id)<=1 % gfile points to file directory

% ---------------------
% 3a) Directory
       server=[];
       if exist(gfile) ~= 2
          disp(['%WARNING: shot_from_gfile: gfile doesnt exist ',gfile])
       end

     else % if length(id)<=1

% ---------------------
% 3b) tree = optional <EFIT01>
       if length(id)<=2 id0= 0; else id0= id(end-2); end
       idd= id0+1:id(end-1)-1;
       if ~isempty(idd)
         tree= gfile(idd); % this is TREE
       else
          tree= 'EFIT01';
          disp([' % NOTE shot_from_gfile:is setting TREE to ',tree]);
       end 
       dr= tree;  

% ---------------------
%  4)  server = optional <d3d>
       if length(id)<=2
         server= 'D3D';
       else
         if length(id)<=3 id0= 0; else id0= id(end-3); end
         server= gfile(id0+1:id(end-2)-1);
       end
% ----------------------
% Could check HERE to see if tree exists in MDS+ 
% Can do for d3d using get_mdsefit_nms
%      fln=[server '.' tree '.' int2str(shot) '.' stms];
      if strmatch(upper(server),'D3D') & ~isempty(shot) & ...
          exist('get_mdsefit_nms')==2
         [efit_nms] = get_mdsefit_nms(shot);
	 if isempty(strmatch(upper(tree),efit_nms))
	   disp([' % CAUTION shot_from_gfile cant find D3D tree in mds+ ' tree])
	   efit_nms
	 end
      end
     end % if length(id)<=1

  end % isempty(id)
 
  return
% =========================================================
% testing
  gfile = '/u/leuer/efit/d3d/shot131498/g131498.02600'
  gfile = '.g123216.03000'
  gfile = '.D3D.EFITRT1.131498.02600'
  gfile = '.D3D.dum.131498.02600'
  gfile = '131498.02600'
  gfile = '.xxx.EFITRT1.131498.02600'

     [dr, nam, ext]= parse_filename(gfile)

   [shot,tms,dr] = shot_from_gfile(gfile) % pointer to file
   [shot,tms,tree,server] = shot_from_gfile(gfile) % pointer to file

