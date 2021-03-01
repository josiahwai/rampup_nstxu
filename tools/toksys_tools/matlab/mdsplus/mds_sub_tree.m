function  mds_ot = mds_sub_tree(mds_in)
%
%  SYNTAX:
%          mds_ot= mds_sub_tree(mds_in)
%
%  PURPOSE:gets mds subtrees for next level down. 
%
%  INPUT: <default>
%    mds_in    = Input mds tree
%
%  OUTPUT:
%    mds_ot   = Output mds subtree names
%
%  WRITTEN BY:  Jim Leuer    ON      3/3/05
%
%  CHANGE LOG:  SMF 20140923 - Changed getnci call to use nid_numbers due to 
%                              fullpath not accepting wildcards.
%
% ==========================================================================
% @(#)mds_sub_tree.m	1.3 07/09/09 

  % Initialize return variable
  mds_ot = char([]);

  % For all mds paths provided as input
  for i = 1:size(mds_in,1)

     % Strip trailing blanks from the mds path 
     mds_inn = deblank(mds_in(i,:));

     % Call GETNCI() to get the NID_NUMBERS 
     mdscmd = ['getnci("\' mds_inn '.*","NID_NUMBER")'];
     [mds_nids,status] = mdsvalue(mdscmd); 
     if ~status
        mdscmd = ['getnci("\\\' mds_inn '.*","NID_NUMBER")'];
        [mds_nids,status] = mdsvalue(mdscmd); 
     end
  
     % Call GETNCI() to get the FULLPATH for all NIDs
     if status 
        mds_fpath = [];
        num_nids = length(mds_nids);
        mds_fpath = cell(1,num_nids);
        for j=1:num_nids
           mdscmd = ['getnci(' num2str(mds_nids(j)) ',"FULLPATH")'];
           mds_fpath{j} = char(mdsvalue(mdscmd));
        end

        % ?
        mds_char = char(mds_fpath);
        id = strmatch(upper(mds_inn), upper(mds_char));
        if ~isempty(id)
           mds_ot = strvcat(mds_ot,mds_char(id,:));
        end
     end

  end

  return

% testing

 mds_in = '\EFIT01::TOP';
 mds_ot = mds_sub_tree(mds_in)
 mds_in = mds_ot
 mds_ot = mds_sub_tree(mds_in)
 mds_in = mds_ot
 mds_ot = mds_sub_tree(mds_in)
 mds_in = mds_ot
 
