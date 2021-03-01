  function [sam,sdiff,s1,s2]= diff_mat(mat1,mat2,del,levels);

% DIFF_MAT compares two structures or two .mat files to see if the 
% variables are the same & have same values
%
% SYNTAX:
%   [sam,sdiff,s1,s2] = diff_mat(mat1,mat2,del (optional),levels (optional))
%   
%   Output arguments are optional. If not set the function will print
%   an output of the comparison to the screen.
%
%   Use the levels argument to set the number of levels of sub-structures
%   to compare (see description below for more info). 
%           
%
%                      diff_mat(mat1,mat2);   % Prints output of comparison
%         [sam,sdiff]= diff_mat(mat1,mat2,del);
%          sam=        diff_mat(mat1,mat2,0); % No approximate comparison
%   [sam,sdiff,s1,s2]= diff_mat(mat1,mat2);   % output variables in struct s1,s2
%
%                      diff_mat(mat1,mat2,del,levels); % Set the number of
%                                                     % levels of substructures
%                                                     % to compare 
%
% INPUT:
%    mat1,mat2= mat1 and mat2 are both .mat file names file names 
%                       to read and compare (example: matlab.mat)
%                       (assumes .mat extension if file doesnt exist)
%               The following types of data can be compared.  All others 
%               will be ignored.
%               'double','sparse','char','int8','uint8','int16',
%               'uint16', 'int32',
%               'uint32','function_handle', and 'struct'. 
%               NOTE: Whether or not a sub-structure is compared is
%               determined by the levels argument.
%
%    del=        Optional evaluation limit for approximate equal evaluation
%                default= eps; if 0 approximate evaluation skipped    
%    levels=     Optional number of levels of substructures to be recursively 
%                evaluated. 
%                   <= -1 (DEFAULT) No limit to the number of sub-structures
%                       to be evaluated.
%                   0.   No substructures will be evaluated.
%                   >=  1. Specified number of substructures will be evaluated. 
%
% OUTPUT:
%    sam  =        1 if mat1 and mat2 identical, 0 if not exactly identical
%    sdiff= structure of differences: (ex. sdiff.sam, sdiff.c1, sdiff.s1and2
%        .sam=     1 if mat1 and mat2 identical, 0 if not exactly identical
%        .c1=      variable names in .mat1 [char array]
%        .c2=      variable names in .mat2 [char array]
%
%        .s1not2=  list of variables in mat1 but not in mat2
%        .s2not1=  list of variables in mat2 but not in mat1
%        .s1and2=  list of variable in mat1 and mat2;
%         Below lists are subsets of .s1and2 which can be compared numerically
%         .s1eq2 =  list of variables that are identically 0 [norm()==0]
%         .s1neq2=  list of variables that are not identically 0 [norm()~=0]
%         .n1neq2=  norm of variables in list s1neq2
%
%         .s1apr2=  approximate variable list: norm(mat2-mat1)/norm(mat1)<=del]
%         .s1napr2= not approx. variable list: norm(mat2-mat1)/norm(mat1)>del]
%         .n1napr2= fractional norm of variables in list s1napr2
%
%         .s1nev2=  list of variables that cannot be evaluated (struc, cell...)
%
%     s1= structure containing all mat1 variables
%     s2= structure containing all mat2 variables
%
% NOTE: All matricies with different sizes are marked with "inf" norm's
%
% NOTE: use 'diff_mat' with 'struc_to_ws' to get 'neq' variables into env
%       [sam,sdiff,s1,s2]= diff_mat_load('good.mat','bad.mat');
%                          struct_to_ws(s1,[],[],[],sdiff.s1neq2,'_1') % ?_1
%                          struct_to_ws(s2,[],[],[],sdiff.s1neq2,'_2') % ?_2
%           
%
% RESTRICTIONS: Fields of type logical and cell will not be compared.

% =======================================
% Jim Leuer 7-9-2004 Leuer@fusion.gat.com
% =======================================
% Modified by Brian Sammuli on Dec. 8, 2008 to allow comparison of structures.

% default input


  if nargin <= 1
     disp('%Caution: diff_mat needs atleast 2 arguments')
     help diff_mat
     return
  elseif nargin <=2
     del= eps; %machine accuracy
     levels = -1; % Don't bother looking at any sub-structures;
  elseif nargin <=3
      levels = -1; % Don't bother looking at any sub-structures;
  end



  if exist(mat1)~=2 & exist([mat1 '.mat'])==2
        mat1= [mat1 '.mat'];
  end

  if exist(mat2)~=2 & exist([mat2 '.mat'])==2
        mat2= [mat2 '.mat'];
  end

  if exist(mat1)~=2
     disp(['% DIFF_MAT couldnt find file: ' mat1])
     return
  end
  if exist(mat2)~=2
     disp(['% DIFF_MAT couldnt find file: ' mat2])
     return
  end
% ========================
% start load
  s1= load(mat1,'-mat');
  s2= load(mat2,'-mat');


  [sam, sdiff] = diff_struct(s1,s2,del,levels);
  

% =======================  
% output
% =======================  

  if nargout<=0
    diff_disp(sdiff,del,mat1,mat2);    
  end    
  
  return
% =======================
  

% =======================
% testing
%  a1= [1,2,3];
%  b1= [1,2;3,4];
%  b2= [b1;b1];
%  b3(:,:,1)=b1;
%  b3(:,:,2)=2*b1;
%  c1= ['dog'];
%  c2= [c1;'cat'];
%  a1= [1,2,3];
%  save mat1 a1 b1 b2 b3 c1 c2
%  a1= -[1,2,3];
%  save mat2 a1 b1 b2 b3 c1 c2
%  save mat1 a1 b1
%  save mat2 b1 a1 
%  
%  diff_mat('mat1','mat2');
%  sam= diff_mat('mat1','mat2')
%  [sam,sd]= diff_mat('mat1','mat2')
%  
%  [sam,sd]= diff_mat(mat1,mat2)
%   
%  mat1='sim_112510_04100.mat';
%  mat2='sim_118526_04000.mat';
%  diff_mat(mat1,mat2);
%  [sam,sd]=diff_mat(mat1,mat2);

%  b1= [1,2;3,4];

% [sam,sd,s1,s2]= diff_mat1(mat1,mat2);

% Modifications:
% jal 4/28/06 '-mat'
