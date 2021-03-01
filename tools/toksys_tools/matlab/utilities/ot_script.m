 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:
%         ot_script
%
%  PURPOSE:  writes variable identified by string vnam into base area
%
%  INPUT: 
%     vnam=	name of variable to write to the base area
%               ex: vnam='vac_objects'; vnam='time'; vnam= 'ip'; ...
%
%  OUTPUT:
%
%    variable identified by string vnam is set in base area
%
%  EXAMPLE USE in Function to allow call of function as script
%
%   if nargout<=0
%     vnam= 'rbs'; ot_script; 
%     vnam= 'minimize_i'; ot_script
%  end
%
%    NOTE: This script when used inside a function allows the function to act
%          like a script for output to base. If variable doesnt exist it isnt set
%
%  RESTRICTIONS:
%       vnam must be a character variable
%       presently only puts output variables in base (not calling routine)
%
% Caution: overwrites variable identified in vnam in base
%
% See also in_script
 
%  METHOD: Good Template for converting a script to a dual function/script type 
%
%  WRITTEN BY:  Jim Leuer 	ON 	02sept2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if exist('vnam')~=1 or ischar(vnam)
       disp('%ERROR _script needs an input string variable vnam');
    else
%      str= ['assignin(''caller'',' '''' vnam '''' ',' vnam ');']
      str= ['assignin(''base'',' '''' vnam '''' ',' vnam ');']
      eval(str);
    end % exist('vnam')

% todo: figure out how to make work with calling function rather than base.
