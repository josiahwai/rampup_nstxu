 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:
%         in_script
%
%  PURPOSE:  reads variable identified by string vnam into function from calling
%            area or from base area, in this order
%
%  INPUT: 
%     vnam=	name of variable to bring in from calling or base area
%               ex: vnam='vac_objects'; vnam='time'; vnam= 'ip'; ...
%
%  OUTPUT:
%
%    variable identified by string vnam is set in routine from calling or base area
%    dumm=     1 if variable set, 0 if variable not found and not set
%
%  EXAMPLE USE in Function to allow call of function as script
%
%  if nargin <= 0
%    vnam= 'rbs'; in_script; if ~dumm return; end
%    vnam= 'minimize_i'; in_script
%  end
%
%    NOTE: This script when used inside a function allows the function to act
%          like a script for input. If variable doesnt exist it does not set
%
%  RESTRICTIONS:
%       vnam must be a character variable
%
% Caution: overwrites variable identified in vnam in routine and variable dumm
%
% see also ot_script
 
%  METHOD: Good Template for converting a script to a dual function/script type 
%
%  WRITTEN BY:  Jim Leuer 	ON 	02sept2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if exist('vnam')~=1 or ischar(vnam)
       disp('%ERROR in_script needs an input string variable vnam');
    else
      dumm= 1;
      evalin('caller',[vnam ';'],'dumm=0;');
      if dumm
        str= [vnam '= evalin(''caller'',[vnam '';''])'];
	eval(str);
      else                           % look in base operating space
        dumm= 1;
        evalin('base',[vnam ';'],'dumm=0;');
        if dumm
          str= [vnam '= evalin(''base'',[vnam '';'']);'];
 	  eval(str);
        else % Not found in calling or base
          disp(['%NOTE: in_script: No variable ' vnam ' in caller or base'])
        end
      end
    end % exist('vnam')
