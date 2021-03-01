%
% Generate a list of mfile dependencies
% in the form: 
% 
% target : dep1 dep2 dep3 depN
% dep1 :
% dep2 :
% dep3 :
% depN :
%
% INPUTS:
%   target - name of the makefile target
%   mcommand - mfile or function used to generate the target
%   depfilename - name of the file where the generated dependencies will be dumped.
%
% This follows the methodology outlined here:
% http://scottmcpeak.com/autodepend/autodepend.html
%
% In particular, here why you have the lines listing dependencies as targets with no prereqs:
% What is needed is a way to say that a particular prerequisite file, if missing, should be treated as
% changed (so the target will be rebuilt), but not cause an error. GNU make has an obscure feature that 
% does just this: if a file (1) appears as a target in a rule with no prerequisites and no commands, and
% (2) that file does not exist and cannot be remade, then make will rebuild anything which depends on that file and not report an error. 
% (See Chapter 4 of the make manual, "Rules without Commands or Prerequisites".)


function gendep(target,mcommand,depfilename)
    f = fopen(depfilename,'w');


    [fList,pList] = matlab.codetools.requiredFilesAndProducts(mcommand);
    
    targ_str = sprintf('%s : %s\n', target, strjoin(fList,' '));
    deps_str = sprintf('%s :\n', strjoin(fList, ' :\n'));
    fprintf(f, targ_str);
    fprintf(f, deps_str);

        



return
