  function set_save_v6;

% set_save_v6 sets default "save" to produce Version 6 compatable .mat files
%
% SYNTAX:
%         set_save_v6
%            

 
% Jim Leuer 22aug2007 Leuer@fusion.gat.com
% see
% http://www.mathworks.com/support/solutions/data/1-1NGF71.html?solution=1-1NGF71

  vr= ver;
  vr(1).Release;
  if ~strcmp(vr(1).Release,'(R14SP2)')
%    disp(['% CAUTION: set_save_v6 designed to run on Matlab7 Release 14 (R14SP2) not: ' ...
%         vr(1).Release]);
    return
  end

    javaobj = com.mathworks.services.Prefs;
    javaMethod('setBooleanPref',javaobj,'SaveMatfileAsUnicode',0);
    javaMethod('setBooleanPref',javaobj,'SaveMatfileCompressed',0);

  return

% testing

  set_save_v6
  dum= [pi pi*2];
  save dum dum

  load dum.mat
  whoss dum
