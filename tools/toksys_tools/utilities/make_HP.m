function make_HP(varargin)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX: make_HP
%
%  PURPOSE:  Customizes the specified template makefile. Customized makefile 
%  can be used to build an executable from the generated model code. Used
%  only by SIMULINK realtime code generation.
%
%  Modified from make_rtw.m provided with Matlab. 
%
%  INPUT:
%	varargin
%
%  OUTPUT:
%
%  RESTRICTIONS:
%
%  METHOD:  
%       MAKE_RTW first invokes the Target Language Compiler to generate the
%       code and then invokes the language specific make procedure.
%
%  WRITTEN BY:  Mike Walker 	ON 	11/14/97
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------%
% Get model name and build arguments %
%------------------------------------%

if nargin > 0 & length(varargin{1}) > 4 & varargin{1}(1:4)=='mdl:'
  modelName = varargin{1}(5:end);
  firstBuildArg = 2;
else
  modelName = bdroot;
  firstBuildArg = 1;
end
if isempty(modelName) 
  error('Unable to obtain current model');
end
hModel = get_param(modelName,'handle');

buildArgs = '';
for i=firstBuildArg:nargin
  buildArgs = [buildArgs, varargin{i}, ' '];
end
if length(buildArgs) > 0
  buildArgs(end) = [];
end

deleteRTWFile = strcmp(get_param(hModel,'RTWRetainRTWFile'),'off');

fprintf(1,'\n### Starting RTW build procedure for model: %s\n',modelName);

%--------------------%
% Platform specifics %
%--------------------%

%
% Use development sandbox, if TMW_V5_SANDBOX environmental variable is set, and
% the directories rtw, simulink/include and extern/include in the sandbox 
% exist. Otherwise use matlabroot.
%
tmwV5Sandbox = getenv('TMW_V5_SANDBOX');
if ~isempty(tmwV5Sandbox) & ...
   exist(fullfile(tmwV5Sandbox,'rtw'))==7 & ...
   exist(fullfile(tmwV5Sandbox,'simulink', 'include'))==7 & ...
   exist(fullfile(tmwV5Sandbox,'extern', 'include'))==7 ...
  rtwroot = fullfile(tmwV5Sandbox,'rtw');
  disp(['### Using rtwroot = ',rtwroot]);
else
  rtwroot = fullfile(matlabroot, 'rtw');
end

%---------------------------------------%
% Get target type from specified solver %
%---------------------------------------%

switch (get_param(hModel,'Solver')) 
  case {'FixedStepDiscrete', 'ode1', 'ode2', 'ode3', 'ode4', 'ode5'}
    targetType    = 'RT';
  case {'VariableStepDiscrete', 'ode45', 'ode23', 'ode113', 'ode15s', 'ode23s'}
    targetType    = 'NRT';
    error(['No support for variable step solvers in this release']);
  otherwise
    error(['Unhandled solver: ',get_param(hModel,'Solver'), ...
           ' in make_rtw.m']);
end


%-----------------------------------------%
% When needed, add Language to dialog box %
%-----------------------------------------%

language    = 'C';
languageDir = 'c';


%----------------------------%
% Get the system target file %
%----------------------------%
targetDir = [];  % Don't know it yet

systemTargetFile = deblank(get_param(hModel,'RTWSystemTargetFile'));
systemTargetFile = fliplr(deblank(fliplr(systemTargetFile)));

if isempty(systemTargetFile)
  error('No system target file specified');
end

k = find(isspace(systemTargetFile)==1);
if ~isempty(k) 
  tlcArgs          = systemTargetFile(k(1):end);
  systemTargetFile = systemTargetFile(1:k(1)-1);
else
  tlcArgs = [];
end

%---------------------------------------%
% Get RTWVerbose setting (default is 1) %
%---------------------------------------%
if length(tlcArgs) >= length('-aRTWVerbose=0') & ...
      ~isempty(findstr('-aRTWVerbose=0',tlcArgs))
  RTWVerbose = 0;
else
  RTWVerbose = 1;  % default
end

fid = fopen(systemTargetFile,'r');
if fid == -1
  %
  % Load the system target file from the default location:
  % 
  %   <rtwroot>/<languageDir>/<targetDir>/<systemTargetFile>
  %
  [file,targetDir]=GetSysTargetFileOrTMF(rtwroot,languageDir,targetDir, ...
					 systemTargetFile);
  if ~isempty(file)
    fid = fopen(file,'r');
    if fid ~= -1
      systemTargetFile = file;
    end
  end
else
  systemTargetFile = which(systemTargetFile);
end

if fid == -1
  error(['Unable to locate system target file: ',systemTargetFile]);
end


%--------------------------------------------------------------------------%
% Get the TargetType (RT or NRT) and Language from the system target file  %
%     %assign TargetType  = "<RT | NRT>" 	  		           %
%     %assign Language    = "<lang>"                                       %
%--------------------------------------------------------------------------%

tlcTargetType = [];  % unknown so far
tlcLanguage   = [];  % unknown so far

while (1)
  line = fgetl(fid);
  if ~isstr(line), break; end

  %
  % Use sscanf to throw away white space, i.e. match "token1 token2 = value"
  %
  [result,count] = sscanf(line, '%s%s%1s%s%1s');

  if count == 4 & length(result) > 19 & ...
     result(1:19) == '%assignTargetType="' & result(end) == '"'

    % Found:    %assign TargetType = "value"

    if ~isempty(tlcTargetType)
      fclose(fid);
      error(['Duplicate TargetType variable defined in ', systemTargetFile]);
    end
    result([1:19, end]) = [];
    tlcTargetType = result;
    if ~strcmp(tlcTargetType,'RT') & ~strcmp(tlcTargetType,'NRT')
      fclose(fid);
      error(['TargetType defined in ',systemTargetFile,' must be RT or NRT']);
    end

  elseif count == 4 & length(result) > 17 & ...
         result(1:17) == '%assignLanguage="' & result(end) == '"'

    % Found:    %assign Language = "value"

    if ~isempty(tlcLanguage)
      fclose(fid);
      error(['Duplicate Language variable defined in ', systemTargetFile]);
    end

    result([1:17, end]) = [];
    tlcLanguage = result;

  end

end

fclose(fid);

if isempty(tlcTargetType)
  error(['TargetType not specified in ', systemTargetFile]);
end

if ~strcmp(tlcTargetType,targetType)
  error(['The system target file is configured for ',tlcTargetType, ...
         ' whereas ', modelName, ' is configured for ',targetType]);
end


if isempty(tlcLanguage)
  error(['Language not specified in ', systemTargetFile]);
end

if ~strcmp(tlcLanguage,language)
  error(['The system target file is configured for language ',tlcLanguage, ...
         ' whereas ', modelName, ' is congigured for ', language]);
end


%-------------------%
% Setup TLC options %
%-------------------%

tlcOpts.targetType       = targetType;
tlcOpts.inlineParameters = get_param(hModel,'RTWInlineParameters');

%
% S-functions: 
%
comp = computer;
ispc = strcmp(comp(1:2),'PC');

sfcnBlks = find_system(hModel,'LookUnderMasks','on','FollowLinks','on', ...
                       'BlockType','S-Function');

sfcns           = {};
sfcnsInc        = '';
noninlinedSFcns = {};
add_fft         = 0;  % For DSP Blockset

for i = 1:length(sfcnBlks)
  sfcnBlk  = sfcnBlks(i);
  if strcmp(get_param(sfcnBlk, 'Tag'), 'Stateflow S-Function')
     stateflowSfcn = [modelName, '_rtw'];
     stateflowFlag = 1;
  else
     stateflowFlag = 0;
  end
  sfcn     = deblank(get_param(sfcnBlk,'FunctionName'));
  sfcn     = fliplr(deblank(fliplr(sfcn)));

  if ispc, sfcn = lower(sfcn); end;

  foundDup = 0;
  for j = 1:length(sfcns)
    if strcmp(sfcn,sfcns(j))
      foundDup = 1;
      break;
    end
  end

  if ~foundDup
    sfcns{end+1} = sfcn;

    sfcnLoc = which(sfcn);          % Get S-function location and extension
    if isempty(sfcnLoc)
      sfcnExt = mexext; % assume mex extension
    else
      if exist(sfcn) == 3, sfcnExt = mexext;, else, sfcnExt = 'm'; end
      idx = findstr(sfcnLoc,[sfcn,'.',sfcnExt]);
      sfcnLoc(idx(end):end) = [];
    end

    % See if we need to inline the S-function

    sfcn_tlc_file = [sfcn, '.tlc'];
    inlineSFcn    = ~isempty(dir(sfcn_tlc_file));

    if ~inlineSFcn & ~isempty(sfcnLoc) & ... 
       ~isempty(dir(fullfile(sfcnLoc, sfcn_tlc_file)))
      inlineSFcn = 1;
      sfcnsInc = [sfcnsInc, ' -I', sfcnLoc];
    end

    if ~inlineSFcn & strcmp(sfcnExt,'m')
      % Don't add m-file S-function to the list of items to be compiled
      inlineSFcn = 1;
    end

    if ~isempty(sfcnLoc)                  % Look in <sfcnLoc>/tlc_<lang> 
      sfcn_tlc_dir = fullfile(sfcnLoc, ['tlc_', languageDir]);
      if exist(sfcn_tlc_dir) == 7
        sfcnsInc = [sfcnsInc, ' -I', sfcn_tlc_dir];
        if ~isempty(dir(fullfile(sfcn_tlc_dir, sfcn_tlc_file)))
          inlineSFcn = 1;
        end
      end
    end

    if ~inlineSFcn
      if stateflowFlag == 0 
        noninlinedSFcns{end+1} = sfcn;
        if strcmp(sfcn,'sfft') | strcmp(sfcn,'srifft') | strcmp(sfcn,'srfft') | strcmp(sfcn,'srfftc')
    
          add_fft = 1;
        end
      else
        noninlinedSFcns{end+1} = stateflowSfcn;
      end
    end
  end
end

if add_fft
  noninlinedSFcns{end+1} = 'fft';
end

tlcOpts.sfcnsInc = sfcnsInc;

%------------------------------------------%
% Build the Stateflow target, if necessary %
%------------------------------------------%
 
BuildStateflowTarget( hModel );

%-----------------------------------------------------------------------%
% Invoke the Target Language Compiler to generate the specific language %
%-----------------------------------------------------------------------%

eval(['[rtwFile,modules]=tlc_',languageDir, ...
      '(modelName,rtwroot,systemTargetFile,tlcOpts,tlcArgs);']);

if strcmp(get_param(hModel,'RTWGenerateCodeOnly'),'on')
  if deleteRTWFile, delete(rtwFile); end;
  disp(['### Successful completion of RTW build procedure for model: ', ...
        modelName]);
  return;
end


%---------------------------%
% Get the template makefile %
%---------------------------%

templateMakefile = deblank(get_param(hModel,'RTWTemplateMakefile'));
templateMakefile = fliplr(deblank(fliplr(templateMakefile)));
if isempty(templateMakefile)
  error('No template makefile specified');
end

if exist(templateMakefile) ~= 2 
  [file,targetDir]=GetSysTargetFileOrTMF(rtwroot,languageDir,targetDir, ...
			                 templateMakefile);
  if ~isempty(file)
    templateMakefile = file;
  else
    error(['Unable to locate template makefile: ',templateMakefile]);    
  end
end

%-------------------------------------------------------------------------%
% Get the following from the model.rtw file:                              %
%      solver    - can't use get_param because we can default to discrete %
%      fixedStep - fixed step size                                        %
%      ncstates  - number of continuous states                            %
%      numst     - number of sample times                                 %
%      tid01eq   - '0' or '1'                                             %
%-------------------------------------------------------------------------%

fid = fopen(rtwFile,'r');
if fid == -1,
  error(['Unable to open ',rtwFile]);
end

solver    = [];
ncstates  = [];
numst     = [];
tid01eq   = [];
stepsize  = [];

while (1)
  line = fgetl(fid); if ~isstr(line), break; end

  if isempty(solver)
    if length(line) > 8 & line(1:8) == '  Solver'
       line(1:8) = [];
       solver = sscanf(line,'%s');
    end
  else 
    [parsedLine,count] = sscanf(line,'%s%g%1s');
    if count == 2
      parsedLine = sscanf(line,'%s%s%1s');
      if isempty(stepsize)
	if length(parsedLine) > 9 & parsedLine(1:9) == 'FixedStep' ... 
                                        & parsedLine(10:10) ~= 'O'
	  parsedLine(1:9) = [];
	  stepsize = parsedLine;
        end
      elseif isempty(tid01eq)
        if length(parsedLine) > 7 & parsedLine(1:7) == 'TID01EQ'
          parsedLine(1:7) = [];
          tid01eq = parsedLine;
        end
      elseif isempty(ncstates)
	if length(parsedLine) > 13 & parsedLine(1:13) == 'NumContStates'
	  parsedLine(1:13) = [];
	  ncstates = parsedLine;
	end
      elseif length(parsedLine) > 14 & parsedLine(1:14) == 'NumSampleTimes'
	parsedLine(1:14) = [];
	numst = parsedLine;
	break;
      end
    end
  end
end

fclose(fid);
if deleteRTWFile, delete(rtwFile); end;

if isempty(numst)
  error(['Unexpected contents in ',rtwFile]);
end

%----------------------------------------------%
% Invoke the language specific build procedure %
%----------------------------------------------%

k = findstr(systemTargetFile,filesep);
if ~isempty(k)
  buildOpts.sysTargetFile = systemTargetFile(k(end)+1:end);
else 
  buildOpts.sysTargetFile = systemTargetFile;
end
buildOpts.noninlinedSFcns = noninlinedSFcns;
buildOpts.solver          = solver;
buildOpts.tid01eq         = tid01eq;
buildOpts.ncstates        = ncstates;
buildOpts.numst           = numst;
buildOpts.modules         = modules;
buildOpts.stepsize        = stepsize;

eval(['rtw_',languageDir, ...
      '(modelName,rtwroot,templateMakefile,buildOpts,buildArgs)']);
disp(['### Successful completion of RTW build procedure for model: ', ...
      modelName]);

%end make_rtw


%-----------------%
% Local functions %
%-----------------%

% Function: GetSysTargetFileOrTMF =============================================
% Abstract:
%	Look in <rtwroot>/<langDir>/<targetDir> for file.
%
function [fileOut,targetDirOut] = GetSysTargetFileOrTMF(rtwroot,langDir, ...
                                                        targetDir,fileIn)

  fileOut      = [];
  targetDirOut = targetDir;

  if ~isempty(targetDir)

    file = fullfile(rtwroot, langDir, targetDir, fileIn);
    if exist(file) == 2
      fileOut = file;
    end

  else

    targetDirs = dir(fullfile(rtwroot, langDir));
    for i=1:length(targetDirs)
      if targetDirs(i).isdir
        targetDir = targetDirs(i).name;
        if ~strcmp(targetDir,'.') & ...
           ~strcmp(targetDir,'..') & ...
           ~strcmp(targetDir,'src') & ...
           ~strcmp(targetDir,'libsrc') & ...
           ~strcmp(targetDir,'lib') & ...
           ~strcmp(targetDir,'tlc') 
          file = fullfile(rtwroot, langDir, targetDir, fileIn);
          if exist(file) == 2
            fileOut      = file;
            targetDirOut = targetDir;
            break;
          end
        end
      end
    end
  end

%end GetSysTargetFileOrTMF



% Function: BuildStateflowTarget ==============================================
% Abstract:
%      Check for a Stateflow target named "rtw".  If nonexistent, create a
%      default target from the "sfun" target.
%
function BuildStateflowTarget( hModel )
 
  stateflowBlks = find_system(hModel,'LookUnderMasks','on', ...
			      'MaskType','Stateflow');
  if ~isempty(stateflowBlks)
     % Detected a Stateflow block => let's build RTW Stateflow target
     if ~sf('License','coder')
       error(['To build RTW with Stateflow blocks requires a valid ', ...
	      'Stateflow Coder license.']);
     end
     chartId = get_param(stateflowBlks(1),'UserData');
     machineId = sf('get',chartId,'.machine');
     targets = sf('TargetsOf',machineId);
     % See if there is RTW Stateflow target for this machine
     rtwTarget = sf('find',targets,'.name','rtw');
     sfunTarget = sf('find',targets,'.name','sfun');
     switch length(rtwTarget)
     case 0 % No RTW targets => create one
	disp('### Creating default Stateflow RTW target.');
	rtwTarget = sf('new','target'...
	   ,'.name','rtw'...
	   ,'.description','Default RTW target.'...
	   ,'.customCode',sf('get',sfunTarget,'.customCode')...
	   ,'.codeCommand','sfc -rtw'...
	   ,'.codeFlags',''...
	   ,'.makeCommand',''...
	   ,'.linkNode.parent',machineId...
	   ,'.document',''...
	);
     case 1
     otherwise
	warning('Multiple Stateflow RTW targets. Using the first one.');
	rtwTarget = rtwTarget(1);
     end
     status = sf('Private','targetman','rtwcode',rtwTarget,1);
     if status~=0
	error('Problem building Stateflow RTW target.');
     end
  end
 
%end BuildStateflowTarget

%[eof] make_rtw.m
