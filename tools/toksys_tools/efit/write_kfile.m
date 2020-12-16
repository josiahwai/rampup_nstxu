function write_kfile(IN1,INWANT,INS,efitin)
%  USAGE:  write_kfile(eq) % Creates k file from TokSys equilibrium
%  Alternatively, create k file with exact content of variables:
%  write_kfile(IN1,INWANT,INS,efitin)
%
%  PURPOSE: write k-files
%
%  INPUTS:  eq, TokSys equilibrium with fields
%  List shows: field name, alternative names, description
%  Note that some alternatives aren't exactly the same things
%shot, a shot number
%time, a time in seconds
%plasma, rog, cpasma, MEASURED total plasma current [A]
%expmp2, bp MEASURED signals from magnetic probes
%coils, fl, MEASURED signals from flux loops
%btor, bzero, toroidal field at rcentr, rzero
%rcentr, rzero, radius where btor, bzero is measured
%psibit, physics units per digitizer bit for flux loops
%bitmpi, physics units per digitizer bit for probes
%bitip, physics units per digitizer bit for rogowski
%bitfc, physics units per digitizer bit for F coil currents
%bitec, physics units per digitizer bit for E coil currents
%tgamma, tangent for measured pitch angles
%sgamma, uncertainties in tgamma
%fitdelz, boolean flag


%ISHOT, a shot number
%ITIME, a time in ms
%PLASMA, measured total plasma current [A]
%EXPMP2, measured signals from magnetic probes
%COILS, measured signals from flux loops
%BTOR, toroidal field at RCENTR
%RCENTR, radius where BTOR measured
%PSIBIT, physics units per digitizer bit for flux loops
%BITMPI, physics units per digitizer bit for probes
%BITIP, physics units per digitizer bit for rogowski
%BITFC, physics units per digitizer bit for F coil currents
%BITEC, physics units per digitizer bit for E coil currents
%tgamma, tangent for measured pitch angles
%sgamma, uncertainties in tgamma
%fitdelz, boolean flag

%MORE VARIABLES
%timeu, some other time, default 0
%qvfit, ?, default 0
%fwtsi, fitting weights for flux loops
%fwtcur, fitting weight for rogowski?, default 1
%limitr, ?, default -nw?
%fwtmp2, fitting weights for magnetic probes
%kffcur, number of knots in spline for ffprim
%kppcur, number of knots in spline for pprime
%fwtqa, ?, default 0
%ierchk, ?, default 1
%fwtbp, ?, default 0
%serror, ?, default 3e-2
%nextra, ?, default 2
%scrape, ?, default 4e-2
%itrace, ?, default 1
%xltype, ?, default 0
%rcentr, rzero, radius where btor, bzero is measured
%
%

%  OUTPUTS: kfiles to disk
%
%  RESTRICTIONS: 
%

%  WRITTEN BY:  Anders Welander ON 2016-05-31
%
%  MODIFICATION HISTORY:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 1
  struct_to_ws(IN1);

  if ~exist('ishot','var') & exist('shot','var')
    ishot = shot;
  end

  tims = round(1000000*time)/1000;
  ss = num2str(ishot); while length(ss)<6, ss = ['0' ss]; end
  tt = num2str(fix(tims)); while length(tt)<5, tt = ['0' tt]; end
  if mod(tims,1), tt = [tt '_' num2str(mod(tims,1))]; end
  filename = ['k' ss '.' tt];
  
  in1{ 1,1} = 'ishot';
  in1{ 1,2} = 'shot';
  in1{ 2,1} = 'itime';
  in1{ 2,2} = 'tims';
  in1{ 3,1} = 'plasma';
  in1{ 3,2} = 'rog';
  in1{ 3,3} = 'cpasma';
  in1{ 4,1} = 'expmp2';
  in1{ 4,2} = 'bp';
  in1{ 5,1} = 'coils';
  in1{ 5,2} = 'fl';
  in1{ 6,1} = 'btor';
  in1{ 6,2} = 'bzero';
  in1{ 7,1} = 'rcentr';
  in1{ 7,2} = 'rzero';
  in1{ 8,1} = 'psibit';
  in1{ 9,1} = 'bitmpi';
  in1{10,1} = 'bitip';
  in1{11,1} = 'bitfc';
  in1{12,1} = 'bitec';
  in1{13,1} = 'tgamma';
  in1{14,1} = 'sgamma';

  inwant{1,1} = 'fitdelz';
  
  
  fid = fopen(filename,'w');
  
  fprintf(fid,['&IN1 ' 10]);
  for i = 1:size(in1,1)
    S = upper(in1{i,1});
    for j = 1:size(in1,2)
      s = in1{i,j};
      if ~isempty(s) & exist(s,'var')
        dum = eval(s);
	for k = 1:numel(dum)
	  if k == 1
	    V = [S ' = '];
	  else
	    V = 32+zeros(1,length(S)+3);
	  end
	  if dum(k) > 0
	    x = ' ';
	  else
	    x = '';
	  end
	  m = num2str(dum(k));
	  if isa(dum,'logical')
	    if dum(k)
	      m = 'T';
	    else
	      m = 'F';
	    end
	  end
	  fprintf(fid,[V x m ',' 10]);
	end
	break % Don't look for alternatives
      end
    end
  end
  fprintf(fid,['/' 10]);
  
  fprintf(fid,['&INWANT ' 10]);
  for i = 1:size(inwant,1)
    S = upper(inwant{i,1});
    for j = 1:size(inwant,2)
      s = inwant{i,j};
      if ~isempty(s) & exist(s,'var')
        dum = eval(s);
	for k = 1:numel(dum)
	  if k == 1
	    V = [S ' = '];
	  else
	    V = 32+zeros(1,length(S)+3);
	  end
	  if dum(k) > 0
	    x = ' ';
	  else
	    x = '';
	  end
	  m = num2str(dum(k));
	  if isa(dum,'logical')
	    if dum(k)
	      m = 'T';
	    else
	      m = 'F';
	    end
	  end
	  fprintf(fid,[V x m ',' 10]);
	end
	break % Don't look for alternatives
      end
    end
  end
  fprintf(fid,['/' 10]);
    
  fclose(fid);
  
else
  % Create namelist files for variables 

end
