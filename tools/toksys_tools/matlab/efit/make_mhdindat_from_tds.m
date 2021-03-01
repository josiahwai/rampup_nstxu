function make_mhdindat_from_tds(tok_data_struct,options)
%  SYNTAX:  make_mhdindat_from_tds(tok_data_struct,options)
%
%  PURPOSE:  Creates EFUND input file mhdin.dat from Toksys tok_data_struct
%
%  INPUT:
%     tok_data_struct
%     options = structure
%        fcid = Connections between fcoils (see cccirc in build_tokamak_system) [1:nc]
%        use_vv = 1:Add VV elements as coils 0:Don't [0]
%        vvid = VV element conenctions (see vvid in build_tokamamk_system) [1:nv]
%
%  OUTPUT: mhdin.dat file
% 
%  RESTRICTIONS:  

%  WRITTEN BY:  Nick Eidietis 	ON 	2016-05-18
%                 based upon Ander Welanders make_file_mhdin_dat4east
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   if exist('options','var')
      struct_to_ws(options);
   end
   struct_to_ws(tok_data_struct);
   if ~exist('use_vv','var')
      use_vv = 0;
   end


   bpnam=' ';
   for j=1:size(bpnames,1)
      bpnam = [bpnam '''' deblank(bpnames(j,:)) ''','];
      if j/5==round(j/5), bpnam=[bpnam 10 32]; end
   end
   xmp2='  ';
   for j=1:size(bpnames,1)
      s = num2str(bpdata(2,j)); while length(s)<5, s=[s '0']; end
      xmp2 = [xmp2  s ', '];
      if j/6==round(j/6), xmp2=[xmp2 10 32 32]; end
   end
   ymp2=' ';
   for j=1:size(bpnames,1)
      s = num2str(bpdata(1,j));
      if bpdata(1,j)==0, s='0.';end
      if bpdata(1,j)>=0, s=[' ' s];end    
      while length(s)<6, s=[s '0']; end
      ymp2 = [ymp2  s ','];
      if j/6==round(j/6), ymp2=[ymp2 10 32]; end
   end
   amp2='';
   for j=1:size(bpnames,1)
      s = [num2str(bpdata(3,j)) ','];
      if bpdata(3,j)>=0, s=[' ' s];end    
      if abs(bpdata(3,j))<10, s=[' ' s];end    
      if abs(bpdata(3,j))<100, s=[' ' s];end    
      while length(s)<7, s=[s ' ']; end
      amp2 = [amp2  s];
      if j/6==round(j/6), amp2=[amp2 10 32]; end
   end
   smp2='';
   for j=1:size(bpnames,1)
      s = [num2str(bpdata(4,j)) ','];
      if bpdata(4,j)>=0, s=[' ' s];end    
      if abs(bpdata(4,j))<10, s=[' ' s];end    
      while length(s)<7, s=[s ' ']; end
      smp2 = [smp2  s];
      if j/6==round(j/6), smp2=[smp2 10]; end
   end
   patmp2=[num2str(length(tok_data_struct.bpnames)),'*0.0'];

   flnam=' ';
   for j=1:size(flnames,1)
      flnam = [flnam '''' deblank(flnames(j,:)) ''','];
      if j/5==round(j/5), flnam=[flnam 10 32]; end
   end
   rsi='  ';
   for j=1:size(flnames,1)
      s = num2str(fldata(2,j)); while length(s)<5, s=[s '0']; end
      rsi = [rsi  s ', '];
      if j/6==round(j/6), rsi=[rsi 10 32 32]; end
   end
   zsi=' ';
   for j=1:size(flnames,1)
      s = num2str(fldata(1,j));
      if fldata(1,j)==0, s='0.';end
      if fldata(1,j)>=0, s=[' ' s];end    
      while length(s)<6, s=[s '0']; end
      zsi = [zsi  s ','];
      if j/6==round(j/6), zsi=[zsi 10 32]; end
   end


   fcdat = '';  
   for ii=[1:size(fcdata,2)]
      fcdat = strvcat(fcdat,sprintf('%-12.4f%-12.4f%-12.4f%-12.4f%-12.4f%-12.4f\n',fcdata([2 1 4 3 5 6],ii)));
   end
   if use_vv
      for ii=[1:nv]
         fcdat = strvcat(fcdat,sprintf('%-12.4f%-12.4f%-12.4f%-12.4f%-12.4f%-12.4f\n',vvdata([2 1 4 3 5 6],ii)));
      end
   end


   % a holds the characters that go in the file
   a =  [ ...
      ' $IN3' 10 ...
      ' RLEFT=',num2str(min(rg)),'     RRIGHT=',num2str(max(rg)) 10 ...
      ' ZBOTTO=',num2str(min(zg)),'    ZTOP=',num2str(max(zg)), ...
      ' NSMP2=25' 10 ...
      ' IGRID=1 IFCOIL=1 IVESEL=0 ISLPFC=1 IECOIL=0 IACOIL=0' 10];

   if ~exist('fcid','var')
      fcid = [1:tok_data_struct.nc];
   end
   if ~use_vv
      vvid = [];
   end
   if ~exist('vvid','var')
      vvid = [1:tok_data_struct.nv];
   end

   [FCID, FCTURN, TURNFC] = make_fc_strings(tok_data_struct,fcid,vvid);

   
   a = [a, ... 
      FCID ...
      FCTURN ...
      TURNFC ...       
      ' LPNAME=' 10 ...
      flnam ...
      '' 10 ...
      ' MPNAM2=' 10 ...
      bpnam ...
      '' 10 ...
      ' XMP2=' 10 ...
      xmp2 ...
      '' 10 ...
      ' YMP2=' 10 ...
      ymp2 ...
      '' 10 ...
      ' RSI=' 10 ...
      rsi ...
      '' 10 ...
      ' ZSI=' 10 ...
      zsi ...
      '' 10 ...
      '' 10 ...
      ' AMP2=' 10 ...
      amp2 ...
      '' 10 ...
      ' SMP2= ' 10 ...
      smp2 ...
      '' 10 ...
      [' PATMP2  =  ' patmp2 10] ...
      '' 10 ...
      '' 10 ...
      ' $' 10 ...
      ];


   f=fopen('mhdin.dat','w');
   fwrite(f,a);
   fwrite(f,fcdat');
   fclose(f);
   disp('Wrote file mhdin.dat')
end


function [fcidStr, fcturnStr, turnfcStr] = make_fc_strings(tds,fcid,vvid)
   allid  = [fcid max(fcid)+vvid];
   allTurn = [tds.ccnturn' ones(size(vvid))];
   %uniqueGroups = unique(abs(fcid));
   uniqueGroups = unique(abs(allid));
   nGroup = length(uniqueGroups);
   nTurnGroup = zeros(1,nGroup);
   for ii=1:nGroup
      idxGroup = find(abs(allid) == uniqueGroups(ii));
      nTurnGroup(ii) = sum(allTurn(idxGroup));
   end
   for ii=1:length(allid)
      fcturn(ii) = sign(allid(ii))*allTurn(ii)/nTurnGroup(abs(allid(ii)));
   end

   fcidStr = ' FCID=';
   for ii=allid(:)'
      fcidStr = [fcidStr num2str(ii) '.,'];
   end
   fcidStr = [fcidStr 10];
   
   fcturnStr = ' FCTURN=';
   for ii=1:length(fcturn)
      if fcturn(ii) < 1
         fcturnStr = [fcturnStr num2str(abs(fcturn(ii))) ','];
      else
         fcturnStr = [fcturnStr num2str(abs(fcturn(ii))) '.,'];
      end
   end
   fcturnStr = [fcturnStr 10];

   turnfcStr = ' TURNFC=';    
   for ii=1:nGroup
      if ismember(uniqueGroups(ii),fcid)
         turnfcStr = [turnfcStr num2str(nTurnGroup(ii)) '.,'];
      else
         turnfcStr = [turnfcStr '1.,'];
      end
   end
   turnfcStr = [turnfcStr 10];
end


