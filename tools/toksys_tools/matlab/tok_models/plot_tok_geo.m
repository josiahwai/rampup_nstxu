function plot_tok_geo(tok_data_struct,options,equil_data)
 %
%  USAGE:  plot_tok_geo(tok_data_struct,options,equil_data)
%
%	Examples: 
%		plot_tok_geo(tok_data_struct,options,east_system.equil_data)
%		plot_tok_geo(tok_data_struct,options,psizr)
%
%  PURPOSE: Plot the basic tokamak geometry.
%
%  INPUTS: [default]
%    tok_data_struct = structure containing tokamak geometry data
%    options = structure containing any of the following options:
%	iblackbg = flag: [0]= white figure background, 1=black figure background
%       ipltpsi  = flag: 1= contour equilibrium psizr, >1= #contours to plot
%				(default = 1 if psizr exists, else 0)
%       ipltB    = flag: 1= contour equilibrium |B|, >1= #contours to plot [0]
%       ipltJ    = flag: 1= contour equilibrium jphi, >1= #contours to plot [0]
%     (Only one of ipltpsi, ipltB, ipltJ may be specified; all need equil_data.)
%	ipltlim  = flag: 1= plot limiter switch [1]
%	ipltfl   = flag: 1= plot Flux Loops [0]
%	ipltbp   = flag: 1= plot B-Probe indices [0]
%	ipltgap  = flag: 1= plot gap locations [0]
%	ilabeleq = flag: 1= contour label flux values [0]
%	ilabelfc = flag: 1= label PF's with indices, 
%                       >1 => use name labels, with fontsize=ilabelfc [0]
%	ilabelcc = flag: 1= label PF coil STATES with indices (cannot also set ilabelfc) 
%       Pcc      = matrix that computes PF currents from states (required if ilabelcc=1)
%	ilabelvv = flag: 1= label VV elements with indices [0]
%	ilabelfl = flag: 1= label FL's with indices [0]
%                       >1 => use name labels, with fontsize=ilabelfl [0]
%	ilabelbp = flag: 1= label BP's with indices,
%                       >1 => use name labels, with fontsize=ilabelbp [0]
%	ilabelgap= flag: 1= label gaps with indices [0]
%	idxvv    = indices of vacuum vessel to plot (default = all)
%	idxfl    = indices of flux loops to plot (default all, used only if ipltfl)
%	idxbp    = indices of Bprobes to plot (default all, used only if ipltbp)
%	vvgroup  = grouping vector for vacuum vessel elements, used only if 
%			ilabelvv=1 to define vessel element labels
%       ileftright: bit 1 = plot left side, bit 0 = plot right side (default)
%    equil_data  = equilibrium psizr from EFIT OR structure containing
%		  psizr and psibnd (e.g. equil_data in *_system struct
%		  produced by build_*_sys.m script. 
%
%  OUTPUTS:
%	Plot of tokamak geometry
 
%  RESTRICTIONS:
%
%  VERSION:  
%  @(#)plot_tok_geo.m	1.22 09/09/11
%
%  WRITTEN BY:  Dave Humphreys  ON	3/27/04
%
%  MODIFICATION HISTORY:
%    jal08aug07 - improve ploting features and fix bugs
%    MLW - made device independent using generic data structures 4/20/06
%    DAH - added test for structure containing psizr and psibnd, so will
%		plot plasma boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Complicated logic here to handle all possibilities: Each iplt* variable may
% or may not be set.  If it is set, then only illegal combination is two or
% more of them being set > 0 (=0 is OK).  Also need to handle:
% - If none of the iplt* variables > 1, but 3rd argument is input, then
%   set ipltpsi = 1.
% - For backward compatibility, interpret iplteq as ipltpsi.

  ipltB= 0;  ipltJ= 0;  
  if nargin > 1
    if(isnumeric(options) & options==0)
       return;
    end
    if(~isempty(options))
       struct_to_ws(options);
    end
    if(exist('ilabelfc')~=1)
        ilabelfc = 0;
    end
    if(exist('ilabelcc')~=1)
        ilabelcc = 0;
    end
    if(ilabelfc & ilabelcc)
       wait('ERROR plot_tok_geo: cannot set both ilabelfc and ilabelcc')
       return;
    end
    if(ilabelcc & exist('Pcc')~=1)
       wait('ERROR plot_tok_geo: Pcc must be provided if ilabelcc is set')
       return;
    end
    if(exist('ipltpsi')~=1 & ipltB==0 & ipltJ==0)
      if(exist('equil_data')==1)
        ipltpsi = 1; 
      else
        ipltpsi = 0;
      end
    end
    if(exist('iplteq')==1 & exist('ipltpsi')~=1)
      ipltpsi = iplteq;
    end
    if(exist('ipltpsi')~=1)
      ipltpsi = 0;
    end
    plot_options= (ipltpsi>0)+(ipltB>0)+(ipltJ>0);
    if(plot_options > 1)
      wait('ERROR plot_tok_geo: only one of ipltpsi, ipltB, or ipltJ allowed')
      return;
    elseif(plot_options > 0 & exist('equil_data')~=1)
      wait('ERROR plot_tok_geo: specifying contour option requires equil_data')
      return;
    end
  else
    ipltpsi= 0;  
    ipltB= 0;  
    ipltJ= 0;  
  end
  if(exist('ilabeleq')~=1)
     ilabeleq= 0; 
  end
  if ~exist('ileftright','var'), ileftright = 1; end

  struct_to_ws(tok_data_struct);

  if nargin > 2
      if isstruct(equil_data) % then extract psizr, psibnd from structure
         if isfield(equil_data,'psizr') % test if psizr exists in structure
            psizr = equil_data.psizr;
	    if isfield(equil_data,'rg') rgefit= equil_data.rg; end
	    if isfield(equil_data,'zg') zgefit= equil_data.zg; end
         end
         if isfield(equil_data,'psibnd')
            psibnd = equil_data.psibnd;
         elseif isfield(equil_data,'psibry')
            psibnd = equil_data.psibry;
            disp('WARNING plot_tok_geo: psibnd is missing - substituting psibry')
         end
         if isfield(equil_data,'nbbbs')
            rbbbs= equil_data.rbbbs(1:equil_data.nbbbs);
            zbbbs= equil_data.zbbbs(1:equil_data.nbbbs);
         end
      else
         psizr = equil_data;    %assume equil_data is the psizr array
      end
  end

if(exist('vvgroup')~=1)
   vvgroup = [1:size(vvdata,2)]';
end
%Plotting parameters:
   dumfl = 0.02; % amount to displace FL index label rel to FL if ilabelfl=1

if exist('iblackbg')~=1
   iblackbg=0;
end

if exist('ipltlim')~=1
   ipltlim=1;
end

% ========================================================================
gcf;

if iblackbg % force background to be black in gcf
  whitebg(gcf,'k');
else
  whitebg(gcf,'w');
end 
                                               
if exist('fcdata')==1&~isempty(fcdata)
  if exist('idxfc')~=1 | isempty(idxfc)
     idxfc = 1:size(fcdata,2);
  end
  if mod(ileftright,4)>1 % Plot left side
    plot_boxx(-fcdata(2,idxfc),fcdata(1,idxfc),fcdata(4,idxfc),fcdata(3,idxfc),...
	'r',-fcdata(5,idxfc),-fcdata(6,idxfc))
    hold on
  end
  if mod(ileftright,2) % Plot right side (default)
    plot_boxx(fcdata(2,idxfc),fcdata(1,idxfc),fcdata(4,idxfc),fcdata(3,idxfc),...
	'r',fcdata(5,idxfc),fcdata(6,idxfc))
  end
end
hold on
if exist('ecdata')==1&~isempty(ecdata)
  if mod(ileftright,4)>1 % Plot left side
    plot_boxx(-ecdata(2,:),ecdata(1,:),ecdata(4,:),ecdata(3,:),'g')
  end
  if mod(ileftright,2) % Plot right side (default)
    plot_boxx(ecdata(2,:),ecdata(1,:),ecdata(4,:),ecdata(3,:),'g')
  end
end
if exist('vvdata')==1&~isempty(vvdata)
  if exist('idxvv')~=1 | isempty(idxvv)
     idxvv = 1:size(vvdata,2);
  end
  if mod(ileftright,4)>1 % Plot left side
    plot_box(-vvdata(2,idxvv),vvdata(1,idxvv),vvdata(4,idxvv),vvdata(3,idxvv),...
	'g',-vvdata(5,idxvv),-vvdata(6,idxvv))
  end
  if mod(ileftright,2) % Plot right side (default)
    plot_box(vvdata(2,idxvv),vvdata(1,idxvv),vvdata(4,idxvv),vvdata(3,idxvv),...
	'g',vvdata(5,idxvv),vvdata(6,idxvv))
  end
end
   
% Convention for limdata is 2 x n, each column = (z,r).
if(size(limdata,2)==2) limdata = limdata'; end;	% for backward-compabitility
if ipltlim & ~isempty(limdata)
  if mod(ileftright,4)>1 % Plot left side
    plot(-limdata(2,:),limdata(1,:),'b')
  end
  if mod(ileftright,2) % Plot right side (default)
    plot(limdata(2,:),limdata(1,:),'b')
  end
end

axis('image')
axis0=axis;  %save the axis settings for blowups and things...

if(length(ipltpsi) > 1 | ipltpsi > 0)
   if ~exist('rgefit'), rgefit=rg; end
   if ~exist('zgefit'), zgefit=zg; end
   [nzeq,nreq]= size(psizr);
   if nreq~=length(rgefit)
     wait(['%Caution: Plasma Model/equil R-grid different: ', ...
            int2str(length(rgefit)) '/' int2str(nreq),' resizing for plot only'])
     rgefit= linspace(rgefit(1),rgefit(end),nreq);
   end
   if nzeq~=length(zgefit)
     wait(['%Caution: Plasma Model/equil Z-grid different: ', ...
            int2str(length(zgefit)) '/' int2str(nzeq),' resizing for plot only'])
     zgefit= linspace(zgefit(1),zgefit(end),nzeq)';
   end        
   if(length(ipltpsi>1))
      cntrs = ipltpsi;
   elseif ipltpsi>=2 
      cntrs = linspace(min(psizr(:)),max(psizr(:)),ipltpsi);
   else
      cntrs = linspace(min(psizr(:)),max(psizr(:)),30);
   end
   if mod(ileftright,4)>1 % Plot left side
     [c,h] = contour(-rgefit,zgefit,psizr,cntrs,'r:');
   end
   if mod(ileftright,2) % Plot right side (default)
     [c,h] = contour(rgefit,zgefit,psizr,cntrs,'r:');
   end
   if(ilabeleq), clabel(c,h); end
   if exist('psibnd')==1
     if mod(ileftright,4)>1 % Plot left side
       [c,h]=contour(-rgefit,zgefit,psizr,[psibnd psibnd],'r');
       if(ilabeleq), clabel(c,h); end
     end
     if mod(ileftright,2) % Plot right side (default)
       [c,h]=contour(rgefit,zgefit,psizr,[psibnd psibnd],'r');
       if(ilabeleq), clabel(c,h); end
     end
   end
   if exist('rbbbs')==1 & exist('zbbbs')==1
     if mod(ileftright,4)>1 % Plot left side
       plot(-rbbbs,zbbbs,'k:')
     end
     if mod(ileftright,2) % Plot right side (default)
       plot(rbbbs,zbbbs,'k:')
     end
   end
end

if(ipltB > 0)
   if ~exist('rgefit'), rgefit=rg; end
   if ~exist('zgefit'), zgefit=zg; end
   [nzeq,nreq]= size(psizr);
   [dfdr,dfdz]= gradient(psizr,rgefit(2)-rgefit(1),zgefit(2)-zgefit(1));
   brr= -dfdz./rgg/(2*pi); % br= 1/(2pi*r)*df/dz
   bzz= +dfdr./rgg/(2*pi);
   bbb= sqrt(brr.^2 + bzz.^2);
   bbb = bbb*1e4;
   ncont= 30;
   if ipltB>=2 ncont= ipltB; end
   [c,h] = ...
     contour(rgefit,zgefit,bbb,linspace(min(bbb(:)),max(bbb(:)),ncont),'r:');
   if(ilabeleq), clabel(c,h); end
end

if(ipltJ > 0)
   if ~exist('rgefit'), rgefit=rg; end
   if ~exist('zgefit'), zgefit=zg; end
   jphi = equil_data.jphi;
   ncont= 30;
   if ipltJ>=2 ncont= ipltJ; end
   [c,h] = ...
     contour(rgefit,zgefit,jphi,linspace(min(jphi(:)),max(jphi(:)),ncont),'r:');
   if(ilabeleq), clabel(c,h); end
end

if exist('subtitle')
 title([upper(tokamak) ' Geometry - ' subtitle],'FontSize',13)
else
 title([upper(tokamak) ' Geometry'],'FontSize',13)
end

xlabel('R [m]')
ylabel('Z [m]')

if(~exist('idxfl','var'))
   idxfl = 1:size(fldata,2);
end
if exist('ipltfl')==1&ipltfl&exist('fldata')
  plot(fldata(2,idxfl),fldata(1,idxfl),'mo')
end
if(~exist('idxbp','var'))
   idxbp = 1:size(bpdata,2);
end
if exist('ipltbp')==1&ipltbp&exist('bpdata') & length(idxbp)>0
  plot(bpdata(2,idxbp),bpdata(1,idxbp),'cs')
  RMPI=bpdata(2,idxbp)'; ZMPI=bpdata(1,idxbp)'; 
  tltMPI=pi*bpdata(3,idxbp)'/180; mpilen=bpdata(4,idxbp)';
%Plotting from *vertical*:
%   h=quiver(RMPI-0.5*mpilen.*sin(tltMPI),ZMPI-0.5*mpilen.*cos(tltMPI), ...
%      mpilen.*sin(tltMPI),mpilen.*cos(tltMPI),0);
%Plotting from outboard horizontal:
   h=quiver(RMPI-0.5*mpilen.*cos(tltMPI),ZMPI-0.5*mpilen.*sin(tltMPI), ...
      mpilen.*cos(tltMPI),mpilen.*sin(tltMPI),0);
end

if exist('ipltgap')==1&ipltgap&exist('gapdata')
  plot(gapdata(2,:),gapdata(1,:),'m*')
end

if exist('ilabelfc')==1&ilabelfc & exist('fcdata')
  if(ilabelfc==1)
     for ii=1:size(fcdata,2)
       dum = 0.5*min(fcdata(4,ii),fcdata(3,ii))*0;
       h=text(fcdata(2,ii)+dum,fcdata(1,ii)+dum,int2str(ii))*0;
       set(h,'FontSize',12,'Color','r');
     end
  elseif(ilabelfc>1 & exist('fcnames'))
     for ii=1:size(fcdata,2)
       dum = 0.5*min(fcdata(4,ii),fcdata(3,ii))*0;
       h=text(fcdata(2,ii)+dum,fcdata(1,ii)+dum,fix_undscr(deblank(fcnames(ii,:))));
       set(h,'FontSize',ilabelfc,'Color','r');
     end
  end
end

if exist('ilabelcc')==1&ilabelcc & exist('fcdata')
     if(isfield(tok_data_struct,'ecnturn'))
        nec = length(ecnturn);
     else
        nec = 0;
     end
     for ii=1:size(fcdata,2)
       dum = 0.5*min(fcdata(4,ii),fcdata(3,ii))*0;
       jj = find(Pcc(nec+ii,:));
       if(length(jj)==1)
          h=text(fcdata(2,ii)+dum,fcdata(1,ii)+dum,int2str(jj));
          set(h,'FontSize',12,'Color','r');
       end
     end
end

if exist('ilabelvv')==1&ilabelvv&exist('vvdata')
  for ii=idxvv
    dum = 0.5*min(vvdata(4,ii),vvdata(3,ii))*0;
    h=text(vvdata(2,ii),vvdata(1,ii),int2str(vvgroup(ii,1)),'hori','ce');
    set(h,'FontSize',10,'Color','g');
  end
end

if exist('ilabelfl')==1&ilabelfl&exist('fldata')
  if(ilabelfl==1)
     for ii=1:length(idxfl)
       dum = 0.5*min(fldata(4,idxfl(ii)),fldata(3,idxfl(ii)))*0;
       h=text(fldata(2,idxfl(ii))+dum,fldata(1,idxfl(ii))+dum,int2str(idxfl(ii)));
       set(h,'FontSize',10,'Color','m');
     end
  elseif(ilabelfl>1 & exist('flnames'))
     for ii=1:length(idxfl)
       dum = 0.5*fldata(4,idxfl(ii))*0;   %use mpilen for displacement of label text
       h=text(fldata(2,idxfl(ii))+dum,fldata(1,idxfl(ii))+dum,fix_undscr(deblank(flnames(idxfl(ii),:))));
       set(h,'FontSize',ilabelfl,'Color','m');
     end
  end
end

if exist('ilabelbp')==1&ilabelbp&exist('bpdata')
  if(ilabelbp==1)
     for ii=1:length(idxbp)
       dum = 0.5*bpdata(4,idxbp(ii))*0;   %use mpilen for displacement of label text
       h=text(bpdata(2,idxbp(ii))+dum,bpdata(1,idxbp(ii))+dum,int2str(idxbp(ii)));
       set(h,'FontSize',10,'Color','c');
     end
  elseif(ilabelbp>1 & exist('bpnames'))
     for ii=1:length(idxbp)
       dum = 0.5*bpdata(4,idxbp(ii))*0;   %use mpilen for displacement of label text
       h=text(bpdata(2,idxbp(ii))+dum,bpdata(1,idxbp(ii))+dum,fix_undscr(deblank(bpnames(idxbp(ii),:))));
       set(h,'FontSize',ilabelbp,'Color','c');
     end
  end
end

if exist('ilabelgap')==1&ilabelgap&exist('gapdata')
  for ii=1:size(gapdata,2)
    dum = 0.5*min(gapdata(4,ii),gapdata(3,ii))*0;
    h=text(gapdata(2,ii)+dum,gapdata(1,ii)+dum,int2str(ii));
    set(h,'FontSize',10,'Color','m');
  end
end


hold off
