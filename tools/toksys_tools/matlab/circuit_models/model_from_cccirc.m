function model = model_from_cccirc(vacuum_objs,cccirc,vvcirc,vvgroup,Mmat,Rvec,iplcirc,verbose)
%
%  SYNTAX:  model = model_from_cccirc(vacuum_objs,cccirc,vvcirc,vvgroup,Mmat,Rvec,iplcirc,verbose)
%
%  PURPOSE:  Create circuit model using cccirc input.
%
%  INPUT:
%    vacuum_objs = vacuum model objects (often called tok_data_struct)
%    cccirc      = coil conductor currents circuit description
%    vvcirc      = vacuum vessel conductors circuit description
%    vvgroup     = vacuum vessel element groupings
%    Mmat        = mutual inductances (matrix) for all conductors
%    Rvec        = resistances (vector) for all conductors
%    iplcirc     = if 1, include plasma circuit (in which case Mmat, Rvec must contain plasma elements), else 0
%    verbose     = (optional)
%
%  OUTPUT:  data struct model containing:
%    Pxx  =
%    Pcc  =
%    Mhat =
%    Vhat =
%    Rhat =
% 
%  RESTRICTIONS:  Modeling with cccirc no longer supports use of Lckt_extra or Rckt_extra.
 
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	2/4/19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   if(nargin<8)
      verbose=0;
   end
   if(isempty(vvcirc))
      vvcirc = [1:vacuum_objs.nv]';
   end
   if(isempty(vvgroup))
      vvgroup = [1:vacuum_objs.nv]';
   end

% Check for consistency of specified vvgroup and vvcirc
      temp = vvgroup(:,1);
      temp = setdiff(unique(temp),0);
      if(length(temp) ~= max(temp))
         wait('ERROR model_from_cccirc: vvgroup must define contiguous numbered groups')
         return;
      end
         temp1 = setdiff(union(temp,abs(vvcirc)),0);
         if(length(temp1)~=length(temp))
            wait('ERROR model_from_cccirc: vvcirc and vvgroup inconsistent')
            return;
         end

      resv = vacuum_objs.resv;
      if(size(vvgroup,2)>2)
        wait('ERROR model_from_cccirc: vvgroup must be 1 or 2 column matrix')
        return;
      end
      if(size(vvgroup,1)~=length(resv))
         wait(['ERROR model_from_cccirc: vvgroup col. size must match '...
                'size of resv = ' int2str(length(resv))])
         return;
      end
      nvx = max(abs(vvgroup(:,1)));      %Find reduced size of VV
   nvx = length(setdiff(unique(abs(vvcirc)),0));

   ncc = vacuum_objs.nc;
   nvv = vacuum_objs.nv;
   if(isfield(vacuum_objs,'ecdata'))
      ecdata  = vacuum_objs.ecdata;
   else
      ecdata = [];
   end
   if(isempty(ecdata))
      necoils=0;
      ecnturn = [];
   else
      ecnturn = vacuum_objs.ecnturn;
      necoils = length(ecnturn);
   end

   if(~isempty(cccirc)) 	% Circuit connections selected by cccirc
      ncx = max(abs(cccirc));
      if verbose>0
         disp('Warning: CC circuit connections being modified by cccirc!!!')
      end
      Pcc = zeros(ncc,ncx);  %Pcc maps from conn-circs to orig unconn-cir

      nturn = [ecnturn;vacuum_objs.fcnturn];
      for ii=1:ncx
        idx0=find(abs(cccirc)==ii);
        idx1=find(cccirc==ii);
        idx2=find(cccirc==-ii);
        if(vacuum_objs.iterminal)
           turns_ratio1 = 1;
           turns_ratio2 = 1;
        else
           turns_ratio1 = nturn(idx1)/sum(nturn(idx0));
           turns_ratio2 = nturn(idx2)/sum(nturn(idx0));
        end
        if ~isempty(idx1)
           Pcc(idx1,ii)=ones(length(idx1),1).*turns_ratio1;
        end
        if ~isempty(idx2)
           Pcc(idx2,ii)=-ones(length(idx2),1).*turns_ratio2;
        end
      end
   else
      Pcc = eye(ncc);
   end

   if(0) 	%   if(make_netlist)
      temp_file = [tokamak '_netlist.dat'];
      [netlist,inode1,inode2] = cccirc_to_netlist(cccirc,vacuum_objs,temp_file);
   end

%Pxx maps from full state vector *with* connected circs to unconn-circs:

   Pvv = eye(nvv);  
     if verbose > 1
        disp('VV circuit connections being modified by vvgroup...')
     end
     if(size(vvgroup,2)==1)       % distribute current according to resistance
        calc_vvfrac = 1;
        vvfrac = zeros(length(vvgroup),1);
     elseif(size(vvgroup,2)==2)   % user-specified sharing of current
        calc_vvfrac = 0;
        vvfrac = vvgroup(:,2);
     end
     vvid = vvgroup(:,1);
     nvg =max(abs(vvid));	% nvg = number of grouped vv elts
     Pvv = zeros(nvv,nvg);  
  
     if(iplcirc)
        resv = Rvec(ncc+1:end-1);
     else
        resv = Rvec(ncc+1:end);
     end

     for ii=1:nvg
        idx1=find(vvid==ii); 
        if ~isempty(idx1)
           if(calc_vvfrac)
              sum_rinv = sum(1./resv(idx1));
              vvfrac(idx1) = 1./resv(idx1)/sum_rinv;
           elseif(abs(sum(vvfrac(idx1))-1)>1e-3);
              fprintf(['WARNING model_from_cccirc: vessel group current ' ...
                        'sharing does not sum to 1\n  vessel indices = ']);
              fprintf('%d ',idx1);
              fprintf('\n');
           end
           Pvv(idx1,ii)=vvfrac(idx1);
        end
        idx1=find(vvid==-ii);
        if ~isempty(idx1)
           if(calc_vvfrac)
              sum_rinv = sum(1./resv(idx1));
              vvfrac(idx1) = 1./resv(idx1)/sum_rinv;
           elseif(abs(sum(vvfrac(idx1))-1)>1e-3);
              fprintf(['WARNING model_from_cccirc: vessel group current ' ...
                        'sharing does not sum to 1\n  vessel indices = ']);
              fprintf('%d ',idx1);
              fprintf('\n');
           end
           Pvv(idx1,ii)=-vvfrac(idx1);
        end
     end

      temp = zeros(nvg,nvx);
      for ii=1:nvx
        idx1=find(vvcirc==ii);
        idx2=find(vvcirc==-ii);
        if ~isempty(idx1)
           temp(idx1,ii)=ones(length(idx1),1);
        end
        if ~isempty(idx2)
           temp(idx2,ii)=-ones(length(idx2),1);
        end
      end
      Pvv = Pvv * temp;
   Pxx = [[Pcc zeros(ncc,nvx)]; [zeros(nvv,ncx) Pvv]];

   if(iplcirc)
      Pxx = [[Pxx zeros(size(Pxx,1),1)]; [zeros(1,size(Pxx,2)) 1]];
   end

   Mhat = Pxx'*Mmat*Pxx;
   Rhat = Pxx'*diag(Rvec)*Pxx;

% error checking:
   temp1 = min(Rhat); [mm,ii] = min(temp1);
   if mm<0
      disp('ERROR model_from_cccirc: ')
      disp(['Consolidation of resistances has ' ...
        'resulted in NEGATIVE resistance value for circuit ' int2str(ii)])
      return;
   end

   nrows = size(Mhat,1);
   Vhat=zeros(nrows,ncx); Vhat(1:ncx,1:ncx)=eye(ncx);  %drive only circuits...

   model = struct( ...
      'Pxx',Pxx, ...
      'Pcc',Pcc, ...
      'ncx',ncx, ...
      'Mhat', Mhat, ...
      'Rhat', Rhat, ...
      'Vhat', Vhat);
%      'rxx', rxx, ...

