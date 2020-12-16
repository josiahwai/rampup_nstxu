function vst_data = plasma_response_vst(eq_file,tokamak,tok_data_struct,ichooseq)
%
%  USAGE:  vst_data = plasma_response_vst(equil_data,iefit,tokamak,tok_data_struct,ichooseq)
%
%  PURPOSE: Calculate plasma response objects using vst in corsica
%
%  INPUTS:
%    eq_file = name of an efit gfile or a saved corsica equilibrium
%    tokamak = device to construct model objects for (e.g. 'NSTX','KSTAR',etc)
%    vac_objs = structure containing: 
%		mcc, mvv, mcv, mpc, mpv, resc, resv, zg, rg, ecdata
%     		(Note that imks and iterminal in this structure define the
%		 units of the data objects to be produced. imks=1 gives MKS
%		 units, otherwise units are MA,uH,uOhms.  iterminal=1 gives
%		 terminal mode, 0 gives lumped.)
%    ichooseq= equilibrium file type:
%	1 = efit_gfile
%	2 = corsica generated flat files
%	3 = saved corsica equilibrium
%
%  OUTPUTS:
%     vst = structure containing:
%
%  METHOD:  The matlab script generates a basis script file and invokes corsica
%	    with this file to generate plasma response objects with vst.
%	    After that execution returns to matlab.
%
%  VERSION: @(#)plasma_response_vst.m	1.2 06/25/09

%  WRITTEN BY:  Anders Welander  ON	9/29/08
%
%  MODIFICATION HISTORY:
%	
%  NOTES:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if exist('tokamak') ~= 1, tokamak = 'D3D'; end
  if exist('tok_data_struct') ~= 1
    switch upper(tokamak)
      case 'DIII-D'
        load /m/GAtools/tokamaks/d3d/make/after_ADP/d3d_obj_mks_struct_6565.mat
      case 'D3D'
        load /m/GAtools/tokamaks/d3d/make/after_ADP/d3d_obj_mks_struct_6565.mat
      case 'ITER'
        load /m/GAtools/tokamaks/iter/make/standard/iter_obj_mks_struct.mat
      case 'EAST'
        load /m/GAtools/tokamaks/east/make/east_obj_mks_struct.mat
    end
  end
  struct_to_ws(tok_data_struct);
  if exist('ichooseq') ~= 1, ichooseq = 3; end
  
  fid = fopen('/m/GAtools/matlab/tok_models/create_objects.bas');
  a=fread(fid);
  fclose(fid);
  a = char(a');
  if ichooseq == 1
     a = strrep(a,'EQUILIBRIUM',['read d3.bas' 10 '  d3("' eq_file '",0)']);
  end
  if ichooseq == 3
     a = strrep(a,'EQUILIBRIUM',['read ' eq_file 10 '  run']);
  end
  b = ['nflux = ' num2str(nfl+nbp*2) 10];
  for j=1:nfl, b = [b '  rflux(' num2str(j) ') = ' num2str(fldata(2,j)) 10]; end % r of flux loops
  for j=1:nbp, b = [b '  rflux(' num2str(j+nfl) ') = ' num2str(bpdata(2,j)+1e-3*sin(pi/180*bpdata(3,j)),9) 10]; end % r+ of B probes
  for j=1:nbp, b = [b '  rflux(' num2str(j+nfl+nbp) ') = ' num2str(bpdata(2,j)-1e-3*sin(pi/180*bpdata(3,j)),9) 10]; end % r- of B probes
  for j=1:nfl, b = [b '  zflux(' num2str(j) ') = ' num2str(fldata(1,j)) 10]; end % r of flux loops
  for j=1:nbp, b = [b '  zflux(' num2str(j+nfl) ') = ' num2str(bpdata(1,j)-1e-3*cos(pi/180*bpdata(3,j)),9) 10]; end % z+ of B probes
  for j=1:nbp, b = [b '  zflux(' num2str(j+nfl+nbp) ') = ' num2str(bpdata(1,j)+1e-3*cos(pi/180*bpdata(3,j)),9) 10]; end % z- of B probes
  a = strrep(a,'DIAGNOSTICS',b);
  fid = fopen('temp_create_objects.bas','w');
  fwrite(fid,a);
  fclose(fid);
  !/d/caltrans/vcaltrans/bin/caltrans temp_create_objects.bas
  disp('Back in matlab')
   
% Load VST objects:
  if upper(tokamak) == 'D3D' | upper(tokamak) =='DIII-D', nc = 18; end

  load linduct.flat; linduct = reshape(linduct,nc+nv,nc+nv);
  mss = linduct([nv+1:nc+nv 1:nv],[nv+1:nc+nv 1:nv])*(.4e-6*pi);
  for j=1:nc, mss(j,:)=mss(j,:)*fcnturn(j); mss(:,j)=mss(:,j)*fcnturn(j); end

  load linduct2diag.flat; linduct2diag = reshape(linduct2diag,nc+nv,nfl+2*nbp)';
  diag = linduct2diag(:,[nv+1:nc+nv 1:nv])*(.4e-6*pi);
  for j=1:nc, diag(:,j)=diag(:,j)*fcnturn(j); end
  vst_data.mls = diag(1:nfl,:);
  for j=1:nbp
     vst_data.gbc(j,:) = (diag(nfl+j,:)-diag(nfl+nbp+j,:))/(2*pi*bpdata(2,j)*2e-3);
  end

  load lmatrix.flat; lmatrix = reshape(lmatrix,nc+nv,nc+nv);
  mxs = lmatrix([nv+1:nc+nv 1:nv],[nv+1:nc+nv 1:nv])*(.4e-6*pi);
  for j=1:nc, mxs(j,:)=mxs(j,:)*fcnturn(j); mxs(:,j)=mxs(:,j)*fcnturn(j); end

  load lmatrix2diag.flat; lmatrix2diag = reshape(lmatrix2diag,nc+nv,nfl+2*nbp)';
  diag = lmatrix2diag(:,[nv+1:nc+nv 1:nv])*(.4e-6*pi);
  for j=1:nc, diag(:,j)=diag(:,j)*fcnturn(j); end
  vst_data.dfldis = diag(1:nfl,:);
  for j=1:nbp
     vst_data.dbpdis(j,:) = (diag(nfl+j,:)-diag(nfl+nbp+j,:))/(2*pi*bpdata(2,j)*2e-3) - vst_data.gbc(j,:);
  end

  if upper(tokamak) == 'D3D' | upper(tokamak) =='DIII-D'
     nc = nc+2; % Adding back the E coils
     mss = [zeros(2,nc+nv); [zeros(nc+nv-2,2) mss]];
     mxs = [zeros(2,nc+nv); [zeros(nc+nv-2,2) mxs]];
     vst_data.dfldis = [zeros(nfl,2) vst_data.dfldis];
     vst_data.mls = [zeros(nfl,2) vst_data.mls];
     vst_data.dbpdis = [zeros(nbp,2) vst_data.dbpdis];
     vst_data.gbc = [zeros(nbp,2) vst_data.gbc];
  end
    
  vst_data.xmats = mxs - mss;

  % No file littering
  !rm .corsica-command-line temp_create_objects.* 
  !rm linduct.flat linduct2diag.flat lmatrix.flat lmatrix2diag.flat
  
  return
  
  % Tests (need to run efit_gfile='g130105.04500', EPS_switch_pos=1,gridsize=6565, build_d3d_sys first)
  
  figure(1),clf,hold
  plot(mss(5,:))
  plot(tok_data_struct.mcc(5,:),'g')
  plot(vst_data.xmats(5,:),'r')
  
  figure(2),clf,hold
  plot(d3d_system.mxx(5,1:48))
  plot(d3d_system.xmatx(5,1:48),'r')
  plot(tok_data_struct.mcc(5,1:20),'g')

  figure(3),clf,hold
  plot(tok_data_struct.mlc(5,1:20),'r+')
  plot(21:48,tok_data_struct.mlv(5,1:28),'r+')
  plot(vst_data.mls(5,:),'b')

  figure(4),clf,hold
  plot(tok_data_struct.gbc(:,3),'r+')
  plot(vst_data.gbc(:,3),'b')

  figure(5),clf,hold
  plot(d3d_system.dbpdis(:,15),'r+')
  plot(vst_data.dbpdis(:,15),'b')

  figure(6),clf,hold
  plot(d3d_system.dfldis(:,15),'r+')
  plot(vst_data.dfldis(:,15),'b')

















