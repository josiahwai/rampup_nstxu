function scaled_system = scale_equil_response(tok_system,scaleip)
 %
%  SYNTAX:  scaled_system = scale_equil_response(tok_system,scaleip)
%
%  PURPOSE:  Modify the response model contained in tok_system for an
%	equilibrium scaled by the factor scaleip.
%
%  INPUT:
%	tok_system = system model of the type built by build_tokamak_system
%	scale_ip = scalar to multiply by currents in equilibrium that was
%			used when generating tok_system
%
%  OUTPUT:
%	scaled_system = system model of the same type representing an
%			equilibrium scaled by the factor scaleip
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	7/1/2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)scale_equil_response.m	1.4 04/30/11

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES:
% (1) Scaling derived from running build_tokamak_system with scaled objects:
%   temp1.mat = original objects
%   temp2.mat = scaled by -1
%   temp3.mat = scaled by -2

% (2) Checked that this function reproduces result from scaled build_tokamak_system:
%	scaled_system = scale_equil_response(sys1,-2);
% Everything identical except:
%   fmat: identical except 4th column, where normalized difference < 1e-4
%   hmat: identical except 4th column, where normalized difference < 1e-5

% (3) FINALLY, checked that physics makes sense for all scalings:
%    CONFIRMED, except objects relating to betap, li (incl. hmat, fmat) not checked
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mony objects are unchanged:
 scaled_system = tok_system;

% scale those objects that do change:
 scaled_system.Ieq = scaled_system.Ieq * scaleip;
 scaled_system.cc0 = scaled_system.cc0 * scaleip;
 scaled_system.equil_data.brsp = scaled_system.equil_data.brsp * scaleip;
 scaled_system.equil_data.cc = scaled_system.equil_data.cc * scaleip;
 scaled_system.equil_data.cpasma = scaled_system.equil_data.cpasma * scaleip;
 scaled_system.equil_data.jphi = scaled_system.equil_data.jphi * scaleip;
 scaled_system.equil_data.pcurrt = scaled_system.equil_data.pcurrt * scaleip;
 scaled_system.equil_data.psibnd = scaled_system.equil_data.psibnd * scaleip;
 scaled_system.equil_data.psibry = scaled_system.equil_data.psibry * scaleip;
 scaled_system.equil_data.psimag = scaled_system.equil_data.psimag * scaleip;
 scaled_system.equil_data.psirz = scaled_system.equil_data.psirz * scaleip;
 scaled_system.equil_data.psizr = scaled_system.equil_data.psizr * scaleip;
 scaled_system.equil_data.ssibry = scaled_system.equil_data.ssibry * scaleip;
 scaled_system.equil_data.ssimag = scaled_system.equil_data.ssimag * scaleip;
 scaled_system.fmat = scaled_system.fmat * scaleip;
 scaled_system.rzrig_data.cphi = scaled_system.rzrig_data.cphi * scaleip;
 scaled_system.rzrig_data.dFrdbetap = scaled_system.rzrig_data.dFrdbetap * scaleip^2;
 scaled_system.rzrig_data.dFrdip =  scaled_system.rzrig_data.dFrdip * scaleip;
 scaled_system.rzrig_data.dFrdis = scaled_system.rzrig_data.dFrdis * scaleip;
 scaled_system.rzrig_data.dFrdli = scaled_system.rzrig_data.dFrdli * scaleip^2;
 scaled_system.rzrig_data.dFrdr = scaled_system.rzrig_data.dFrdr * scaleip^2;
 scaled_system.rzrig_data.dFrdz = scaled_system.rzrig_data.dFrdz * scaleip^2;
 scaled_system.rzrig_data.dFzdis = scaled_system.rzrig_data.dFzdis * scaleip;
% scaled_system.rzrig_data.dFzdip = scaled_system.rzrig_data.dFzdip * scaleip;
 scaled_system.rzrig_data.dFzdr = scaled_system.rzrig_data.dFzdr * scaleip^2;
 scaled_system.rzrig_data.dFzdz = scaled_system.rzrig_data.dFzdz * scaleip^2;
 scaled_system.rzrig_data.dbpdr = scaled_system.rzrig_data.dbpdr * scaleip;
 scaled_system.rzrig_data.dbpdz = scaled_system.rzrig_data.dbpdz * scaleip;
 scaled_system.rzrig_data.dcdr = scaled_system.rzrig_data.dcdr * scaleip;
 scaled_system.rzrig_data.dcdz = scaled_system.rzrig_data.dcdz * scaleip;
 scaled_system.rzrig_data.dfcdr = scaled_system.rzrig_data.dfcdr * scaleip;
 scaled_system.rzrig_data.dfcdz = scaled_system.rzrig_data.dfcdz * scaleip;
 scaled_system.rzrig_data.dfldr = scaled_system.rzrig_data.dfldr * scaleip;
 scaled_system.rzrig_data.dfldz = scaled_system.rzrig_data.dfldz * scaleip;
 scaled_system.rzrig_data.dflvdr = scaled_system.rzrig_data.dflvdr * scaleip;
 scaled_system.rzrig_data.dflvdz = scaled_system.rzrig_data.dflvdz * scaleip;
 scaled_system.rzrig_data.dfpdrC = scaled_system.rzrig_data.dfpdrC * scaleip;
 scaled_system.rzrig_data.dfpdzC = scaled_system.rzrig_data.dfpdzC * scaleip;
 scaled_system.rzrig_data.dfsdr = scaled_system.rzrig_data.dfsdr * scaleip;
 scaled_system.rzrig_data.dfsdz = scaled_system.rzrig_data.dfsdz * scaleip;
 scaled_system.rzrig_data.dfvdr = scaled_system.rzrig_data.dfvdr * scaleip;
 scaled_system.rzrig_data.dfvdz = scaled_system.rzrig_data.dfvdz * scaleip;
 scaled_system.rzrig_data.drdip = scaled_system.rzrig_data.drdip / scaleip;
 scaled_system.rzrig_data.drdis = scaled_system.rzrig_data.drdis / scaleip;
 scaled_system.rzrig_data.dzdis = scaled_system.rzrig_data.dzdis / scaleip;
 scaled_system.rzrig_data.dzdip = scaled_system.rzrig_data.dzdip / scaleip;
 scaled_system.rzrig_data.ip0 = scaled_system.rzrig_data.ip0 * scaleip;

% Objects for output equation - scaling depends on what kind of output.

% In cmat, responses will stay the same for: currents, flux loops, Bprobes,
% Rogowskis, and loop voltages, but will change for R and Z outputs
 idxR = strmatch('R',scaled_system.outputs,'exact');
 idxZ = strmatch('Z',scaled_system.outputs,'exact');
 scaled_system.cmat([idxR;idxZ],:) = scaled_system.cmat([idxR;idxZ],:) / scaleip;

 scaled_system.hmat = scaled_system.hmat * scaleip;
 scaled_system.hmat([idxR;idxZ],:) = scaled_system.hmat([idxR;idxZ],:) / scaleip;

 scaled_system.yeq = scaled_system.yeq * scaleip;
 scaled_system.yeq(idxR) = scaled_system.yeq(idxR)/scaleip;
 scaled_system.yeq(idxZ) = scaled_system.yeq(idxZ)/scaleip;

if 0	% Can we convert equil_data to a scaled equivalent?? Should we?
 scaled_system.source_equil_data = tok_system.equil_data;
 scaled_system.equil_data = tok_system.equil_data;
 scaled_system.equil_data.cpasma = scaleip*tok_system.equil_data.cpasma;
 scaled_system.equil_data.brsp = scaleip*tok_system.equil_data.brsp;
 scaled_system.equil_data.pcurrt = scaleip*tok_system.equil_data.pcurrt;
 scaled_system.equil_data.psirz = scaleip*tok_system.equil_data.psirz;
 scaled_system.equil_data.cc = scaleip*tok_system.equil_data.cc;
 scaled_system.equil_data.psizr = scaleip*tok_system.equil_data.psizr;
 scaled_system.equil_data.psimag = scaleip*tok_system.equil_data.psimag;
 scaled_system.equil_data.psibry = scaleip*tok_system.equil_data.psibry;
 scaled_system.equil_data.psibnd = scaleip*tok_system.equil_data.psibnd;
 scaled_system.equil_data.ssibry = scaleip*tok_system.equil_data.ssibry;
 scaled_system.equil_data.ssimag = scaleip*tok_system.equil_data.ssimag;
 scaled_system.equil_data.jphi = scaleip*tok_system.equil_data.jphi;
 %scaled_system.equil_data.ffprim = scaleip*tok_system.equil_data.ffprim;
end
