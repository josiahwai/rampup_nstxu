%
%  USAGE:   gs_edge_response
%
%  PURPOSE: Calculate how xedge responds to changes in rbbbs, zbbbs
%
%  INPUTS: nedge, fedge, (from gs_trace_edge)
%          DX = diff(rbbbs)', DY = diff(zbbbs)'
%          nz, grid variable
%
%  OUTPUTS: dxedgedrbbbs0, dxedgedzbbbs0, derivatives w.r.t. point where xedge = 0
%           dxedgedrbbbs1, dxedgedzbbbs1, derivatives w.r.t. point where xedge = 1

%  METHOD:  
	
%  VERSION @(#)gs_edge_response.m	1.1 10/07/13
%
%  WRITTEN BY:  Anders Welander  ON	3/12/13
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dxedgedrbbbs0 = zeros(nbbbs_max,1);
dxedgedzbbbs0 = zeros(nbbbs_max,1);
dxedgedrbbbs1 = zeros(nbbbs_max,1);
dxedgedzbbbs1 = zeros(nbbbs_max,1);

for j = 1:nedge
  k = iedge(j);
  if fedge(j) == -nz | fedge(j) == nz
    dxedgedrbbbs0(j) = -(1-(redge(j)-rbbbs(k))/DX(k))/DX(k);
    dxedgedrbbbs1(j) = (rbbbs(k)-redge(j))/DX(k)^2;
    dxedgedzbbbs0(j) = 0;
    dxedgedzbbbs1(j) = 0;
  elseif fedge(j) == -1 | fedge(j) == 1
    dxedgedrbbbs0(j) = 0;
    dxedgedrbbbs1(j) = 0;
    dxedgedzbbbs0(j) = -(1-(zedge(j)-zbbbs(k))/DY(k))/DY(k);
    dxedgedzbbbs1(j) = (zbbbs(k)-zedge(j))/DY(k)^2;
  end
end
dxedgedrbbbs0(nedge+1) = dxedgedrbbbs0(1);
dxedgedrbbbs1(nedge+1) = dxedgedrbbbs1(1);
dxedgedzbbbs0(nedge+1) = dxedgedzbbbs0(1);
dxedgedzbbbs1(nedge+1) = dxedgedzbbbs1(1);
