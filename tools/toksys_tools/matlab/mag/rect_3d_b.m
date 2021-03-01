  function [bx,by,bz]= rect_3d_b.m(x,y,z,a,b,c)

% PURPOSE:
% This function generates B-field in carthesian 3-space [X] == [x,y,z]
% produced by a uniform current density rectangle of dimensions [A] == [a,b,c]
% contining uniform current in the +Z direction and carrying unit current I=1A.
% The coordinate system is centered about the geometric center of rectangle.
% Note: All Units are MKS (meters, Amps, Tesla)
%
% SYNTAX:
%        [bx,by,bz]= rect_3d_b.m(x,y,z,a,b,c)
%
% INPUT:
%       [x,y,z] == [X]  vectors of coordinates, one point per row
%       [a,b,c] == [A]  vectors of rectangle dimensions, one rectangle per row
%
% Note: X and a should be the same size.
%       If A is a single row its rows are expanded to the rows of x
%       If X is a single row its rows are expanded to the rows of a
%
% OUTPUT:
%       [bx,by,bz] == [B] vectors of B-fields at point corresponding to x & a
%
% NOTE: Vectorized so can handle vectors of same length and scalers
% NOTE: Produces double precision solution every where except eps from edges
% CAUTION: Needs additional Validations at edge discontinuities (Part II)

% Jim Leuer, General Atomics, 12-97
% Uses NO Matlab_5 capabilities; coded for speed and space save

% ------------------------------------------------------------
% initialization & check input data integrity

  prog_name= 'rect_3d_b.m';

  if nargin <= 5
    disp(['%ERROR ',prog_name,': Must have at least 6 arguments']);
    eval(['help ',prog_name,]) % print out help from prog_name
    return
  end

  nnx= size(x);  nny= size(y);  nnz= size(z);

  if nnx(2) > 1
    disp(['%ERROR ',prog_name,': x,y,z must be scalor or vector not matrix'])
    help rect_3d_b.m
    return
  end

  nna= size(a);  nnb= size(b);  nnc= size(c);

  if nna(2) > 1
    disp(['%ERROR ',prog_name,': a,b,c must be scalor or vector not matrix'])
    help rect_3d_b.m
    return
  end

% Compare Overall length of x,y,z and make sure same as a,b,c {use x and a}
  if nnx(1) ~= nna(1)
    if nnx(1) < nna(1)                       % fill in x matrix rows to fit a
       x= [x;ones(nna(1)-nnx(1),1)*x(end,:)];
       y= [y;ones(nnb(1)-nny(1),1)*y(end,:)];
       z= [z;ones(nnc(1)-nnz(1),1)*z(end,:)];
       nnx(1)= nna(1); nny(1)= nnb(1); nnz(1)= nc(1);
    else                                   % fill in a matrix rows to fit x
       a= [a;ones(nnx(1)-nna(1),1)*a(end,:)];
       b= [b;ones(nny(1)-nnb(1),1)*b(end,:)];
       c= [c;ones(nnz(1)-nnc(1),1)*c(end,:)];
       nna(1)= nnx(1); nnb(1)= nny(1); nnc(1)= nnz(1);
    end
  end

  n= nna(1); % all input argument should now be same length of vectors
% ---------------------------------------------------------------------------
% PART I
% Standard computation
% ---------------------------------------------------------------------------

  pm= [0.5,-0.5]; % plus minus 0.5

  bx= zeros(n,1);
  by= zeros(n,1);
  bz= zeros(n,1);

% ---------------------------------------------------------------------------
% fast, compact calculation -- but could produce NAN & INF results
  for ii=1:2
    xx= x+pm(ii)*a;
    for jj=1:2
      yy= y+pm(jj)*b;
      for kk=1:2
        zz= z+pm(kk)*c;
        ss= sqrt(xx.^2+yy.^2+zz.^2);
        ty1= log(zz+ss);
        tx1=       xx.*ty1;
        ty1=       yy.*ty1;
        tx1= tx1 + zz.*log(xx+ss);
        ty1= ty1 + zz.*log(yy+ss);
        tx1= tx1 - yy.*atan(xx.*zz./(yy.*ss));
        ty1= ty1 - xx.*atan(yy.*zz./(xx.*ss));
        sgn= 1-2*mod(ii+jj+kk,2); % quick way to get (-1)^(i+j+k)
        bx= bx - sgn*tx1;
        by= by + sgn*ty1;
      end
    end
  end

% ---------------------------------------------------------------------------
% PART II
% check for nan & inf and carefully recalculate
% ---------------------------------------------------------------------------
%
% NOTE: this approximation gives single precision results around eps of edges

  idd= find(isnan(bx+by)|~finite(bx+by));
  if ~isempty(idd)
   nn= length(idd);
   disp(['rect_3d_b.m correcting INF & NAN entries: ',int2str(nn)]) 
   sml= sqrt(2*eps)*max([a(idd)';b(idd)';c(idd)'])';
   bx(idd)= zeros(nn,1);
   by(idd)= zeros(nn,1);
   for ii=1:2
    xx= x(idd)+pm(ii)*a(idd);
    for jj=1:2
      yy= y(idd)+pm(jj)*b(idd);
      for kk=1:2
       zz= z(idd)+pm(kk)*c(idd);
       ss= sqrt(xx.^2+yy.^2+zz.^2);
       tzz= log(zz+ss);

       id= find(yy~=0 & finite(tzz));
       ty1= zeros(nn,1);
       ty1(id)= yy(id).*tzz(id);
       tx3= zeros(nn,1);
       tx3(id)= yy(id).*atan(xx(id).*zz(id)./(yy(id).*ss(id)));

       id= find(xx~=0 & finite(tzz));
       tx1= zeros(nn,1);
       tx1(id)= xx(id).*tzz(id);
       ty3= zeros(nn,1);
       ty3(id)= xx(id).*atan(yy(id).*zz(id)./(xx(id).*ss(id)));

       tzz= log(xx+ss);
       id= find(finite(tzz));
       tx2= zeros(nn,1);
       tx2(id)= zz(id).*tzz(id);

       tzz= log(yy+ss);
       id= find(finite(tzz));
       ty2= zeros(nn,1);
       ty2(id)= zz(id).*tzz(id);

       sgn=   ii+jj+kk;
       sgn=   1-2*mod(sgn,2);

       bx(idd)= bx(idd) - sgn*(tx1+tx2-tx3);
       by(idd)= by(idd) + sgn*(ty1+ty2-ty3);
     end
    end
   end
  end % if

% ---------------------------------------------------------------------------
% convert to b field: B= 4*pi*1e-7*H

  ss= 1.0e-7./(a.*b);
  bx= ss.*bx;
  by= ss.*by;

  return
