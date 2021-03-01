  function [bx,by,bz]= filr_3d_b(x,y,z,c,d)

% PURPOSE:
% Generates B-field in carthesian 3-space [X] == [x,y,z] from straight filament
% Filament has finite radius [d] and length [c] in z direction.
% Inside r<d the field is proportional to r rather than 1/r for filament.
% Filament has uniform current density and carries unit current (I=1A) 
% in the +z direction.
% The coordinate system is centered about the geometric center of filament.
% Note: All Units are MKS (meters, Amps, Tesla)
%
% SYNTAX:
%        [bx,by,bz]= filr_3d_b(x,y,z,c,d)
%
% INPUT:
%       [x,y,z] == [X]  vectors of coordinates, one point per row
%       [c] ==          vector  of filament lengths one filament per row
%       [d] ==          vector  of filament radius one filament per row
%
%  Note: If x,y,z & c,d are matrix or vectors they must be same size.
%        Any of the components can be a scalar. 
%
% OUTPUT:
%       [bx,by,bz] == [B] vectors of B-fields at point X
%
% NOTE: use fil_3b_b if [X] is outside filament for faster execution 
%

% Jim Leuer, General Atomics, 1-98

% ---------------
% initialization stuff


  if nargin <= 3
    prog_name= 'filr_3d_b';
    disp(['%ERROR ',prog_name,': Must have at least 4 arguments']);
    eval(['help ',prog_name,]) % print out help from prog_name
    return
  elseif nargin <= 4
    prog_name= 'filr_3d_b';
    disp(['%Caution ',prog_name,': no filament rad given, using fil_3d_b']);
    [bx,by,bz]= fil_3d_b(x,y,z,c);
    return
  end

% ---------------------------------------------------------------------------
% compute magnetics using the least memory and fastest computation

% Constants
  r2= x.^2+y.^2;
  d= d.^2; % caution reuse d to reduce storage
  
% ---------------
% Calculation for points outside filament radius

  id= find(r2 >= d); % find points outside filament

  if ~isempty(id)
    d(id)= r2(id); % this takes care of inside filament
  end

% ---------------
% Calculation (points inside filament radius use d2, outside use r2)

  z= z - 0.5*c; % Caution reuse of z to lower mem requirements
  bx= z./sqrt(r2+z.^2); % 1st term 
  z= z + c;
  bx= bx - z./sqrt(r2+z.^2); % 2nd term
  bx= 1e-7*bx./d;     % Inside= I*mu/(4pi*d^2); outside= I*mu/(4pi*r^2)
  by= -x.*bx;
  bx= y.*bx;
  bz= zeros(size(bx));           % no axial current

% ---------------
% Set any NAN to zero since this occures at r= 0 for which => b= 0

  id=find(isnan(bx+by));
  if ~isempty(id)
    bx(id)= zeros(size(id));
    by(id)= bx(id);
  end

  return

% ------------------------------------------------------------
% TEST Stuff Below:
%  bx= 1*bx./r2/(4*pi);      % fake constant this is H

% test input parameters

%  if exist('fil_3d_b.diary') == 2 
%     ! rm fil_3d_b.diary
%  end
%  format compact
%  diary fil_3d_b.diary
%  nargin= 4;

%  x= 0:.01:1; y= 0:.02:2; z= 0:.02:2; c=ones(size(x)); d=ones(size(x));
%  x= 0:.01:2; y= zeros(size(x)); z= 0:.01:2; c=ones(size(x)); d=ones(size(x));
%  x= 0:.01:2; y= zeros(size(x)); z= zeros(size(x));
%  x= 0:.01:2; y= zeros(size(x)); z= ones(size(x));
%  c=ones(size(x)); d=ones(size(x));
%  [bx,by,bz]= filr_3d_b(x,y,z,c,d);
%  bx= bx*1e+6; by= by*1e+6;
%  b= sqrt(bx.^2+by.^2);
%  [bx1,by1,bz1]= fil_3d_b(x,y,z,c);
%  bx1= bx1*1e+6; by1= by1*1e+6;
%  b1= sqrt(bx1.^2+by1.^2);
%  r= sqrt(x.^2+y.^2);
%  figure(1)
%%  clf
%  plot([r,r],[b,b1])
%  hold on
%  grid on
%  axis([0,2,0,.4])
%
%
%  format short
%  disp(' Test Output of fil_3d_b:')
%  disp('     x         y         z         c')
%  disp([x,y,z,c]);
%
%% ---------------
% print results below:
% 
%  b= sqrt(bx.^2+by.^2);
%  format long
%  disp('           bx                by                  bz              B')
%  disp([bx,by,bz,b]);
%  format short
%
%  diary off
%
