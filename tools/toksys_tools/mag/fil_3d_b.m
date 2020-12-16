  function [bx,by,bz]= fil_3d_b(x,y,z,c)

% PURPOSE:
% Generates B-field in carthesian 3-space [X] == [x,y,z] 
% produced by a straight filament with length [c] in the z direction
% carrying Unit current (I=1A) in the +z direction
% The coordinate system is centered about the geometric center of filament.
% Note: All Units are MKS (meters, Amps, Tesla)
%
% SYNTAX:
%        [bx,by,bz]= fil_3d_b(x,y,z,c)
%
% INPUT:
%       [x,y,z] == [X]  coordinates of field point for calculation of B
%       [c] ==          length of filament (scalar or same size as X)
%
% Note: x,y,z or c can be a scalar or matrix but all matricies must be same
%       size. 
%       Any of the components can be a scalar. 
%
% Caution: if x and y are both zero answer is NAN (use fil_3d_b_rad)
%
% OUTPUT:
%       [bx,by,bz] == [B] B-field at point X
%
% NOTE: see filr_3d_b for filaments with finite radius computation inside rad.

% Jim Leuer, General Atomics, 12-97

% ---------------
% initialization stuff

  if nargin <= 3
    prog_name= 'fil_3d_b';
    disp(['%ERROR ',prog_name,': Must have at least 4 arguments']);
    eval(['help ',prog_name,]) % print out help from prog_name
    return
  end

% ---------------------------------------------------------------------------
% compute magnetics using the least memory and fastest computation

% ---------------
% start computation

% Constants
  r2= x.^2+y.^2;

% ---------------
  z= z - 0.5*c; % Caution reuse of z to lower mem requirements
  bx= z./sqrt(r2+z.^2); % 1st term 
  z= z + c;
  bx= bx - z./sqrt(r2+z.^2); % 2nd term
  bx= 1e-7*bx./r2;         % constant I*mu/(4pi)= 1*4pi*1e-7/(4pi)
  by= -x.*bx;
  bx= y.*bx;
  bz= zeros(size(bx));           % no axial current

% ---------------
% 
  return
% ---------------
%  bx= 1*bx./r2/(4*pi);      % fake constant this is H

% test input parameters

%  if exist('fil_3d_b.diary') == 2 
%     ! rm fil_3d_b.diary
%  end
%  format compact
%  diary fil_3d_b.diary

%  nargin= 4;
%  x= [1;-1];  y= [2;2];  z= [3;3]; c=[1;1];
%  x= x(1); y= y(1); z= z(1); c= c(1); 

%  format short
%  disp(' Test Output of fil_3d_b:')
%  disp('     x         y         z         c')
%  disp([x,y,z,c]);

% ---------------
% print results below:
 
%  b= sqrt(bx.^2+by.^2);
%  format long
%  disp('           bx                by                  bz              B')
%  disp([bx,by,bz,b]);
%  format short

%  diary off
