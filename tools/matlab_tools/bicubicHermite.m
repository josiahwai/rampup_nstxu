function [z0, dz0dx, dz0dy] = bicubicHermite(xg, yg, zg, x0, y0)
%
% BICUBICHERMITE
%
%   Perform interpolation of 2-D gridded data using bicubic Hermite splines
%
% USAGE: bicubicHermite.m
%
% INPUTS:
%
%   xg......coordinates of grid points in first dimension  (n x 1)
%
%   yg......coordinates of grid points in second dimension (m x 1)
%
%   zg......values of the underlying function at the grid points (m x n)
%
%   x0......coordinate of query point in first dimension
%
%   y0......coordinate of query point in second dimension
%
% OUTPUTS
%
%   z0......interpolated value at the query point
%
%   dz0dx...derivative of interpolant at query point w-r-t first dim
%
%   dz0dy...derivative of interpolant at query point w-r-t second dim
%
% METHOD: This function implements the method of cubic convolution 
%         interpolation as described in the following paper:
%         R. G. Keys. "Cubic Convolution Interpolation for Digital
%         Image Processing." IEEE 1981.
%
% AUTHOR: Patrick J. Vail
%
% DATE: 06/09/2017
%
% MODIFICATION HISTORY:
%   Patrick J. Vail: Original File 06/09/2017
%
%..........................................................................

% Grid point weighting matrix for cubic convolution (Keys, IEEE 1981)

mx = 1/2 * [0 2 0 0; -1 0 1 0; 2 -5 4 -1; -1 3 -3 1];

% Convert xg and yg to column vectors if necessary

if size(xg,1) == 1
    xg = xg';
end
if size(yg,1) == 1
    yg = yg';
end

dx = xg(2) - xg(1);
dy = yg(2) - yg(1);

% Convert zg to (m x n) if necessary

nx = length(xg);
ny = length(yg);

if size(zg,1) == 1 || size(zg,2) == 1
    zg = reshape(zg, ny, nx);
end

% Locate four grid points in each direction around the query point
    
ix = find(xg < x0, 1, 'last');
iy = find(yg < y0, 1, 'last');
   
ii1 = ny*(ix-2) + [iy-1 iy iy+1 iy+2];
ii2 = ny*(ix-1) + [iy-1 iy iy+1 iy+2];
ii3 = ny*(ix+0) + [iy-1 iy iy+1 iy+2];
ii4 = ny*(ix+1) + [iy-1 iy iy+1 iy+2];
     
F = [zg(ii1(1)) zg(ii1(2)) zg(ii1(3)) zg(ii1(4)); ...
     zg(ii2(1)) zg(ii2(2)) zg(ii2(3)) zg(ii2(4)); ...
     zg(ii3(1)) zg(ii3(2)) zg(ii3(3)) zg(ii3(4)); ...
     zg(ii4(1)) zg(ii4(2)) zg(ii4(3)) zg(ii4(4))  ...
];

% Normalize the x and y intervals

tx = (x0 - xg(ix))/dx;
ty = (y0 - yg(iy))/dy;

% Interpolate to find value at the query point

b0 = [1 tx tx^2 tx^3]*mx*F(:,1);
b1 = [1 tx tx^2 tx^3]*mx*F(:,2);
b2 = [1 tx tx^2 tx^3]*mx*F(:,3);
b3 = [1 tx tx^2 tx^3]*mx*F(:,4);

b0_r = ([0 1 2*tx 3*tx^2]/dx)*mx*F(:,1);
b1_r = ([0 1 2*tx 3*tx^2]/dx)*mx*F(:,2);
b2_r = ([0 1 2*tx 3*tx^2]/dx)*mx*F(:,3);
b3_r = ([0 1 2*tx 3*tx^2]/dx)*mx*F(:,4);

z0 = [1 ty ty^2 ty^3]*mx*[b0 b1 b2 b3]';

dz0dx = [1 ty ty^2 ty^3]*mx*[b0_r b1_r b2_r b3_r]';
dz0dy = ([0 1 2*ty 3*ty^2]/dy)*mx*[b0 b1 b2 b3]';

end
