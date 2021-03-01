 %gridit:
%  Script to put grid lines at selected values of x,y.
% Usage:
%   xgrid=linspace(1.7,.01,1.74); ygrid=linspace(0,.5,6); gridit;
% Inputs: define these variables (as desired):
%   xgrid = vector of grid lines corr. to x-axis points
%   ygrid = vector of grid lines corr. to y-axis points

   atmp=axis;
   hold on

% Do x-axis grid lines:
   if exist('xgrid')&~isempty(xgrid)
    for ii=1:length(xgrid)
     plot([xgrid(ii) xgrid(ii)], atmp(3:4),'w:')
    end
   end

% Do y-axis grid lines:
   if exist('ygrid')&~isempty(ygrid)
    for ii=1:length(ygrid)
     plot(atmp(1:2),[ygrid(ii) ygrid(ii)],'w:')
    end
   end

                               
