function [zero_x_set,nearest_idx] = zero_crossing(x,y,plot_zeros)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  [zero_x_set,nearest_idx] = zero_crossing(x,y,plot_zeros)
%
%  PURPOSE:  For a function defined by (x,y) data pairs, estimate the value(s) 
%	of x at which the function crosses y=0 (assuming linear interpolation
% 	between data points).
%
%  INPUT:
%	x,y = set of data pairs defining function
%	plot_zeros = (optional) flag to plot zeros on top of function
%
%  OUTPUT:
%	zero_x_set = x values at which y=0
%	nearest_idx = index of point nearest to the zero crossing
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	8/1/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
debug=0;

if(nargin<3)
   plot_zeros=0;
end

idxplus = find(y>0);
idxminus = find(y<0);
zero_idx_set = setdiff([1:length(x)],union(idxplus,idxminus));
zero_x_set = x(zero_idx_set);
nearest_idx = zero_idx_set;

y1 = zeros(size(y));
y1(idxplus) = 1;
y1(idxminus) = -1;
y1(zero_idx_set)= 1;

if debug
   next_figure
   plot(x,y,'x-')
   hold on
   plot(x,y1,'ro-')
   hold off
end

% use diff to find all sign changes

ydiff = diff(y1);
idxplus = find(ydiff>0);	% change from -1 to +1
idxminus = find(ydiff<0);	% change from +1 to -1

if debug
  next_figure
  plot(ydiff,'*')
end

for k=union(idxplus,idxminus)
   mult = abs(y(k)/(y(k)-y(k+1)));
   zero_x = x(k) + mult * (x(k+1)-x(k));
   zero_x_set = [zero_x_set zero_x];
   if(mult<0.5)
      nearest_idx = [nearest_idx k];
   else
      nearest_idx = [nearest_idx k+1];
   end
end

if(plot_zeros) 
   next_figure
   plot(x,y)
   hold on
   grid on
   plot(zero_x_set,zeros(size(zero_x_set)),'x')
   hold off
end
