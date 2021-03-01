function [xi,yi,found] = intersect_bdry(x,y,xstart,ystart,xend,yend,bounded)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  [xi,yi,found] = intersect_bdry(x,y,xstart,ystart,xend,yend,bounded)
%
%  PURPOSE:  Given a set of connected line segments (such as a closed
%	boundary) and the end points of an independent line segment, 
%	find the intersection of the independent segment
%	and the boundary.  The first intersection found will be returned
%	so for reliability the segment should intersect the boundary
%	in at most one place.
%
%  INPUT:
%    x,y           = (x,y) of boundary (vectors)
%    xstart,ystart = (x,y) of start of independent segment
%    xend,yend     = (x,y) of end of independent segment
%    bounded = set to 1 to require the intersection to be between endpoints. 
%	       set to 0 to allow intersection anywhere along the infinite lines 
%	 	defined by extending the independent line segment and a segment
%		in the boundary. (In this case, the intersection will be found
%		for the first segment in the bdry which is not parallel to the
%		independent segment.
%
%  OUTPUT:
%    xi,yi = (x,y) coordinate of intersection
%    found = 0, if no intersection found.
%	   = index of second point on intersecting bdry segment, otherwise.
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	9/9/97 based on J.Ferrons intersect.pro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equation of line segment.
%
   if(xstart == xend) 	% line segment has infinite slope.
      ptg = 0;
   else 		% line segment has finite slope.
      mg = (ystart - yend)/(xstart - xend);    %slope of line segment.
      bg = ystart - mg*xstart;                 %intercept of line segment.
      ptg = 1;
   end
%
%----------------------------------------------------------------------
%Find the boundary segment that intersects the line segment.
%The line segment should be short enough that it intersects the bdy
%at only one point.
%
   nbdry = length(x);		% number of boundary points
   found = 0;
   xi = 0.0;
   yi = 0.0;
   j = 1;
%
   while( (found == 0) & (j < nbdry-1))
      j = j+1;
%
%----------------------------------------------------------------------
%Intersection point between this boundary segment and the line segment.
%Unless the segments are parallel, they intersect somewhere.
%
      pt = 0;          	% flag to indicate if an intersection was found.
      if(x(j-1) == x(j)) 	% bdy segment has infinite slope.
%
         if(ptg ~= 0)        	% line segment slope not infinite also
            xi = x(j);
            yi = mg*xi + bg;
            pt = 1;                    %an intersection was found.
         end
%
      else 			% bdy segment has finite slope.
%
         mb = (y(j-1) - y(j))/(x(j-1) - x(j));  %slope of bdry line segment
         bb = y(j) - mb*x(j);                   %intercept of bdry line segment
%
         if(ptg == 0) 		% line segment has infinite slope.
            xi = xstart;
            yi = mb*xi + bb;
            pt = 1;                           %an intersection was found.
%
         elseif(mg ~= mb)
%both segments have finite slope and the slopes are not equal.
            xi = (bb - bg)/(mg - mb);        %x where two segments intersect
            yi = mb * xi + bb;               %y where two segments intersect
            pt = 1;
         end
      end
%
%----------------------------------------------------------------------
% If the segments intersect, do they intersect on the portion of the line
% between the defined end points.
%
      if(pt == 1) 
%
         if(bounded == 1) 
%
% do the two lines intersect on the boundary segment?
%
            db = sqrt((x(j)-x(j-1))^2 + (y(j)-y(j-1))^2); % length of segment
            d1 = sqrt((x(j)-xi)^2 + (y(j)-yi)^2);    % distance of intersection 
					 	     % from end point 2.
            d2 = sqrt((x(j-1)-xi)^2 + (y(j-1)-yi)^2);% distance of intersection
					 	     % from end point 1.
%
% Do the two lines intersect on the line segment?
%
            dc = sqrt((xstart-xend)^2 + (ystart-yend)^2);
            dc1 = sqrt((xstart-xi)^2 + (ystart-yi)^2);
            dc2 = sqrt((xend-xi)^2 + (yend-yi)^2);
%
%If yes both times, we are done.
%
            if( d1 <= db & d2 <= db & dc1 <= dc & dc2 <= dc) 
	       found = j;
            end
         else
            found = j;
         end
      end
   end
