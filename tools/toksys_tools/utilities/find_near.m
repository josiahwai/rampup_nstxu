function idx = find_near(vec,t0);
% function idx=find_near(vec,t0)
%   Function to find vec index with vec value nearest to t0 (eg find index 
%   nearest to time t0 in a time vector vec). t0 can be a vector as well. Then 
%   returns vector of nearest value indices corr to t0 values.

if length(t0)==1
  idx = find(abs(vec-t0)==min(abs(vec-t0)));
else
 idx=zeros(length(t0),1);
 for ii=1:length(t0)
  tmp = find(abs(vec-t0(ii))==min(abs(vec-t0(ii))));
  idx(ii) = tmp(1);
 end
end


