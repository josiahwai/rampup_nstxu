  function [mn,i,j]= min_min(array)

% min_min.m function returns the minimum value of an array and the i,j index's
%
% SYNTAX:  [mn,i,j]= min_min(array);
%
% INPUT:
%       array=    input 2-d array
%
% OUTPUT:
%       mn=       scaler of maximum value
%       i=        row of maximum value
%       j=        column of maximum value
%
% Caution: only tested on 2-d arrays (not on vectors or scalars)

% Jim Leuer 4-17-97

   [arr1,ic]= min(array);
   [mn,j]=    min(arr1);
   i=         ic(j);

  return
