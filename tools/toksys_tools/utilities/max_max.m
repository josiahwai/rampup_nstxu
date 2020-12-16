  function [mx,i,j]= max_max(array)

% max_max.m function returns the maximum value of an array and the i,j index's
%
% SYNTAX:  [mx,i,j]= max_max(array);
%
% INPUT:
%       array=    input 2-d array
%
% OUTPUT:
%       mx=       scaler of maximum value
%       i=        row of maximum value
%       j=        column of maximum value
%
% Caution: only tested on 2-d arrays (not on vectors or scalars)

% Jim Leuer 4-17-97

   [arr1,ic]= max(array);
   [mx,j]=    max(arr1);
   i=         ic(j);

  return
