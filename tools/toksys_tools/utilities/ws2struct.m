function WStruct=ws2struct()
%   @(#)ws2struct.m       1.1       14/02/14
%  USAGE: ws2struct
%
%  PURPOSE: This function allows to save all the variables from the 'caller' workspace into a struct array
%
%  INPUTS:
%
%  OUTPUTS: 
%
%  RESTRICTIONS: 
%
%  METHOD:  
%Example:
% a='LALALA'
% b=[1:12:258]
% c={'cell1', 'cell2', 'cell3'}
% d=768
% e=true(3)
% theworkspace=ws2struct();
% theworkspace = 
% 
%     a: 'LALALA'
%     b: [1x22 double]
%     c: {'cell1'  'cell2'  'cell3'}
%     d: 768
%     e: [3x3 logical]
%
%  VERSION 
%
%  WRITTEN BY:  Matthew J. Lanctot on March 19 2013
%
%  MODIFICATION HISTORY:
% 	2013-03-19	Created
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

WSVARS = evalin('caller', 'who');
for wscon=1:size(WSVARS,1)
    thisvar=evalin('caller', WSVARS{wscon});
    THEWORKSPACE.(WSVARS{wscon})=thisvar;
end

WStruct=THEWORKSPACE;


