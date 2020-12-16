function [num,den]=bld_filter(type,corner)
 %
% SYNTAX: [num,den]=bld_filter(type,corner)
%
% PURPOSE:  Build a 4-pole Bessel or Chebyshev filter.  
% (From FAX from Frequency Devices 1/21/94.)
%
% INPUTS:
%  type = one of:
%		'Bessel' = bessel filter
%		'Chebyshev2' = .2dB Chebychev filter
%		'Chebyshev5' = .5dB Chebychev filter
%  corner = corner frequency for filter (Hz)
%
% OUTPUTS:
%   num, den = numerator and denominator polys in transfer function 
%		describing filter

% WRITTEN BY:  Mike Walker	ON	1/21/94
%
%   			Copyright 1993 General Atomics.  
% Unpublished - rights reserved under the Copyright Laws of the United States.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @(#)bld_filter.m	1.2 06/08/11

if strcmp(type,'Bessel')
   wp1 = 2*pi*(1.430172)*corner;
   wp2 = 2*pi*(1.603358)*corner;
   q1 = .521935;
   q2 = .805538;

%elseif type(1:10) == 'Chebyshev2'
elseif strcmp(type,'Chebyshev2')
   wp1 = 2*pi*(.701109)*corner;
   wp2 = 2*pi*(1.09483)*corner;
   q1 = .645897;
   q2 = 2.435013;

%elseif type(1:10) == 'Chebyshev5'
elseif strcmp(type,'Chebyshev5')
   wp1 = 2*pi*(.597002)*corner;
   wp2 = 2*pi*(1.031270)*corner;
   q1 = .705110;
   q2 = 2.940554;
   
else
   fprintf(['type ' type ' is not available']);
   return
end;      
   
num = wp1^2 * wp2^2;
den = [1 (wp1/q1 + wp2/q2) (wp1^2 + wp2^2 + wp1*wp2/(q1*q2)) ...
			        (wp1*wp2^2/q1 + wp2*wp1^2/q2) wp1^2*wp2^2];



