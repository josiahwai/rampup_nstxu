function modvalue = modulo(number,base)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  modvalue = modulo(number,base)
%
%  PURPOSE: Calculates remainder when number is divided by base value.
%
%  INPUT:
%	number = number to reduce modulo base
%	base   = base of modular arithmetic
%
%  OUTPUT:
%	modvalue = number reduced modulo base
%
%  RESTRICTIONS:  This works for "all" floating point numbers and for every
%  pair of integers I've tried, but since matlab does only floating point 
%  arithmetic I expect that some combination of integers will give the wrong
%  answer.  Some additional logic needs to be added sometime to address this.

%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	4/3/96
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp = floor(number/base);

modvalue = number - temp*base;


