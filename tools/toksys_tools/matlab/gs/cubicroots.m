function [x1, x2, x3, alfa, beta] = cubicroots(a, b, c, d, n)
%
%  USAGE:   [x1, x2, x3] = cubicroots(a, b, c, d)
%
%  PURPOSE: Find the roots to the cubic equation:
%           a*x^3 + b*x^2 + c*x + d = 0
%           Can also be used with a = 0 and/or b = 0
%
%  INPUTS:  a, b, c, d - may be 1 arbitrary size in a mix with scalars
%
%  OUTPUTS: The roots such that:
%           (x-x1)*(x-x2)*(x-x3) = x^3 + b/a*x^2 + c/a*x + d/a
%           If a == 0 then x3 is nan since polynomial is 2:nd degree
%           If also b == 0 then x2 will also be nan
%           If also c == 0 then x1 will also be nan
%
%  NOTES:   Faster than the function roots, in particular when several
%           equations solved at once. Also more accurate than roots
%           Solutions are in the same order as from function roots
%           If a ~= 0, x1 is always real
%           cubicroots(a,b,c,d,n) iterates n times to reduce tiny errors

%           annoying small error with this call
%           [x1, x2, x3] = cubicroots(a,b,0,0)
%           solution must work for complex a, b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
%
%  WRITTEN BY:  Anders Welander  ON	3/17/14
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
  n = 0;
end

alfa = (c./a - b.^2./(3*a.^2))/3;
beta = 2*b.^3./(27*a.^3) - c.*b./(3*a.^2) + d./a;

e = 4*alfa.^3./beta.^2;
s = sqrt(1+e)-1;
i = abs(s) < 1e-4;
s(i) = e(i)/2-e(i).^2/8+e(i).^3/16; % -e(i).^4*5/128;
z = s.*beta/2;

p = z.^(1/3);
q = -alfa./p;

i = alfa == 0;
q(i) = beta(i).^(1/3);

y = p + q;

y(beta == 0) = 0;

x1 = y-b/3./a;
sq = sqrt((b./a+x1).^2/4+d./a./x1);
x2 = -(b./a+x1)/2 + sq;
x3 = -b./a-x1-x2;

for k = 1:n
  d1 = 3*a.*x1.^2+2*b.*x1+c;
  d2 = 3*a.*x2.^2+2*b.*x2+c;
  d3 = 3*a.*x3.^2+2*b.*x3+c;

  e1 = a.*x1.^3+b.*x1.^2+c.*x1+d;
  e2 = a.*x2.^3+b.*x2.^2+c.*x2+d;
  e3 = a.*x3.^3+b.*x3.^2+c.*x3+d;

  x1(d1~=0) = x1(d1~=0) - e1(d1~=0)./d1(d1~=0);
  x2(d2~=0) = x2(d2~=0) - e2(d2~=0)./d2(d2~=0);
  x3(d3~=0) = x3(d3~=0) - e3(d3~=0)./d3(d3~=0);
end

i = abs(imag(x1)) < 1e-9*abs(real(x1));
x1(i) = real(x1(i));

i = abs(imag(x2)) < 1e-9*abs(real(x2));
x2(i) = real(x2(i));

i = abs(imag(x3)) < 1e-9*abs(real(x3));
x3(i) = real(x3(i));

i = logical(ones(size(x1)));

% The case where a = 0   =>   x^2 + c/b*x + d/b = 0  

i = i & a == 0;
x3(i) = nan;
s(i) = sqrt(1-4*b(i).*d(i)./c(i).^2)-1;
x1(i) = s(i).*c(i)./b(i)/2;
x2(i) = (-2-s(i)).*c(i)./b(i)/2;


% The case where a = 0 & b = 0   =>   c*x + d = 0  

i = i & b == 0;
x2(i) = nan;
x1(i) = -d(i)./c(i);


% The case where a = 0 & b = 0 & c = 0   =>   d = 0  

i = i & c == 0;
x1(i) = nan;


% The case where c = 0 & d = 0  =>   a*x^3 + b*x^2  = 0  

i = c == 0 & d == 0 & a ~= 0;
x1(i) = -b(i)./a(i);
x2(i) = 0;
x3(i) = 0;
e(i) = angle(x1(i));
i = i & (e > pi/3 | e < -pi/3);
x2(i) = x1(i);
x1(i) = 0;
i = i & (e > 2*pi/3 | e < -2*pi/3);
x3(i) = x2(i);
x2(i) = 0;



% Notes from web site:
% http://au.answers.yahoo.com/question/index?qid=20070509175706AA5KWEa

% Solving cubic equations - Cardano's method:

% A cubic equation has the form ax³ + bx² + cx + d = 0. 
% Now a must be nonzero, or else this is in fact a quadratic equation
% and can be solved using the methods for that. 
% So divide by a, to get:

% x³ + b/a x² + c/a x + d/a = 0

% Now, we make a change of variables. 
% Let y=x+b/(3a). Then x=y-b/(3a), so we have:

% (y-b/(3a))³ + b/a (y-b/(3a))² + c/a (y-b/(3a)) + d/a = 0

% Expanding:

% y³ - b/a y² + b²/(3a²) y - b³/(27a³) + b/a y² - 2b²/(3a²) y + b³/(9a³) + c/a y - cb/(3a²) + d/a

% Simplifying:

% y³ + (c/a - b²/(3a²)) y + (2b³/(27a³) - cb/(3a²) + d/a)

% Let's look at what we have done here:
% we have taken the general cubic equation and,
% through our change of variables,
% transformed it into a cubic equation with no quadratic coefficient.
% This type of cubic equation is known as a depressed cubic.
% If we can find solutions for the depressed cubic,
% we can then reverse the substitution we just made and
% thereby solve the general cubic.
% So let us solve the depressed cubic.
% But first, we introduce some new variables
% to keep the notation from becoming unwieldy.
% Let alfa=c/a - b²/(3a²) and beta=2b³/(27a³) - cb/(3a²) + d/a.

% So we have:

% y³ + alfa*y + beta = 0

% Now, if alfa=0, then we have y³+beta=0, which means that y³=-beta,
% so this equation has solutions y=(-beta)^(1/3)*[cos(v3) + 1i*sin(v3)],
% where v3 = 2*pi/3*[0 1 2];
% Similarly, if beta=0, then we can easily factor this cubic as y(y²+alfa)
% and it will have solutions y=0, y=sqrt(-alfa), and y=-?(-?).
% So in either case if we are fortunate enough to have either
% alfa or beta equal to zero the cubic can be easily solved,
% so we can direct our attention to the case where both are nonzero.

% In order to solve the cubic in this case,
% we make the assumption that the solution,
% whatever it is, can be written as the sum of two other numbers p and q
% (this will actually help us). Now, making the substitution y=p+q:

% (p+q)³ + alfa(p+q) + beta = 0

% We partially expand the product:

% p³ + 3pq(p+q) + q³ + alfa(p+q) + beta = 0

% Regrouping:

% p³ + q³ + beta + (3pq+alfa)(p+q) = 0

% Now, clearly if we can choose p and q such that
% p³ + q³ + beta = 0 and 3pq+alfa = 0,
% then p+q will be a solution to the original equation.
% So this gives us the system of two equations:

% 3pq+? = 0
% p³ + q³ + ? = 0

% Solve the first equation for q:

% q=-alfa/(3p)

% Substitute into the second equation:

% p³ - alfa³/(27p³) + beta = 0

% Multiply by p³:

% p^6 + beta*p³ - alfa³/27 = 0

% Now, at this point you may wonder what we're doing:
% we started with an equation of degree 3,
% and transformed it into an equation of degree six!
% But note that there are only two terms with a variable,
% and one is the square of the other.
% So now we can make the substitution z=p³ to obtain:

% z² + beta*z - alfa³/27 = 0

% And this is a quadratic equation.
% We know how to solve quadratic equations, so let's do that:

% z = (-beta ± sqrt(beta²+4alfa³/27))/2

% Now, p³=z, so the possible solutions for p are
% p=z^(1/3)*[cos(v3) + 1i*sin(v3)], p=(i?3-1)/2
% Now we must find the values of q.
% The naive solution would be to substitute these values for p
% into the equation 3pq+alfa = 0 and solve, which produces a correct
% (but very ugly) solution of the cubic.
% Fortunately, there is a much nicer method:
% notice that the original system of equations we set up
% were symmetric polynomials in p and q -- in other words,
% any solution for p can also be a solution for q, and vice versa.
% There are six solutions for p, each solution having a partner
% which is the corresponding value of q for that p.
% So we actually already have the corresponding values of q,
% and all we need do is figure out which solution gets paired with which.
% If we designate the solutions to the quadratic equation
% z1 = (-beta + sqrt(beta²+4alfa³/27))/2 and
% z2 = (-beta - sqrt(beta²+4alfa³/27))/2,

% the z1, z2, and alfa are all real, 
% and that the real-valued cube root is taken in all cases,
% then the proper pairings are:

% p1 = cuberoot(z1), q1=cuberoot(z2)
% p2 = (i*sqrt(3)-1)/2*cuberoot(z1) , q2 = (-i*sqrt(3)-1)/2*cuberoot(z2)
% p3 = (-i?3-1)/2 ?zv(1), qv(3) = (i?3-1)/2 ?zv(2)

% Finally, now that we have the values of p and q,
% the solutions of the depressed cubic are 
% yv(1)=pv(1)+qv(1), yv(2)=pv(2)+qv(2), and yv(3)=pv(3)+qv(3),
% and the solutions of the original cubic equation are
% xv(1)=yv(1)-b/(3a), xv(2)=yv(2)-b/(3a), and xv(3)=yv(3)-b/(3a). 


% Notes by Anders Welander

% x1 is always real but what about x2 and x3?
% The original equation is:
% (x-x1)*(x-x2)*(x-x3) = x^3 + b/a*x^2 + c/a*x + d/a
% We find x1 by making sure we have the real root, and we know that:
% x1+x2+x3 = -b/a
% x1*x2*x3 = -d/a
% Therefore
% x3 = -b/a-x1-x2
% x1*x2*(b/a+x1+x2)-d/a = 0
% x1*x2^2 + x1*(b/a+x1)*x2 - d/a = 0
% x2^2 + (b/a+x1)*x2 -d/a/x1 = 0
% x2 = -(b/a+x1)/2 - sqrt((b/a+x1)^2/4+d/a/x1)
% x3 = -b/a-x1-x2
% This guarantees that the imaginary part is always exactly zero
% when the roots are in fact real, and it executes faster this way
% As an extra check, we should have:
% x1*x2+x1*x3+x2*x3-c/a = 0
% This code now deals with numerical precision issues that aren't
% treated on the web site.
% It was tested extensively with around 1e8 calculations and has equal
% or better accuracy than roots and the speed appears to be very roughly
% 25*log10(n) higher with this code compared to roots. This is partly
% because roots can only do one solution at a time.
% For 1e0 calculations, cubicroots is   1 times faster.
% For 1e1 calculations, cubicroots is   4 times faster.
% For 1e2 calculations, cubicroots is  22 times faster.
% For 1e3 calculations, cubicroots is  42 times faster.
% For 1e4 calculations, cubicroots is  77 times faster.
% For 1e5 calculations, cubicroots is 107 times faster.
% For 1e6 calculations, cubicroots is 158 times faster.
% In all these cases the roots calculation took 900 us on average
% but the cubicroots calculations get faster the more that are done at once
% For testing...
% a.*x1.^3 + b.*x1.^2 + c.*x1 + d
% a.*x2.^3 + b.*x2.^2 + c.*x2 + d
% a.*x3.^3 + b.*x3.^2 + c.*x3 + d
% min(min(min(a.*x1.^3 + b.*x1.^2 + c.*x1 + d)))
% min(min(min(a.*x2.^3 + b.*x2.^2 + c.*x2 + d)))
% min(min(min(a.*x3.^3 + b.*x3.^2 + c.*x3 + d)))
% max(max(max(a.*x1.^3 + b.*x1.^2 + c.*x1 + d)))
% max(max(max(a.*x2.^3 + b.*x2.^2 + c.*x2 + d)))
% max(max(max(a.*x3.^3 + b.*x3.^2 + c.*x3 + d)))
