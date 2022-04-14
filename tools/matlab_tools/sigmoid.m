function y = sigmoid(x, scale, shift, plotit)

if ~exist('scale', 'var'), scale = 1; end
if ~exist('shift', 'var'), shift = 0; end
if ~exist('plotit', 'var'), plotit = 0; end

y = 1 ./ (1 + exp(-scale .* (x-shift)));
y = y(:);

if plotit
  plot(x, y)
end