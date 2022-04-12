function s = nansum(x)

s = sum(x(~isnan(x)));
