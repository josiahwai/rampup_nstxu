function error = error_data2(x,data)

% Calculate error between x1 + x2*k and data.

% This routine assumes that data(:,1) is the independent variable
% and data(:,2) the associated dependent variable.

[len width] = size(data);
estimate = x(1) + x(2)*data(:,1);
error = estimate-data(:,2);

