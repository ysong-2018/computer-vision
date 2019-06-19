function filt = blurfilter(sigma, n)
%BLURFILTER Creates a Gaussian blur filter.
%   FILT = BLURFILTER(SIGMA, N) creates a Gaussian blur filter with
%   bandwidth parameter SIGMA, and size N x N. If N isn't specified, it
%   defaults to ceil(6*sigma) + 1.
%
%   The output of this function is the filter FILT.

if ~exist('n', 'var')
  n = ceil(6*sigma) + 1;
end

rad = (n - 1)/2;
[xs, ys] = meshgrid(linspace(-rad, rad, n), linspace(rad, -rad, n));

% Your code here -- (make sure you understand what xs and ys represent!)

[a,b] = size(xs);
M = zeros(a,b);

for i = 1:a
    for j = 1:b
        M(i,j)=exp(-(xs(i,j)^2+ys(i,j)^2)/(2*sigma^2));
    end
end

M = M./(sum(sum(M)));

filt = M;







